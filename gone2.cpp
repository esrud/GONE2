/*
 * GONE2 — Genetic Optimization for Ne Estimation
 *
 * Authors: Enrique Santiago, Carlos Köpke
 */

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "lib/population.hpp"
#include "lib/ped_map.hpp"
#include "lib/tped.hpp"
#include "lib/vcf.hpp"
#include "lib/sample.hpp"
#include "lib/libgone.hpp"
#include "lib/ne_ref.hpp"

#ifdef GONE_NCURSES_TUI
#include "lib/ncurses_tui.hpp"
#endif

static std::string g_progressTmpFile;
static std::terminate_handler g_prevTerminate = nullptr;
// Delete the progress .tmp file if its path was recorded; registered as the
// std::atexit handler.
static void RemoveProgressTmpFile() {
  if (!g_progressTmpFile.empty()) std::remove(g_progressTmpFile.c_str());
}
// Write a possibly multi-line failure reason to the _STATS stream, repeating
// the comment prefix on every line so the file stays fully commented.
// Params: out — destination stream; prefix — comment prefix to repeat;
//   text — the reason, lines separated by '\n'.
static void WriteReasonLines(std::ostream& out, const char* prefix,
                             const std::string& text) {
  size_t start = 0;
  while (start <= text.size()) {
    const size_t brk = text.find('\n', start);
    const size_t end = (brk == std::string::npos) ? text.size() : brk;
    out << prefix << text.substr(start, end - start) << "\n";
    if (brk == std::string::npos) break;
    start = brk + 1;
  }
}

// std::terminate handler: remove the temp file, chain to the previous
// terminate handler (preserving its diagnostic), then abort.
static void TerminateCleanup() {
  RemoveProgressTmpFile();
  if (g_prevTerminate) g_prevTerminate();
  std::abort();
}

// -x -T resampled mix inference. Instead of one Fst pass on a single
// SNP subset, run the mix inference on up to kMaxRounds different
// random SNP subsets (size -s, seeded from -S), accumulate each
// converged round's per-bin reconstructed δ²pob, then:
//   - geo-mean δ²pob per recombination bin (over rounds where it was
//     positive — includes the high-δ²/would-be-negative-Ne rounds, so
//     the average is not survivor-biased; only over-corrected δ²≤0
//     rounds are dropped),
//   - mean per-bin SNP-pair count over those rounds (not the sum: the bin
//     compaction thresholds downstream are absolute pair counts),
//   - median Fst / ps / m across converged rounds,
//   - leave the aggregated δ²pob in sampleInfo->d2[] for the GA below.
// Stops early once kTargetConvergedRounds rounds have converged (so each
// bin has up to that many d²pob samples to geo-mean), else after kMaxRounds.
// When neither -s nor -i subsamples the data every round would see the same
// SNPs and individuals — only in a different order — so a single pass is run
// and the geo-mean/median reduce to that pass's own values.
// Params: popInfo — genotype source; popInfoMix — metapop scratch state;
//   sampleInfo — receives the aggregated per-bin δ²pob, raw d², and
//   structure transform; params — run configuration. If no round converges,
//   twoPassGA and mixStructureEstimated are cleared and the reason recorded
//   by the structure search (params->mixFailReason) is reported to stderr,
//   the progress file and the _STATS file.
static void RunResampledMix(PopulationInfo* popInfo,
                            PopulationInfoMix* popInfoMix,
                            SampleInfo* sampleInfo, AppParams* params) {
  int wantIndi = params->numSample > 0 ? params->numSample : popInfo->numIndi;
  if (popInfo->haplotype == 2 && params->numSample > 0) wantIndi *= 2;
  const bool singlePass =
      (params->numSNPs <= 0 || params->numSNPs >= popInfo->numLoci) &&
      wantIndi >= popInfo->numIndi;
  const int kMaxRounds            = singlePass ? 1 : 20;
  const int kTargetConvergedRounds = singlePass ? 1 : 10;

  std::vector<int>      cntD2(MAXBINS, 0);
  std::vector<long int> sumNxc(MAXBINS, 0);
  std::vector<char>     binUsable(MAXBINS, 0);
  std::vector<double>   xcGrid(MAXBINS, 0.0);
  std::vector<double>   FstList, mList;
  std::vector<int>      psList;
  std::vector<double>   logSumRaw(MAXBINS, 0.0);
  std::vector<double>   logSumScale(MAXBINS, 0.0);
  std::vector<double>   sumOffset(MAXBINS, 0.0);
  int    binMaxSeen   = 0;
  int    usableCount  = 0;
  int    converged    = 0;

  for (int round = 0;
       round < kMaxRounds && converged < kTargetConvergedRounds; ++round) {
    const double base = std::max(
        static_cast<double>(converged) / kTargetConvergedRounds,
        static_cast<double>(round) / kMaxRounds);
    const double nextFloor = std::max(
        static_cast<double>(converged) / kTargetConvergedRounds,
        static_cast<double>(round + 1) / kMaxRounds);
    params->progress.SetTaskWindow(base, nextFloor - base);
    params->progress.SetStatusDetail(
        "Resampled -x: round " + std::to_string(round + 1) + ", " +
        std::to_string(converged) + " / " +
        std::to_string(kTargetConvergedRounds) + " converged");

    std::mt19937 gr(static_cast<unsigned>(params->semilla) + round + 1);
    for (int j = 0; j < popInfo->numLoci; ++j) sampleInfo->shuffledLoci[j] = j;
    std::shuffle(sampleInfo->shuffledLoci,
                 &(sampleInfo->shuffledLoci[popInfo->numLoci]), gr);
    if (popInfo->haplotype != 2) {
      for (int j = 0; j < popInfo->numIndi; ++j) sampleInfo->shuffledIndi[j] = j;
      std::shuffle(sampleInfo->shuffledIndi,
                   &(sampleInfo->shuffledIndi[popInfo->numIndi]), gr);
    } else {
      const int ndips = popInfo->numIndi / 2;
      for (int j = 0; j < ndips; ++j) sampleInfo->shuffledIndi[j] = j;
      std::shuffle(sampleInfo->shuffledIndi,
                   &(sampleInfo->shuffledIndi[ndips]), gr);
      for (int j = ndips - 1; j >= 0; --j) {
        sampleInfo->shuffledIndi[j * 2]     = sampleInfo->shuffledIndi[j] * 2;
        sampleInfo->shuffledIndi[j * 2 + 1] = sampleInfo->shuffledIndi[j * 2] + 1;
      }
    }

    sampleInfo->numLoci =
        params->numSNPs > 0 ? params->numSNPs : popInfo->numLoci;
    if (sampleInfo->numLoci > popInfo->numLoci)
      sampleInfo->numLoci = popInfo->numLoci;
    if (params->numSample > 0) {
      sampleInfo->sampleSizeIndi = params->numSample;
      if (popInfo->haplotype == 2) sampleInfo->sampleSizeIndi *= 2;
    } else {
      sampleInfo->sampleSizeIndi = popInfo->numIndi;
    }
    if (sampleInfo->sampleSizeIndi > popInfo->numIndi)
      sampleInfo->sampleSizeIndi = popInfo->numIndi;

    sampleInfo->numSegLoci      = 0;
    popInfo->avgNumIndiAnalyzed = 0;
    CalculateFrequencies(popInfo, sampleInfo);
    CalculateF(sampleInfo, popInfo);
    if (params->haplotype == 0 || params->haplotype == 2)
      CalculateHeterozigosity(popInfo, sampleInfo);
    popInfo->hetEspAll = sampleInfo->hetEspAll;

    CalculateD2Parallel(popInfo, sampleInfo, params);

    params->mixRoundProbe = true;
    params->mixConverged  = true;
    CalculateD2ParallelFst(popInfoMix, popInfo, sampleInfo, params);
    params->mixRoundProbe = false;
    if (params->mixFatal) break;
    if (!params->mixConverged) continue;
    ++converged;

    FstList.push_back(popInfo->Fst);
    mList.push_back(popInfo->m);
    psList.push_back(popInfo->ps);
    if (sampleInfo->binMax > binMaxSeen) binMaxSeen = sampleInfo->binMax;

    for (int b = 0; b < sampleInfo->binMax; ++b) {
      if (sampleInfo->nxc[b] <= 0) continue;
      const double d2pob = sampleInfo->mixBinD2pob[b];
      const double raw   = sampleInfo->d2[b];
      const double sc    = sampleInfo->structScale[b];
      const double of    = sampleInfo->structOffset[b];
      if (!(d2pob > 0.0) || !std::isfinite(d2pob)) continue;
      if (!(raw   > 0.0) || !std::isfinite(raw))   continue;
      if (!(sc    > 0.0) || !std::isfinite(sc))    continue;
      if (!std::isfinite(of)) continue;
      logSumRaw[b]   += std::log(raw);
      logSumScale[b] += std::log(sc);
      sumOffset[b]   += of;
      ++cntD2[b];
      sumNxc[b] += sampleInfo->nxc[b];
      xcGrid[b]  = sampleInfo->xc[b];
      const double c   = sampleInfo->xc[b];
      const double c2  = c * c;
      const double c12 = (1.0 - c) * (1.0 - c);
      const double Ne  = (1 + c2 - 2.2 * d2pob * c12) /
                         (2 * d2pob * (1 - c12));
      if (Ne > 0 && !binUsable[b]) { binUsable[b] = 1; ++usableCount; }
    }
  }

  if (converged == 0) {
    params->twoPassGA = false;
    params->mixStructureEstimated = false;
    std::string reason = params->mixFailReason.empty()
        ? std::string("-x analysis could not converge. No Ne estimate is"
                      " produced.")
        : params->mixFailReason;
    const std::string warn =
        (params->mixFatal
             ? std::string("WARNING: ")
             : singlePass
                 ? std::string("WARNING: the -x structure search did not"
                               " converge. ")
                 : "WARNING: none of the " + std::to_string(kMaxRounds) +
                       " resampled -x rounds converged. ") +
        reason + " Only the _STATS file is written.";
    if (!params->quiet) {
      std::cerr << " " << warn << "\n";
    }
    params->progress.SetWarning(warn);
    params->progress.SetStatusDetail("Mix analysis did not converge");
    popInfo->ps  = 0;
    popInfo->Fst = 0;
    popInfo->m   = 0;
    return;
  }

  auto medianD = [](std::vector<double> v) -> double {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
  };
  auto medianI = [](std::vector<int> v) -> int {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  const double Fst  = medianD(FstList);
  const double mMig = medianD(mList);
  const int    ps   = medianI(psList);
  popInfo->Fst = Fst;
  popInfo->m   = mMig;
  popInfo->ps  = ps;
  params->mixStructureEstimated = true;

  sampleInfo->binMax = binMaxSeen;
  for (int b = 0; b < binMaxSeen; ++b) {
    if (cntD2[b] > 0) {
      sampleInfo->d2[b]           = std::exp(logSumRaw[b]   / cntD2[b]);
      sampleInfo->structScale[b]  = std::exp(logSumScale[b] / cntD2[b]);
      sampleInfo->structOffset[b] = sumOffset[b] / cntD2[b];
      sampleInfo->nxc[b]          = sumNxc[b] / cntD2[b];
      sampleInfo->xc[b]           = xcGrid[b];
    } else {
      sampleInfo->d2[b]           = 0.0;
      sampleInfo->nxc[b]          = 0;
      sampleInfo->structScale[b]  = 1.0;
      sampleInfo->structOffset[b] = 0.0;
    }
  }
  sampleInfo->structuredMix = true;

  if (!params->quiet) {
    if (singlePass) {
      std::cerr << " Single-pass -x (all SNPs and individuals analysed): "
                << usableCount << " usable bin(s).\n";
    } else {
      std::cerr << " Resampled -x: " << converged << " converged round(s), "
                << usableCount << " usable bin(s).\n";
    }
  }
}

// Program entry point: parse the command line, read the input, compute the
// observed d², then run the estimation — the single-population GA, or for
// -x the Stage A/B structured search plus the optional -x -T two-pass GA —
// and write the _GONE2_Ne / _GONE2_d2 / _GONE2_STATS (and -x _mix)
// outputs. Sets up the progress temp-file cleanup handlers.
// Params: argc/argv — command-line arguments.
// Returns: 0 on success; EXIT_FAILURE if the analysis was abandoned (too few
//   recombination bins), in which case only _GONE2_STATS is written, carrying
//   the reason; fatal input errors call exit(EXIT_FAILURE).
int main(int argc, char * argv[]) {
  AppParams params;
  HandleInput(argc, argv, &params);

  std::mt19937 g(params.semilla);
  if (params.mix && params.twoPassGA) {
    params.progress.ConfigurePhases({
        "Reading input files",
        "Preparing sample",
        "Calculating d2",
        "Inferring metapop structure",
        "Fitting N_T",
        "Writing outputs"});
  } else if (params.mix) {
    params.progress.ConfigurePhases({
        "Reading input files",
        "Preparing sample",
        "Calculating d2",
        "Inferring metapop structure",
        "Writing outputs"});
  } else {
    params.progress.ConfigurePhases({
        "Reading input files",
        "Preparing sample",
        "Calculating d2",
        "Estimating Ne",
        "Writing outputs"});
  }
#ifdef GONE_NCURSES_TUI
  if (!params.quiet) {
    StartNcursesTui(&params.progress);
  }
#endif

  if (!params.realNeFile.empty()) {
    std::vector<double> refNe = LoadReferenceNe(params.realNeFile);
    if (!refNe.empty()) {
      params.progress.SetReferenceNe(refNe);
    }
  }

  int j, j3;
  double start, stop;

  std::string fichProgress = params.fileOut + "_GONE_progress.tmp";
  g_progressTmpFile = fichProgress;
  std::atexit(RemoveProgressTmpFile);
  g_prevTerminate = std::set_terminate(TerminateCleanup);

  if (params.numThreads > 0) {
    omp_set_num_threads(params.numThreads);
  }

#ifndef GONE_NCURSES_TUI
  if (!params.quiet) {
    std::cout << " Progress information stored at " << fichProgress << "\n";
    std::cout << " Check it by issuing 'cat " << fichProgress << "'\n";
    std::cout << " Loading data files" << std::endl;
  }
#endif

  start = omp_get_wtime();
  const int kTotalPhases = !params.mix
      ? 5
      : (params.twoPassGA ? 6 : 5);
  params.progress.InitTotalTasks(kTotalPhases, fichProgress.c_str());
  params.progress.SetCurrentTask(0, "Reading input files");
  params.progress.InitCurrentTask(1);
  params.progress.SetStatusDetail(params.fich);
  params.progress.SaveProgress();

  PopulationInfo* popInfo = new PopulationInfo();
  PopulationInfoMix* popInfoMix = new PopulationInfoMix();
  SampleInfo* sampleInfo = new SampleInfo();

  sampleInfo->haplotype = params.haplotype;
  sampleInfo->mix = params.mix;
  popInfo->haplotype = params.haplotype;
  sampleInfo->hayrecentbins = params.hayrecentbins;
  popInfo->hayrecentbins = params.hayrecentbins;
  popInfo->basecallcorrec = 1;

  if (params.ftype == "tped") {
    if (!ReadTped(params.fich, popInfo)) {
      delete popInfo;
      delete popInfoMix;
      delete sampleInfo;
      exit(EXIT_FAILURE);
    }
  } else if (params.ftype == "vcf") {
    if (!ReadVcf(params.fich, popInfo)) {
      delete popInfo;
      delete popInfoMix;
      delete sampleInfo;
      exit(EXIT_FAILURE);
    }
    if (params.cMMb == 0) {
      std::string fichmap = params.fich.substr(0, params.fich.rfind(".")) + ".map";
      if (!ReadMap(fichmap, popInfo)) {
        std::cerr <<
          "\nNo recombination rate has been specified (option -r) "
          "and no map file present (" << fichmap << ")" << std::endl;
        exit(EXIT_FAILURE);
      }
      if (popInfo->Mtot == 0) {
        std::cerr << "\nFound map file \"" << fichmap <<
          "\" but no genetic map was found"
          << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } else if (params.ftype == "ped") {
    std::string fichped = params.fich;
    std::string fichmap = params.fich.substr(0, params.fich.rfind(".")) + ".map";
    if (!ReadFile(fichped, fichmap, popInfo)) {
      delete popInfo;
      delete popInfoMix;
      delete sampleInfo;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cerr << "Unknown file extension." << std::endl;
    exit(EXIT_FAILURE);
  }

  params.progress.SetTaskProgress(1);
  double tProcessFile = omp_get_wtime() - start;
#ifdef GONE_NCURSES_TUI
  (void)tProcessFile;
#endif
#ifndef GONE_NCURSES_TUI
  if (!params.quiet) {
    std::cout << " Reading the input file/s took "
              << std::fixed << std::setprecision(2) << tProcessFile << " sec"
              << std::endl;
  }
#endif

  params.progress.SetCurrentTask(1, "Preparing sample");
  params.progress.InitCurrentTask(4);

  params.progress.SetStatusDetail("Shuffling loci and individuals");
  params.progress.SaveProgress();

  for (j = 0; j < popInfo->numLoci; ++j) {
    sampleInfo->shuffledLoci[j] = j;
  }
  std::shuffle(sampleInfo->shuffledLoci,
               &(sampleInfo->shuffledLoci[popInfo->numLoci]), g);

  if (popInfo->haplotype != 2) {
    for (j = 0; j < popInfo->numIndi; ++j) {
      sampleInfo->shuffledIndi[j] = j;
    }
    std::shuffle(sampleInfo->shuffledIndi,
                 &(sampleInfo->shuffledIndi[popInfo->numIndi]), g);
  } else {
    const int ndips = popInfo->numIndi / 2;
    for (j = 0; j < ndips; ++j) {
      sampleInfo->shuffledIndi[j] = j;
    }
    std::shuffle(sampleInfo->shuffledIndi,
                 &(sampleInfo->shuffledIndi[ndips]), g);
    for (j = ndips - 1; j >= 0; --j) {
      sampleInfo->shuffledIndi[j * 2]     = sampleInfo->shuffledIndi[j] * 2;
      sampleInfo->shuffledIndi[j * 2 + 1] = sampleInfo->shuffledIndi[j * 2] + 1;
    }
  }

  if (params.numSNPs > 0) {
    sampleInfo->numLoci = params.numSNPs;
  } else {
    sampleInfo->numLoci = popInfo->numLoci;
  }

  if (params.numSample > 0) {
    sampleInfo->sampleSizeIndi = params.numSample;
    if (popInfo->haplotype == 2){sampleInfo->sampleSizeIndi *= 2;}
  } else {
    sampleInfo->sampleSizeIndi = popInfo->numIndi;
  }
  if (sampleInfo->sampleSizeIndi > popInfo->numIndi) {
    sampleInfo->sampleSizeIndi = popInfo->numIndi;
  }
  
  if (sampleInfo->numLoci > popInfo->numLoci) {
    sampleInfo->numLoci = popInfo->numLoci;
  }

  params.progress.SetStatusDetail("Calculating allele frequencies");
  params.progress.SetTaskProgress(1);
  CalculateFrequencies(popInfo, sampleInfo);

  params.progress.SetStatusDetail("Calculating Hardy-Weinberg deviation");
  params.progress.SetTaskProgress(2);
  CalculateF(sampleInfo, popInfo);

  params.progress.SetStatusDetail("Calculating heterozygosity");
  params.progress.SetTaskProgress(3);
  if ((params.haplotype == 0) || (params.haplotype == 2)){
    CalculateHeterozigosity(popInfo, sampleInfo);
  }
  params.progress.SetStatusDetail("Sample preparation complete");
  params.progress.SetTaskProgress(4);

  popInfo->hetEspAll = sampleInfo->hetEspAll;

  params.progress.SetCurrentTask(2, "Calculating d2");
  params.progress.SetStatusDetail("Preparing d2 calculation");
  params.progress.SaveProgress();
#ifndef GONE_NCURSES_TUI
  if (!params.quiet) {
    std::cout << " Measuring d²" << std::endl;
  }
#endif

  CalculateD2Parallel(popInfo, sampleInfo, &params);
 
  std::stringstream salida2;
  salida2
      << params.haplotype
      << "\t# Phase (0:unphased diploids; 1:haploids; 2:phased diploids; 3:low_coverage)\n";
  salida2 << sampleInfo->numIndiEff
          << "\t# Sample size (x2 in phased diploids)\n";
  salida2 << sampleInfo->f << "\t# Hardy-Weinberg deviation in the sample\n";
  salida2 << sampleInfo->binExtra << "\t# Numero de bins extras\n";
  for (j3 = 0; j3 <= sampleInfo->binMax; ++j3) {
    if (sampleInfo->nxc[j3] > 0) {
      salida2 << std::fixed << std::setprecision(0) << sampleInfo->nxc[j3]
              << "\t" << std::fixed << std::setprecision(10)
              << sampleInfo->xc[j3] << "\t" << sampleInfo->d2[j3] << "\n";
    }
  }
  std::string fichsal2 = params.fileOut + "_d2.txt";
  std::ofstream outputFile2;
  outputFile2.open(fichsal2);
  outputFile2 << salida2.str();
  outputFile2.close();

  if (params.cMMb == 0 && popInfo->Mtot == 0) {
    std::cerr << " There is no recombination map." << std::endl;
    delete popInfo;
    delete popInfoMix;
    delete sampleInfo;
    exit(EXIT_FAILURE);
  }
  if (!params.mix) {
#ifndef GONE_NCURSES_TUI
    if (!params.quiet) {
      std::cout << " Estimating Ne" << std::endl;
    }
#endif
    params.progress.SetCurrentTask(3, "Estimating Ne");
    params.progress.SetStatusDetail("Preparing GA rounds");
    params.progress.SaveProgress();
  }

  gone(&params, fichsal2, argc, argv, popInfo, sampleInfo);
  if (!params.analysisFailReason.empty()) {
    params.twoPassGA = false;
  }
  if (params.mix && params.analysisFailReason.empty()) {
    params.progress.SetCurrentTask(3, "Inferring metapop structure");
    params.progress.SetStatusDetail(
        "Computing inter-/intra-chromosome d² and searching Fst, ps, m");
    params.progress.SaveProgress();
    RunResampledMix(popInfo, popInfoMix, sampleInfo, &params);

    if (params.twoPassGA) {
      params.progress.SetCurrentTask(4,
          "Fitting N_T");
      params.progress.SetStatusDetail(
          "Running GA on Fst-corrected d²(c)");
      params.progress.SaveProgress();
#ifndef GONE_NCURSES_TUI
      if (!params.quiet) {
        std::cout << " Fitting Ne(t)" << std::endl;
      }
#endif
      std::remove(fichsal2.c_str());
      std::stringstream salida3;
      salida3 << params.haplotype
              << "\t# Phase (0:unphased diploids; 1:haploids; "
                 "2:phased diploids; 3:low_coverage)\n";
      salida3 << sampleInfo->numIndiEff
              << "\t# Sample size (x2 in phased diploids)\n";
      salida3 << sampleInfo->f << "\t# Hardy-Weinberg deviation in the sample\n";
      salida3 << sampleInfo->binExtra << "\t# Number of extra bins\n";
      for (j3 = 0; j3 < sampleInfo->binMax; ++j3) {
        if (sampleInfo->nxc[j3] != 0) {
          salida3 << sampleInfo->nxc[j3] << "\t" << sampleInfo->xc[j3] << "\t"
                  << sampleInfo->d2[j3];
          if (sampleInfo->structuredMix) {
            salida3 << "\t" << sampleInfo->structScale[j3] << "\t"
                    << sampleInfo->structOffset[j3];
          }
          salida3 << "\n";
        }
      }
      std::ofstream outputFile3;
      outputFile3.open(fichsal2);
      outputFile3 << salida3.str();
      outputFile3.close();
      params.mix = false;
      gone(&params, fichsal2, argc, argv, popInfo, sampleInfo);
      params.mix = true;
    }
  }
  std::remove(fichsal2.c_str());

  if (params.mix && !params.chartNTMix.empty() && !params.twoPassGA) {
    params.progress.SetNeSnapshot(params.chartNTMix, "N_T");
  }

  params.progress.SetCurrentTask(kTotalPhases - 1, "Writing outputs");

  if (!params.mix &&
      (params.haplotype == 0 || params.haplotype == 2) &&
      sampleInfo->numSegLoci > 0 && sampleInfo->numIndiEff > 1) {
    const double fis = popInfo->f;
    const double fisThreshold = std::max(
        0.02, 2.0 / std::sqrt(static_cast<double>(sampleInfo->numSegLoci)));
    if (fis > fisThreshold) {
      int tailIdx[3] = {-1, -1, -1};
      for (int j = 0; j < sampleInfo->binMax; ++j) {
        if (sampleInfo->nxc[j] <= 0 || sampleInfo->d2[j] <= 0) continue;
        const double cj = sampleInfo->xc[j];
        for (int k = 0; k < 3; ++k) {
          if (tailIdx[k] == -1 || cj > sampleInfo->xc[tailIdx[k]]) {
            for (int m = 2; m > k; --m) tailIdx[m] = tailIdx[m - 1];
            tailIdx[k] = j;
            break;
          }
        }
      }
      double tailSum = 0.0;
      int tailCount = 0;
      for (int k = 0; k < 3; ++k) {
        if (tailIdx[k] != -1) {
          tailSum += sampleInfo->d2[tailIdx[k]];
          ++tailCount;
        }
      }
      const double tailMean = tailCount > 0 ? tailSum / tailCount : 0.0;
      const double nh = 2.0 * sampleInfo->numIndiEff;
      const double nh1 = nh - 1.0;
      const double sampleY = (2.0 * nh - 4.0) / (nh1 * nh1);
      const double tailRatio = sampleY > 0 ? tailMean / sampleY : 0.0;
      std::ostringstream warn;
      warn << std::fixed;
      warn << " WARNING: Fis = " << std::setprecision(4) << fis
           << " exceeds " << fisThreshold
           << " — sample may come from a structured population. "
           << "d2 tail / sampling-noise floor = " << std::setprecision(2)
           << tailRatio
           << (tailRatio > 2.0 ? " (supports structure hypothesis)."
                               : " (does not confirm structure).")
           << " Consider re-running with -x.";
      std::cerr << warn.str() << std::endl;
      params.progress.SetWarning(warn.str());
    }
  }

  params.progress.InitCurrentTask(2);
  params.progress.SetStatusDetail(
      "Writing statistics and cleaning temporary files");
  params.progress.SaveProgress();

  stop = omp_get_wtime() - start;
  std::stringstream salida;
  salida << "# (GONE v2.0)\n";
  salida << "# Command:";
  for (int i = 0; i < argc; ++i) {
    salida << " " << argv[i];
  }
  salida << "\n";
  salida << "# Running time:";
  salida << static_cast<float>(stop) << "sec\n";
  salida << "#\n";
  if (!params.analysisFailReason.empty()) {
    salida << "# ANALYSIS FAILED — NO Ne ESTIMATES WERE PRODUCED.\n";
    WriteReasonLines(salida, "#   ", params.analysisFailReason);
    salida << "#   Only this _STATS file was written; there is no _Ne and no"
              " _d2 file.\n";
    salida << "#   The statistics below describe the input data only.\n";
    salida << "#\n";
  }

  salida << "# PREPROCESSING INFORMATION:\n";
  if (params.haplotype == 0){
      salida << "# Temporal Ne estimation using unphased diploid data under the assumption of\n";
  }
  else if (params.haplotype == 1){
      salida << "# Temporal Ne estimation using haploid data under the assumption of\n";
  }
  else if (params.haplotype == 2){
      salida << "# Temporal Ne estimation using phased diploid data under the assumption of\n";
  }
  else if (params.haplotype == 3){
      salida << "# Temporal Ne estimation using diploid low-coverage data under the assumption of\n";
  }
  if (params.mix){
      salida << "#   a metapopulation model with panmixia within subpopulations.\n";
  }
  else{
      salida << "#   a single population model with panmixia.\n";
  }
  if (params.cMMb > 0){
      salida << "# The marker locations in a physical map are known and the chromosome sizes and the\n";
      salida << "#   total genome size will be calculated using those locations. The physical map will be\n";
      salida << "#   converted into a genetic map using a constant recombination rate across the genome: \n";
      salida << "#   distances in Morgans will be calculated using the physical distances (in Mb).\n";
      salida << "# Ne will be inferred using the recombination rates between loci and the weighted cuadratic\n";
      salida << "#   correlations of alleles of loci pairs.\n";
  } 
  else{
      salida << "# The marker locations in a genetic map are known and the chromosome sizes and \n";
      salida << "#   the total genetic size will be calculated using those locations. The genetic \n";
      salida << "#   distances (in Morgans) between loci pairs will be calculated directly using\n";
      salida << "#   the locations in the genetic map.\n";
      salida << "# Ne will be inferred using the recombination rates between loci and the weighted cuadratic\n";
      salida << "#   correlations of alleles of loci pairs.\n";
  }
  if (params.mix){
      salida << "# The historical Ne series will inferred by assuming that observed correlation between\n";
      salida << "#   sites with recombination rate c reflects the Ne of 1/(2c) generations ago.\n#\n";
  }
  else{
      salida << "# The complete theoretical model will be used to search for the historical Ne series\n";
      salida << "#   that best fit the the observed correlations across bins of recomination.\n#\n";
  }
  salida << "# INPUT PARAMETERS:\n";
  salida << "# Type of genotyping data. 0:unphased diploids, 1:haploids, 2:phased diploids, 3:pseudohaploids:\n";
  salida << std::fixed << std::setprecision(0);
  salida << params.haplotype << "\n";
  if ((params.haplotype == 0) || (params.haplotype == 2)){
    if (params.coverage>1){
      salida << "# Coverage (depth of DNA sequencing):\n";
      salida << std::fixed << std::setprecision(1);
      salida << params.coverage << "\n";
    }
  }
  salida << "# Rate of base call errors:\n";
  salida << std::fixed << std::setprecision(5);
  salida << params.basecallerror << "\n";
  salida << "# Number of chromosomes:\n";
  salida << std::fixed << std::setprecision(0);
  salida << popInfo->numCromo << "\n";
  salida << "# Genome size in Morgans:\n";
  salida << std::fixed << std::setprecision(4);
  salida << popInfo->Mtot << "\n";
  salida << "# Genome size in Mb:\n";
  salida << std::fixed << std::setprecision(2);
  salida << popInfo->Mbtot << "\n";
  salida << "# Total number of individuals in the input file:\n";
  salida << std::fixed << std::setprecision(0);
  if (sampleInfo->haplotype == 2){
    salida << popInfo->numIndi/2 << "\n";
  }
  else{
    salida << popInfo->numIndi << "\n";
  }
  salida << "# Number of individuals included in the analysis:\n";
  if (sampleInfo->haplotype == 2){
    salida << sampleInfo->sampleSizeIndi/2 << "\n";
  }
  else{
    salida << sampleInfo->sampleSizeIndi << "\n";
  }
  salida << "# Average Number of individuals included in the analysis (corrected for missing genotypes):\n";
  salida << std::fixed << std::setprecision(4);
  if (sampleInfo->haplotype == 2){
    salida << popInfo->avgNumIndiAnalyzed/2 <<"\n";
  }
  else{
    salida << popInfo->avgNumIndiAnalyzed <<"\n";
  }
  salida << "# Effective Number of individuals for correlations (corrected for missing genotypes):\n";
  if (sampleInfo->haplotype == 2){
    salida << sampleInfo->numIndiEff/2 <<"\n";
  }
  else{
    salida << sampleInfo->numIndiEff <<"\n";
  }
  salida << "# Number of markers in the input file:\n";
  salida << std::fixed << std::setprecision(0);
  salida << popInfo->numLoci << "\n";
  salida << "# Number of SNPs included in the analysis (only polymorphic and with less than 20% missing data):\n";
  salida << std::fixed << std::setprecision(0);
  salida << sampleInfo->numSegLoci << "\n";
  salida << std::fixed << std::setprecision(8);
  salida << "# Proportion of missing data:\n";
  salida << popInfo->propMiss << "\n";
  salida << "#\n";

  if ((params.haplotype == 0) || (params.haplotype == 2)){
    salida << "# Estimated Fis value of the population (deviation from H-W "
            "proportions):\n";
    salida << popInfo->f << "\n";
    salida << "# Heterozygosity observed in the sample (only polymorphic sites):\n";
    salida << sampleInfo->hetAvg<< "\n";
    salida << "# Heterozygosity observed in the sample (all sites):\n";
    salida << sampleInfo->hetAvgAll<< "\n";
  }
  salida << "# Heterozygosity expected (H-W eq.) in the sample (only polymorphic sites):\n";
  salida << sampleInfo->hetEsp<< "\n";
  salida << "# Heterozygosity expected (H-W eq.) in the sample (all sites):\n";
  salida << sampleInfo->hetEspAll<< "\n";
  if (params.mixRun) {
    salida << "#\n";
    if (!params.mixStructureEstimated) {
      std::string why = params.mixFailReason;
      if (!params.analysisFailReason.empty()) {
        why = params.analysisFailReason;
      } else if (why.empty()) {
        why = "-x analysis could not converge.";
      }
      salida << "# NOT ESTIMATED: no metapopulation structure parameters and"
                " no Ne curve (_Ne / _d2 files)\n";
      salida << "#   were produced. Reason:\n";
      WriteReasonLines(salida, "#   ", why);
      salida << "# Number of subpopulations in the metapopulation:\n";
      salida << "NA\n";
      salida << "# Fst among subpopulations:\n";
      salida << "NA\n";
      salida << "# Migration rate:\n";
      salida << "NA\n";
    } else {
      salida << "# Number of subpopulations in the metapopulation:\n";
      salida << std::fixed << std::setprecision(0) << popInfo->ps << "\n";
      salida << "# Fst among subpopulations:\n";
      salida << std::fixed << std::setprecision(8) << popInfo->Fst << "\n";
      salida << "# Migration rate:\n";
      salida << std::fixed << std::setprecision(8) << popInfo->m << "\n";
    }
  }
  params.progress.SetTaskProgress(2);
#ifdef GONE_NCURSES_TUI
  if (!params.printToStdOut) {
    std::string fichsal = params.fileOut + Gone2OutPrefix(params.mixRun) + "STATS";
    std::ofstream outputFile(fichsal);
    outputFile << salida.str();
  }
  if (!params.quiet) {
    if (params.mix && !params.chartNTMix.empty() && !params.twoPassGA) {
      params.progress.SetNeSnapshot(params.chartNTMix, "N_T");
    }
    params.progress.SetStatusDetail(
        "Done. Outputs written to " + params.fileOut +
        Gone2OutPrefix(params.mixRun) +
        ((!params.analysisFailReason.empty() ||
          (params.mixRun && !params.mixStructureEstimated))
             ? "STATS only"
             : "{Ne,d2,STATS}") +
        ". Press Ctrl-C to exit.");
    params.progress.SaveProgress();
    WaitForSigint();
  }
  StopNcursesTui();
  if (params.printToStdOut && !params.quiet) {
    std::cout << salida.str() << std::endl;
  }
#else
  if (params.printToStdOut) {
    if (!params.quiet) {
      std::string output = salida.str();
      std::cout << output << std::endl;
    }
  } else {
    std::string fichsal = params.fileOut + Gone2OutPrefix(params.mixRun) + "STATS";
    std::ofstream outputFile;
    outputFile.open(fichsal);
    outputFile << salida.str();
    outputFile.close();
  }
#endif

  delete sampleInfo;
  delete popInfo;
  delete popInfoMix;
  return params.analysisFailReason.empty() ? 0 : EXIT_FAILURE;
}
