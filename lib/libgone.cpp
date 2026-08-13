
#include "libgone.hpp"

#include "ne_band.hpp"

#include <algorithm>
#include <limits>
#include <vector>

// GA driver. Loads the d² file, then (in the non-mix path) runs the
// Ne-estimation GA — one pass per config in the combo build — over
// GONE_ROUNDS parallel rounds, geometric-mean-averages the per-round Ne
// curves, keeps the best-fitting pass, pushes live updates to the chart,
// and writes the final _GONE2_Ne / _GONE2_d2 (plus the -E density TSV). In
// mix mode the GA is skipped; the caller handles the structured path.
// If fewer than kMinBinsForGA bins survive loading, nothing is estimated:
// params->analysisFailReason is set and the function returns, leaving the
// caller to report it in _GONE2_STATS.
// Params: params — run configuration and progress sink; fichero — the
//   intermediate d² file to read (its base name drives the outputs); argc/
//   argv — original command line (echoed into _GONE2_STATS); popInfo —
//   population/structure info (Fst signals the -x -T N_T scale); sInfo —
//   sample info carrying the observed d².
void gone(AppParams *params, std::string fichero, int argc, char *argv[],
          PopulationInfo *popInfo, SampleInfo *sInfo) {
  const double clow = 0, chigh = 0.5;

  if (params->hc <= params->lc) {
    std::cerr << " Invalid range of recombination frequencies." << std::endl;
    exit(1);
  }

  xprng.setSeed(params->semilla);

  const std::string fichsal = fichero.substr(0, fichero.length() - 7);
  const std::string g2pref = Gone2OutPrefix(params->mixRun);
  const std::string fichero_sal_NeH  = fichsal + g2pref + "Ne";
  const std::string fichero_sal_d2   = fichsal + g2pref + "d2";
  const std::string fichero_sal_evol = fichsal + g2pref + "evol";
  const std::string fichero_dbg      = fichsal + g2pref + "dbg";
  std::ofstream salida;

  GsampleInfo *sampleInfo = new GsampleInfo();
  sampleInfo->flags     = params->flags;
  sampleInfo->hetEspAll = popInfo->hetEspAll;
  sampleInfo->hetEsp    = popInfo->hetEsp;
  sampleInfo->hetAvgAll = popInfo->hetAvgAll;
  sampleInfo->hetAvg    = popInfo->hetAvg;
  sampleInfo->gaConfig  = MakeDefaultGAConfig();

  if (params->basecallerror == 0) {
    sampleInfo->basecallcorrec = 1;
  } else {
    sampleInfo->basecallcorrec =
        pow((1 - 4 * params->basecallerror * (1 - params->basecallerror)), 2);
  }
  popInfo->basecallcorrec = sampleInfo->basecallcorrec;

  sampleInfo->mix = sInfo->mix;
  if ((sampleInfo->flags & FLAG_RESIZE_BINS) > 0) {
    sampleInfo->sizeBins = params->sizeBins;
  } else {
    sampleInfo->sizeBins =
        params->mixRun ? DEFAULT_BIN_SIZE_MIX : DEFAULT_BIN_SIZE;
  }
  sampleInfo->nBins = kNumBins;

  const int linesRead = ProcessFile(fichero, clow, chigh, sampleInfo);
  if (linesRead < kMinBinsForGA) {
    const std::string what =
        "There are not enough recombination bins to perform the analysis: " +
        std::to_string(linesRead) + " bin(s) available, " +
        std::to_string(kMinBinsForGA) + " required.";
    const bool allLociAnalysed =
        params->numSNPs <= 0 || params->numSNPs >= popInfo->numLoci;
    const long int suggested =
        static_cast<long int>(0.8 * static_cast<double>(sInfo->numSegLoci));
    std::string hint =
        "The bins are drawn from the -l / -u recombination window, so widen"
        " that first.";
    if (allLociAnalysed && suggested >= 2) {
      hint += " Reducing the SNPs (e.g. -s " + std::to_string(suggested) +
              ", 80% of the " + std::to_string(sInfo->numSegLoci) +
              " segregating loci) can add a few bins under -x, but it cannot"
              " help when the window itself is the limit.";
    }
    params->analysisFailReason = what + "\n" + hint;
    std::cerr << " " << what << "\n";
    std::cerr << " " << hint << "\n";
    params->progress.SetWarning("WARNING: " + what + " " + hint);
    params->progress.SetStatusDetail("Not enough recombination bins");
    delete sampleInfo;
    return;
  }
  memcpy(&(sInfo->xc[0]), &(sampleInfo->cVal[0]), linesRead * sizeof(double));
  for (int i = 0; i < linesRead; ++i) {
    sInfo->nxc[i] = static_cast<long int>(sampleInfo->nBin[i]);
  }
  memcpy(&(sInfo->d2[0]), &(sampleInfo->d2cObs[0]), linesRead * sizeof(double));
  sInfo->binMax = linesRead;

  sampleInfo->cValMin = MAX_DOUBLE;

  if (sampleInfo->haplotype == 0) {
    sampleInfo->correccion = 1.0 / Square<double>(1.0 + sampleInfo->fVal);
  }

  sampleInfo->muestraSalida = params->muestraSalida;

  CalculateSumNBins(sampleInfo);
  CalculateAverageNe(sampleInfo);

  const double g = params->ngensampling;
  for (int i = 0; i < sampleInfo->nBins; ++i) {
    const int indx = sampleInfo->indx[i];
    const double ac = 1 - sampleInfo->cVal[indx];
    if (g == 1) {
      sampleInfo->ngensamplingcorrec[indx] = 1;
    } else {
      sampleInfo->ngensamplingcorrec[indx] =
          (g - 2 * ac + 2 * pow(ac, g + 1) - ac * ac * g) /
          (g * g * (ac - 1) * (ac - 1));
    }
  }

  if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
    salida.open(fichero_sal_evol, std::ios::out);
    salida << "Gener\tSCbest\tSCmed1\tnsegbest\tnsegmed\n";
    salida.close();
  }

  if (!params->mix) {
    Pool *pool = new Pool();
    SetInitialPoolParameters(pool, sampleInfo->cValMin);
    pool->gaConfig = sampleInfo->gaConfig;
    PrePopulatePool(pool, sampleInfo);

    if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
      salida.open(fichero_dbg, std::ios::app);
      for (int j = 0; j < pool->parents[0].nSeg; ++j) {
        salida << j << "\t" << pool->parents[0].segBl[j] << "\t"
               << pool->parents[0].NeBl[j] / 2.0 << "\n";
      }
      salida << "\n";
      salida.close();
    }

    std::vector<GAConfig> configs;
#ifdef GONE_COMBO
    configs.push_back(MakeComboTruncKickConfig());
    configs.push_back(MakeComboL2Config());
    configs.push_back(MakeComboL1KickConfig());
#else
    configs.push_back(sampleInfo->gaConfig);
#endif

    double avgD2Pred[kNumLinMax] = {};
    double avgNe[MAXBINS] = {};
    double bestAvgD2Pred[kNumLinMax] = {};
    double bestAvgNe[MAXBINS] = {};
    int    bestGmax2 = 0;
    double bestResidual = std::numeric_limits<double>::infinity();
    std::string bestLabel;
    std::vector<std::vector<double>> passNeSamples;
    std::vector<std::vector<double>> bestNeSamples;

    static const Pool emptyPool = Pool();
    int numThreads = params->numThreads;
    if (numThreads == 0) {
      numThreads = omp_get_max_threads();
    }

    double **tavgNe     = new double *[numThreads];
    double **tavgD2Pred = new double *[numThreads];
    int    *maxNeConta  = new int[numThreads]{};
    for (int z = 0; z < numThreads; ++z) {
      tavgD2Pred[z] = new double[kNumLinMax]{};
      tavgNe[z]     = new double[kNumGenMax]{};
      std::fill_n(tavgNe[z], kNumGenMax, 1.0);
    }
    params->progress.InitCurrentTask(
        static_cast<float>(GONE_ROUNDS * configs.size()));
    std::vector<double> liveLogNe(MAXBINS, 0.0);
    std::vector<int> liveNeCounts(MAXBINS, 0);
    int completedRounds = 0;
    int liveMaxConta = 0;
    const std::string neUnitLabel =
        params->haplotype == 1 ? "Ne_haploids" : "Ne_diploids";

    for (size_t passIdx = 0; passIdx < configs.size(); ++passIdx) {
      const GAConfig& cfg = configs[passIdx];
      sampleInfo->gaConfig = cfg;
      xprng.setSeed(params->semilla);
      std::fill_n(avgD2Pred, kNumLinMax, 0.0);
      std::fill_n(avgNe, MAXBINS, 0.0);
      for (int z = 0; z < numThreads; ++z) {
        std::fill_n(tavgD2Pred[z], kNumLinMax, 0.0);
        std::fill_n(tavgNe[z], kNumGenMax, 1.0);
        maxNeConta[z] = 0;
      }
      std::fill(liveLogNe.begin(), liveLogNe.end(), 0.0);
      std::fill(liveNeCounts.begin(), liveNeCounts.end(), 0);
      if (params->showNeBand) {
        passNeSamples.assign(MAXBINS, std::vector<double>());
      }
      int passCompleted = 0;
      liveMaxConta = 0;
      double passBestResidual = std::numeric_limits<double>::infinity();
      params->progress.ResetBestScore();
      const std::string passPrefix =
          configs.size() > 1
              ? "[" + std::to_string(passIdx + 1) + "/" +
                std::to_string(configs.size()) + " " +
                std::string(cfg.label) + "] "
              : std::string();

#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      Pool *privpool = new Pool();

      for (int _i = tid; _i < GONE_ROUNDS; _i += numThreads) {
        xprng.setSeed(params->semilla + _i);
        *privpool = emptyPool;
        SetInitialPoolParameters(privpool, sampleInfo->cValMin);
        privpool->gaConfig = sampleInfo->gaConfig;
        PrePopulatePool(privpool, sampleInfo);
        if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
          RunDbg(privpool, sampleInfo, fichsal);
        } else {
          Run(privpool, sampleInfo, fichsal);
        }

        Bicho *bestBicho = &privpool->parents[0];
        double bestD2Pred[kMaxD2PredBins] = {};
        CalculaSC(bestBicho, sampleInfo, bestD2Pred);
        for (int i = 0; i < sampleInfo->nBins; ++i) {
          tavgD2Pred[tid][i] += bestD2Pred[i] / GONE_ROUNDS;
        }

#ifdef GA_LOG_ROUNDS
        {
          const std::string fichRounds = fichsal + "_GA_rounds.csv";
#pragma omp critical(ga_log_rounds)
          {
            static bool headerWritten = false;
            std::ofstream out(fichRounds, std::ios::app);
            if (!headerWritten) {
              out << "round,scval,nSeg,ne_gen1,ne_gen10,ne_gen100\n";
              headerWritten = true;
            }
            auto neAtGen = [&](int gen) -> double {
              int seg = 0;
              while (seg < bestBicho->nSeg &&
                     bestBicho->segBl[seg + 1] <= gen) ++seg;
              if (seg >= bestBicho->nSeg) seg = bestBicho->nSeg - 1;
              const double scale = params->haplotype == 1 ? 1.0 : 0.5;
              return std::max(MIN_NE_SIZE, bestBicho->NeBl[seg] * scale);
            };
            out << (_i + 1) << "," << bestBicho->SCval << ","
                << bestBicho->nSeg << "," << neAtGen(1) << ","
                << neAtGen(10) << "," << neAtGen(100) << "\n";
          }
        }
#endif

        int conta = 0;
        for (int i = 0; i < bestBicho->nSeg; ++i) {
          for (int j = bestBicho->segBl[i]; j < bestBicho->segBl[i + 1]; ++j) {
            ++conta;
          }
        }
        if (conta > kNumLinMax) conta = kNumLinMax;
        const int conta2 = conta;
        double sumNe[kNumGenMax] = {0};

        for (int i = 0; i < conta; ++i) sumNe[i] = 1;
        for (int ii = 0; ii < sampleInfo->muestraSalida; ++ii) {
          conta = 0;
          for (int i = 0; i < privpool->parents[ii].nSeg; ++i) {
            for (int j = privpool->parents[ii].segBl[i];
                 j < privpool->parents[ii].segBl[i + 1]; ++j) {
              sumNe[conta] *= pow(privpool->parents[ii].NeBl[i],
                                  1.0 / sampleInfo->muestraSalida);
              ++conta;
              if (conta > privpool->poolParams.gmax[0]) break;
            }
            if (conta > privpool->poolParams.gmax[0]) break;
          }
        }

        for (int i = 0; i < conta2; ++i) {
          sumNe[i] /= 2;
          tavgNe[tid][i] *= pow(sumNe[i], 1.0 / GONE_ROUNDS);
        }
        std::vector<double> roundNe(conta2, MIN_NE_SIZE);
        for (int i = 0; i < conta2; ++i) {
          const double reportedNe =
              params->haplotype == 1 ? sumNe[i] * 2.0 : sumNe[i];
          roundNe[i] = std::max(reportedNe, MIN_NE_SIZE);
        }
        if (conta2 > maxNeConta[tid]) maxNeConta[tid] = conta2;
#pragma omp critical(gone_progress_update)
        {
          ++completedRounds;
          ++passCompleted;
          liveMaxConta = std::max(liveMaxConta, conta2);
          for (int i = 0; i < conta2; ++i) {
            liveLogNe[i] += std::log(roundNe[i]);
            liveNeCounts[i] += 1;
          }
          if (params->showNeBand) {
            for (int i = 0; i < conta2; ++i) {
              passNeSamples[i].push_back(roundNe[i]);
            }
          }
          std::vector<double> liveNe(liveMaxConta, MIN_NE_SIZE);
          for (int i = 0; i < liveMaxConta; ++i) {
            if (liveNeCounts[i] > 0) {
              liveNe[i] = std::exp(liveLogNe[i] / liveNeCounts[i]);
            }
          }
          double roundRes = 0.0, roundN = 0.0;
          for (int i = 0; i < sampleInfo->nBins; ++i) {
            if (sampleInfo->cVal[i] != 0) {
              const double d = sampleInfo->d2cObs[i] - bestD2Pred[i];
              roundRes += sampleInfo->nBin[i] * d * d;
              roundN   += sampleInfo->nBin[i];
            }
          }
          const double roundScore =
              roundN > 0
                  ? roundRes / roundN
                  : std::numeric_limits<double>::infinity();
          if (roundScore < passBestResidual) {
            passBestResidual = roundScore;
          }
          params->progress.SetStatusDetail(
              passPrefix + "Completed GA round " +
              std::to_string(passCompleted) + " of " +
              std::to_string(GONE_ROUNDS));
          if (popInfo->Fst > 0.0 && popInfo->Fst < 1.0) {
            params->progress.SetNeSnapshot(liveNe, "N_T");
          } else {
            params->progress.SetNeSnapshot(liveNe, neUnitLabel);
          }
          if (params->showNeBand) {
            const NeBand band = ComputeNeBand(passNeSamples);
            params->progress.SetNeBand(band.min, band.max);
            params->progress.SetNeDensity(passNeSamples);
          }
          params->progress.SetBestScore(roundScore, passCompleted,
                                        GONE_ROUNDS);
          params->progress.SetTaskProgress(
              static_cast<float>(completedRounds));
        }
      }
      delete privpool;
    }

    std::fill_n(avgNe, MAXBINS, 1.0);
    for (int nt = 0; nt < numThreads; ++nt) {
      for (int i = 0; i < maxNeConta[nt]; ++i) {
        avgNe[i] *= tavgNe[nt][i];
      }
      for (int i = 0; i < sampleInfo->nBins; ++i) {
        avgD2Pred[i] += tavgD2Pred[nt][i];
      }
    }

    {
      int curveLen = 0;
      for (int nt = 0; nt < numThreads; ++nt) {
        if (maxNeConta[nt] > curveLen) curveLen = maxNeConta[nt];
      }
      const int maxSeg = std::min(curveLen, kNumLinMax - 1);
      if (maxSeg >= 1) {
        Bicho synth = {};
        synth.nSeg = maxSeg;
        for (int i = 0; i < maxSeg; ++i) {
          synth.segBl[i] = i;
          synth.NeBl[i] = std::max(MIN_NE_SIZE, avgNe[i] * 2.0);
        }
        synth.segBl[maxSeg] = maxSeg;
        synth.efval = sampleInfo->fVal;
        double predD2[kMaxD2PredBins] = {};
        CalculaSC(&synth, sampleInfo, predD2);
        double res = 0.0, n = 0.0;
        const bool binWeight = (GA_BIN_WEIGHTED != 0);
        for (int i = 0; i < sampleInfo->nBins; ++i) {
          if (sampleInfo->cVal[i] != 0) {
            const double d = sampleInfo->d2cObs[i] - predD2[i];
            const double w = binWeight ? sampleInfo->nBin[i] : 1.0;
            res += w * d * d;
            n   += w;
          }
        }
        passBestResidual =
            (n > 0) ? res / n : std::numeric_limits<double>::infinity();
      }
    }

    if (passBestResidual < bestResidual) {
      bestResidual = passBestResidual;
      bestLabel    = cfg.label ? cfg.label : "";
      bestGmax2    = pool->poolParams.gmax[2];
      std::copy(avgNe, avgNe + MAXBINS, bestAvgNe);
      std::copy(avgD2Pred, avgD2Pred + kNumLinMax, bestAvgD2Pred);
      if (params->showNeBand) bestNeSamples = passNeSamples;
    }
    }

    for (int nt = 0; nt < numThreads; ++nt) {
      delete[] tavgNe[nt];
      delete[] tavgD2Pred[nt];
    }
    delete[] tavgNe;
    delete[] tavgD2Pred;
    delete[] maxNeConta;

    if (configs.size() > 1) {
      std::vector<double> finalNe(bestGmax2 + 1, MIN_NE_SIZE);
      for (int i = 0; i <= bestGmax2 && i < kNumGenMax; ++i) {
        const double scale = params->haplotype == 1 ? 2.0 : 1.0;
        finalNe[i] = std::max(MIN_NE_SIZE, bestAvgNe[i] * scale);
      }
      const bool useNTLabel =
          popInfo->Fst > 0.0 && popInfo->Fst < 1.0;
      if (useNTLabel) {
        params->progress.SetNeSnapshot(finalNe, "N_T");
      } else {
        params->progress.SetNeSnapshot(finalNe, neUnitLabel);
      }
      params->progress.SetStatusDetail(
          "Best fit: " + bestLabel +
          " (bin-weighted d² residual " + std::to_string(bestResidual) + ")");
      params->progress.SetBestScoreLabeled(
          bestResidual, "Best fit", "kept: " + bestLabel);
      std::cerr << " Combo kept: " << bestLabel
                << " (bin-weighted d² residual "
                << std::to_string(bestResidual) << ")\n";
    }

    if (params->showNeBand && !bestNeSamples.empty()) {
      const NeBand band = ComputeNeBand(bestNeSamples);
      params->progress.SetNeBand(band.min, band.max);
      params->progress.SetNeDensity(bestNeSamples);
      std::ofstream dens(fichsal + g2pref + "density.tsv", std::ios::out);
      dens << FormatDensityTsv(bestNeSamples);
      dens.close();
    }

    salida.open(fichero_sal_NeH, std::ios::out);
    if (params->haplotype == 1) {
      salida << "Generation\tNe_haploids\n";
      for (int j = 0; j < bestGmax2 + 1; ++j) {
        const int generacion = j + 1;
        if (generacion < 151) {
          bestAvgNe[j] *= 2;
          salida << generacion << "\t" << std::max(bestAvgNe[j], MIN_NE_SIZE)
                 << "\n";
        }
      }
    } else {
      const bool useNTLabel =
          popInfo->Fst > 0.0 && popInfo->Fst < 1.0;
      if (useNTLabel) {
        salida << "Generation\tN_T\n";
      } else {
        salida << "Generation\tNe_diploids\n";
      }
      for (int j = 0; j < bestGmax2 + 1; ++j) {
        const int generacion = j + 1;
        if (generacion < 151) {
          const double NT =
              std::max(bestAvgNe[j], MIN_NE_SIZE);
          salida << generacion << "\t" << NT << "\n";
        }
      }
    }
    salida.close();

    salida.open(fichero_sal_d2, std::ios::out);
    salida << "c_bin\tnumber_of_SNP_pairs\tObserved_d2\tPredicted_d2\n";
    for (int i = 0; i < linesRead; ++i) {
      if (sampleInfo->cVal[i] != 0) {
        salida << std::fixed << std::setprecision(8) << sampleInfo->cVal[i]
               << "\t" << std::fixed << std::setprecision(0)
               << sampleInfo->nBin[i] << "\t" << std::fixed
               << std::setprecision(8) << sampleInfo->d2cObs[i] << "\t"
               << bestAvgD2Pred[i] << "\n";
      }
    }
    salida.close();

    delete pool;
  }
  delete sampleInfo;
}
