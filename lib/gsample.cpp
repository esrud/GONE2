
#include "gsample.hpp"

thread_local Xoshiro256plus xprng(0);

// Fast integer exponentiation by squaring.
// Params: base — the base; exponent — non-negative integer power.
// Returns: base raised to exponent.
static inline double PowIntSameOrder(double base, int exponent) {
  double result = 1;
  while (exponent > 0) {
    if (exponent & 1) {
      result *= base;
    }
    base *= base;
    exponent >>= 1;
  }
  return result;
}

// Selection sort of sampleInfo->indx[] by descending cVal over all
// linesRead bins. Kept as selection sort for output stability —
// std::sort would tie-break differently, which can shift bin contents
// in CompressSample.
// Params: sampleInfo — holds cVal[] and the indx[] permutation reordered
//   in place; linesRead — number of bins to sort.
void SortSampleBycVal(GsampleInfo* sampleInfo, int linesRead) {
  for (int i = 0; i < (linesRead - 1); ++i) {
    int maxcIdx = i;
    double maxc = sampleInfo->cVal[sampleInfo->indx[maxcIdx]];
    for (int j = i + 1; j < linesRead; ++j) {
      const double cv = sampleInfo->cVal[sampleInfo->indx[j]];
      if (maxc < cv) {
        maxcIdx = j;
        maxc = cv;
      }
    }
    std::swap(sampleInfo->indx[i], sampleInfo->indx[maxcIdx]);
  }
}

// Merge the lowest-count bin with its smaller adjacent neighbour (combining
// counts, harmonic-mean recombination rate, count-weighted d²) and shift
// the trailing entries down by one.
// Params: sampleInfo — per-bin arrays modified in place; linesRead —
//   current bin count; minBinIdx — sorted-order index of the bin to merge.
// Returns: the new bin count (linesRead - 1).
static int MergeSmallestBin(GsampleInfo* sampleInfo, int linesRead,
                            int minBinIdx) {
  int prevBin = minBinIdx - 1;
  int nextBin = minBinIdx + 1;
  if (prevBin < 0) {
    prevBin = nextBin;
  } else if (nextBin < linesRead &&
             sampleInfo->nBin[sampleInfo->indx[prevBin]] >
                 sampleInfo->nBin[sampleInfo->indx[nextBin]]) {
    prevBin = nextBin;
  }

  if (prevBin < minBinIdx) {
    minBinIdx = prevBin;
    prevBin += 1;
  }

  const int aIdx = sampleInfo->indx[minBinIdx];
  const int bIdx = sampleInfo->indx[prevBin];
  const double nA = sampleInfo->nBin[aIdx];
  const double nB = sampleInfo->nBin[bIdx];
  const double cA = sampleInfo->cVal[aIdx];
  const double cB = sampleInfo->cVal[bIdx];
  const double d2A = sampleInfo->d2cObs[aIdx];
  const double d2B = sampleInfo->d2cObs[bIdx];
  const double totalN = nA + nB;

  const double invMean = (nA / cA + nB / cB) / totalN;
  sampleInfo->cVal[aIdx] = 1.0 / invMean;
  sampleInfo->d2cObs[aIdx] = (nA * d2A + nB * d2B) / totalN;
  sampleInfo->nBin[aIdx] = totalN;

  --linesRead;
  for (int i = prevBin; i < linesRead; ++i) {
    const int dst = sampleInfo->indx[i];
    const int src = sampleInfo->indx[i + 1];
    sampleInfo->cVal[dst]   = sampleInfo->cVal[src];
    sampleInfo->d2cObs[dst] = sampleInfo->d2cObs[src];
    sampleInfo->nBin[dst]   = sampleInfo->nBin[src];
  }
  return linesRead;
}

// Locate the bin (in sorted order) with the smallest sample count.
// Params: sampleInfo — per-bin arrays; linesRead — number of bins.
// Returns: the sorted-order index of the bin with the fewest sample pairs.
static int FindMinBin(const GsampleInfo* sampleInfo, int linesRead) {
  int minBinIdx = 0;
  double minBinN = sampleInfo->nBin[sampleInfo->indx[0]];
  for (int i = 1; i < linesRead; ++i) {
    const double n = sampleInfo->nBin[sampleInfo->indx[i]];
    if (minBinN > n) {
      minBinIdx = i;
      minBinN = n;
    }
  }
  return minBinIdx;
}

// Compact the recombination bins by repeatedly merging the smallest one:
// first until every remaining bin has enough sample pairs (≥18 bins kept),
// then down to the target bin count sampleInfo->nBins.
// Params: sampleInfo — per-bin arrays modified in place; linesRead —
//   starting bin count.
// Returns: the final bin count after compaction.
int CompressSample(GsampleInfo* sampleInfo, int linesRead) {
  while (linesRead > 18) {
    const int minBinIdx = FindMinBin(sampleInfo, linesRead);
    if (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] > sampleInfo->sizeBins) {
      break;
    }
    linesRead = MergeSmallestBin(sampleInfo, linesRead, minBinIdx);
  }

  if ((sampleInfo->flags & FLAG_NBINS) == 0) {
    CalculateNBins(sampleInfo, linesRead);
  }

  while (linesRead > sampleInfo->nBins) {
    const int minBinIdx = FindMinBin(sampleInfo, linesRead);
    linesRead = MergeSmallestBin(sampleInfo, linesRead, minBinIdx);
  }
  return linesRead;
}

// Sum the per-bin sample-pair counts into sumNBin and precompute the
// per-bin recombination quantities the scorer reuses (powers of (1-c)²,
// (1+c²), c², the cVal min/max, etc.).
// Params: sampleInfo — reads nBin[]/cVal[]/flags, writes sumNBin and the
//   precomputed per-bin arrays.
void CalculateSumNBins(GsampleInfo* sampleInfo) {
  sampleInfo->sumNBin = 0;
  const bool repFlag = (sampleInfo->flags & FLAG_REP) > 0;
  for (int i = 0; i < sampleInfo->nBins; ++i) {
    const int indx = sampleInfo->indx[i];
    sampleInfo->sumNBin += sampleInfo->nBin[indx];

    const double cv = sampleInfo->cVal[indx];
    if (cv < sampleInfo->cValMin) {
      sampleInfo->cValMin = cv;
    }
    if (cv > sampleInfo->cValMax) {
      sampleInfo->cValMax = cv;
    }
    const double oneMinusCSq = Square<double>(1.0 - cv);
    sampleInfo->oneMinuscValSq[indx] = oneMinusCSq;
    for (int exponent = 0; exponent < kNumGenMax; ++exponent) {
      sampleInfo->oneMinuscValSqPow[indx][exponent] =
          PowIntSameOrder(oneMinusCSq, exponent);
    }
    sampleInfo->cValRep[indx] = repFlag ? 1.0 : oneMinusCSq;
    sampleInfo->onePluscValSq[indx] = 1.0 + Square<double>(cv);
    sampleInfo->cValSq[indx] = Square<double>(cv);
  }
}

// Choose the target number of recombination bins from the smallest bin's
// pair count, clamped to [30, 60] (and never above linesRead). Stored in
// sampleInfo->nBins.
// Params: sampleInfo — writes nBins; linesRead — current bin count.
void CalculateNBins(GsampleInfo* sampleInfo, const int linesRead) {
  double producto = log10(sampleInfo->nBin[sampleInfo->indx[linesRead - 1]]);
  producto = producto - 3.0;
  if (producto < 0.0) {
    producto = 0.0;
  }
  producto = pow(2.0, producto) * 10;
  sampleInfo->nBins = static_cast<int>(producto + 8);
  if (sampleInfo->nBins > linesRead) {
    sampleInfo->nBins = linesRead;
  }
  if (sampleInfo->nBins > 60) {
    sampleInfo->nBins = 60;
  }
  if ((sampleInfo->nBins < 30) && (linesRead > 30)) {
    sampleInfo->nBins = 30;
  }
}

// Estimate a starting average Ne (NeMed) from the observed d² of a few
// mid-range recombination bins, using the genotyping-mode-appropriate d²
// correction. Falls back to a default if no bin yields a usable estimate.
// Params: sampleInfo — reads per-bin d²/corrections, writes NeMed.
void CalculateAverageNe(GsampleInfo* sampleInfo) {
  int conta = 0;
  double dpob, Ne0;
  sampleInfo->NeMed = 0;

  for (int i = 2; i < 8; ++i) {
    const int indx = sampleInfo->indx[i];
    if (sampleInfo->d2cObs[indx] <= 0) continue;

    if ((sampleInfo->haplotype == 1) || (sampleInfo->haplotype == 2)) {
      dpob = (sampleInfo->d2cObs[indx] - sampleInfo->sampleZ3) /
             (sampleInfo->sampleZ2 * sampleInfo->cValRep[indx] *
              sampleInfo->basecallcorrec);
    } else if (sampleInfo->mix) {
      dpob = sampleInfo->d2cObs[indx];
    } else {
      dpob = 4.0 * (sampleInfo->d2cObs[indx] - sampleInfo->sampleZ1) /
             (sampleInfo->cValRep[indx] * sampleInfo->sampleZ2 *
              sampleInfo->basecallcorrec);
    }
    Ne0 = 2 * (sampleInfo->onePluscValSq[indx] /
                   (dpob * 2 * (1 - sampleInfo->oneMinuscValSq[indx])) -
               1.1 * sampleInfo->oneMinuscValSq[indx] /
                   (1 - sampleInfo->oneMinuscValSq[indx]));
    if (Ne0 > 10) {
      sampleInfo->NeMed += Ne0;
      ++conta;
    }
  }
  if (conta > 0) {
    sampleInfo->NeMed = sampleInfo->NeMed / conta;
  } else {
    sampleInfo->NeMed = 2000;
  }
  if (sampleInfo->NeMed < 10) {
    sampleInfo->NeMed = 10;
  }
}

// Read the intermediate d² file (header lines plus per-bin rows of pair
// count, recombination rate, observed d², and optional -x -T structure
// columns) into sampleInfo, keeping only bins with c in [clow, chigh], then
// sort and compact them into the final recombination bins.
// Params: fichero — path to the d² file; clow/chigh — recombination-rate
//   window to keep; sampleInfo — filled in place.
// Returns: the number of bins retained after sorting and compaction.
int ProcessFile(std::string fichero, double clow, double chigh,
                GsampleInfo* sampleInfo) {
  std::ifstream entrada;
  entrada.open(fichero, std::ios::in);
  if (!(entrada.is_open())) {
    std::cerr << " Error opening file " << fichero << "\n";
    exit(1);
  }

  int nLinea = -4;
  double col0, col1, col2;
  std::string line;
  sampleInfo->structuredMix = false;
  for (int i = 0; i < kNumLinMax; ++i) {
    sampleInfo->structScale[i]  = 1.0;
    sampleInfo->structOffset[i] = 0.0;
  }
  while (getline(entrada, line)) {
    if (line.empty()) continue;
    const char first = line[0];
    if (first == '#' || first == '*' || first == '>' || first == '/') {
      continue;
    }
    std::istringstream iss(line);
    if (nLinea == -4) {
      if (!(iss >> col0)) {
        std::cerr << " Format error in file " << fichero << std::endl;
        entrada.close();
        exit(1);
      }
      if (col0 != 0 && col0 != 1 && col0 != 2 && col0 != 3) {
        std::cerr << " Specify a valid code for type of genotyping data "
                     "(0, 1, 2 or 3)."
                  << std::endl;
        entrada.close();
        exit(1);
      }
      sampleInfo->haplotype = static_cast<int>(col0);
      ++nLinea;
    } else if (nLinea == -3) {
      if (!(iss >> col0)) {
        std::cerr << " Format error in file " << fichero << std::endl;
        entrada.close();
        exit(1);
      }
      if (col0 < 2) {
        std::cerr << " Specify a sample size larger than 1." << std::endl;
        entrada.close();
        exit(1);
      }
      sampleInfo->sampleSize = col0;

      const double n  = sampleInfo->sampleSize;
      const double nm = n - 1.0;
      const double nh = 2.0 * n;
      const double nh1 = nh - 1.0;
      const double nh2 = nh - 2.0;

      sampleInfo->sampleX = (Cube<double>(nh2) +
                              8.0 / 5.0 * Square<double>(nh2) +
                              4 * nh2) /
                             (Cube<double>(nh1) + Square<double>(nh1));
      sampleInfo->sampleY = (2.0 * nh - 4.0) / Square<double>(nh1);

      sampleInfo->sampleZ1 = (1.0 + n) / (n * nm);
      sampleInfo->sampleZ2 = 1.0 - 0.2 / nm;
      sampleInfo->sampleZ3 = 1.0 / nm;
      sampleInfo->sampleZ4 = (sampleInfo->haplotype == 2)
                                 ? (1.0 + sampleInfo->fVal) / nm
                                 : 1.0 / nm;
      sampleInfo->sampleZ5 = (n + 1.0) / nm / n;

      ++nLinea;
    } else if (nLinea == -2) {
      if (!(iss >> col0)) {
        std::cerr << " Format error in file " << fichero << std::endl;
        entrada.close();
        exit(1);
      }
      if (col0 < -1 || col0 > 1) {
        std::cerr << " Specify a f value in the range -1 and +1." << std::endl;
        entrada.close();
        exit(1);
      }
      sampleInfo->fValSample = col0;
      ++nLinea;
    } else if (nLinea == -1) {
      if (!(iss >> col0)) {
        std::cerr << " Format error in file " << fichero << std::endl;
        entrada.close();
        exit(1);
      }
      sampleInfo->binExtra = col0;
      ++nLinea;
    } else {
      if (!(iss >> col0 >> col1 >> col2)) {
        std::cerr << " Format error in file " << fichero << std::endl;
        entrada.close();
        exit(1);
      }
      double col3, col4;
      const bool hasStruct = static_cast<bool>(iss >> col3 >> col4);
      if ((col0 < 0) || (col1 < 0) || (col1 > 0.5)) {
        std::cerr << " Wrong data type. Probably some input values are out "
                     "of range."
                  << std::endl;
        entrada.close();
        exit(1);
      }
      if ((col1 >= clow) && (col1 <= chigh) && (col0 > 0)) {
        sampleInfo->nBin[nLinea]   = col0;
        sampleInfo->cVal[nLinea]   = col1;
        sampleInfo->d2cObs[nLinea] = col2;
        sampleInfo->indx[nLinea]   = nLinea;
        if (hasStruct) {
          sampleInfo->structScale[nLinea]  = col3;
          sampleInfo->structOffset[nLinea] = col4;
          sampleInfo->structuredMix        = true;
        }
        ++nLinea;
      }
    }
    if (nLinea == kNumLinMax) break;
  }
  entrada.close();

  if (sampleInfo->haplotype == 0) {
    sampleInfo->fVal =
        (1.0 + sampleInfo->fValSample * (2.0 * sampleInfo->sampleSize - 1.0)) /
        (2.0 * sampleInfo->sampleSize - 1.0 + sampleInfo->fValSample);
  } else if (sampleInfo->haplotype == 2) {
    sampleInfo->fVal =
        (1.0 + sampleInfo->fValSample * (sampleInfo->sampleSize - 1.0)) /
        (sampleInfo->sampleSize - 1.0 + sampleInfo->fValSample);
  } else {
    sampleInfo->fVal = 0;
  }

  SortSampleBycVal(sampleInfo, nLinea);
  nLinea = CompressSample(sampleInfo, nLinea);
  return nLinea;
}
