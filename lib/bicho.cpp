
#include "./bicho.hpp"

#include <cmath>

#ifndef GA_SMOOTH_LAMBDA
#define GA_SMOOTH_LAMBDA 0.0001
#endif

// Fast integer exponentiation by squaring.
// Params: base — the base; exponent — non-negative integer power.
// Returns: base raised to exponent.
inline double powah(double base, int exponent) {
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

// Core of computation:
// Calculates the predicted d²_c for a piecewise-constant Ne history
// described by `bicho` and accumulates the squared residual against
// the observed d²_c stored in `sampleInfo`.
//
// Within one call, every segment-only quantity (Ne12base, Ne122base,
// Ne102base, their powers over the segment width, and b = (1-Ne12^w)/
// (1-Ne12base)) is invariant across recombination bins. We compute
// them once into local arrays and reuse them inside the per-bin loop;
// only the c-dependent powers of (1-c)² remain per-bin work.
// Params: bicho — candidate piecewise-constant Ne history; sampleInfo —
//   observed per-bin d²; d2cPred — if non-null, receives the predicted d²
//   per bin; cutoff — early-out bound (return once the partial score
//   exceeds it, used to prune clearly-worse candidates).
// Returns: the bin-weighted squared-residual score (lower is better).
static double CalculaSCImpl(Bicho* bicho, GsampleInfo* sampleInfo,
                            double* d2cPred, double cutoff) {
  const int nsegmentos = bicho->nSeg;
  const double Necons = bicho->NeBl[nsegmentos - 1];

  int hastasegmento = nsegmentos;
  if (hastasegmento > 1) {
    --hastasegmento;
  }

  int    segExp[kNumLinMax];
  double segNeBl[kNumLinMax];
  double Ne12Base[kNumLinMax];
  double Ne122Base[kNumLinMax];
  double Ne102Base[kNumLinMax];
  double Ne12Pow[kNumLinMax];
  double Ne122Pow[kNumLinMax];
  double Ne102Pow[kNumLinMax];
  double bConst[kNumLinMax];
  double swCoeff[kNumLinMax];

  for (int i = 0; i < hastasegmento; ++i) {
    const int   w   = bicho->segBl[i + 1] - bicho->segBl[i];
    const double Ne = bicho->NeBl[i];
    segExp[i]    = w;
    segNeBl[i]   = Ne;
    Ne12Base[i]  = 1.0 - 2.0 / Ne;
    Ne122Base[i] = 1.0 - 2.2 / Ne;
    Ne102Base[i] = 1.0 - 0.2 / Ne;
    Ne12Pow[i]   = powah(Ne12Base[i],  w);
    Ne122Pow[i]  = powah(Ne122Base[i], w);
    Ne102Pow[i]  = powah(Ne102Base[i], w);
    bConst[i]    = (1.0 - Ne12Pow[i]) / (1.0 - Ne12Base[i]);
    swCoeff[i]   = (1.0 - Ne12Pow[i]) / (2.0 / Ne);
  }

  const double Ne12nSegMinusOne_base = 1.0 - 2.0 / Necons;
  const double oneMinus2_2overNecons = 1.0 - 2.2 / Necons;
  const double Necons_div2 = Necons / 2.0;

  const int haplotype = sampleInfo->haplotype;
  const bool mix      = sampleInfo->mix;

  double score = 0;

  {
    const GAConfig& cfg = sampleInfo->gaConfig;
    if (cfg.smoothLambda > 0.0) {
      double smooth = 0;
      const double truncSq = cfg.smoothCutoff * cfg.smoothCutoff;
      for (int i = 1; i < nsegmentos; ++i) {
        const double a = std::max(bicho->NeBl[i - 1], MIN_NE_SIZE);
        const double b = std::max(bicho->NeBl[i],     MIN_NE_SIZE);
        const double d = std::log10(b) - std::log10(a);
        switch (cfg.smoothKind) {
          case GASmoothKind::L1:
            smooth += std::fabs(d);
            break;
          case GASmoothKind::Truncated: {
            const double dSq = d * d;
            smooth += (dSq < truncSq) ? dSq : truncSq;
            break;
          }
          case GASmoothKind::Quadratic:
          default:
            smooth += d * d;
            break;
        }
      }
      score += cfg.smoothLambda * smooth;
    }
  }

#if GA_BIN_WEIGHTED
  const bool useBinWeight = true;
  double invMeanNBin = 1.0;
  if (useBinWeight) {
    double sumNBinAll = 0;
    for (int ii = 0; ii < sampleInfo->nBins; ++ii) {
      sumNBinAll += sampleInfo->nBin[ii];
    }
    invMeanNBin = (sumNBinAll > 0) ? (sampleInfo->nBins / sumNBinAll) : 1.0;
  }
#endif

  for (int ii = 0; ii < sampleInfo->nBins; ++ii) {
    const double cv = sampleInfo->oneMinuscValSq[ii];
    double p1a = 1, p1b = 1, r1a = 1, s1 = 0;
    double Sd2 = 0, Sw = 0;
    double acuOneMinusCvalSq = 1;

    const double Nec122nSegMinusOne = cv * oneMinus2_2overNecons;

    for (int i = 0; i < hastasegmento; ++i) {
      const double Ne12ancho  = Ne12Pow[i];
      const double Ne122ancho = Ne122Pow[i];
      const double Ne102ancho = Ne102Pow[i];
      const double oneMinuscValSqAncho =
          sampleInfo->oneMinuscValSqPow[ii][segExp[i]];

      const double Ne12base  = Ne12Base[i];
      const double Ne122base = Ne122Base[i];
      const double Nec122    = cv * Ne122base;
      const double Nec102    = cv * Ne102Base[i];

      const double a = (1.0 - Ne122ancho * oneMinuscValSqAncho) / (1.0 - Nec122);
      const double b = bConst[i];

      Sd2 += s1 * p1a * b + p1b / segNeBl[i] * acuOneMinusCvalSq *
                                (Ne12base * b - Nec122 * a) /
                                (Ne12base - Nec122);
      Sw += p1a * swCoeff[i];

      s1 += r1a * acuOneMinusCvalSq / segNeBl[i] *
            (1.0 - Ne102ancho * oneMinuscValSqAncho) / (1.0 - Nec102);
      r1a *= Ne102ancho;
      p1a *= Ne12ancho;
      p1b *= Ne122ancho;
      acuOneMinusCvalSq *= oneMinuscValSqAncho;
    }

    const double aPlus = Nec122nSegMinusOne / (1.0 - Nec122nSegMinusOne);
    const double bPlus = Ne12nSegMinusOne_base * Necons_div2;

    Sd2 += s1 * p1a * Necons_div2 + p1b / Necons * acuOneMinusCvalSq *
                                       (bPlus - aPlus) /
                                       (Ne12nSegMinusOne_base - Nec122nSegMinusOne);
    Sw += p1a * Necons_div2;

    double d2c = Sd2 / Sw;

    d2c *= sampleInfo->ngensamplingcorrec[ii];
    if (haplotype != 1) {
      d2c *= sampleInfo->onePluscValSq[ii];
    }

    double pred;
    switch (haplotype) {
      case 1:
      case 2:
        pred = d2c * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] +
               sampleInfo->sampleZ4;
        break;
      case 3:
        pred = d2c / 4 * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] +
               sampleInfo->sampleZ3;
        break;
      default:
        if (mix) {
          pred = d2c * sampleInfo->basecallcorrec;
        } else {
          pred = (d2c * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] *
                      sampleInfo->sampleX +
                  sampleInfo->cValSq[ii] / Necons * sampleInfo->sampleX +
                  sampleInfo->sampleY) /
                 sampleInfo->correccion;
        }
        break;
    }
    if (sampleInfo->structuredMix) {
      pred = sampleInfo->structScale[ii] * pred + sampleInfo->structOffset[ii];
    }
    if (d2cPred != nullptr) {
      d2cPred[ii] = pred;
    }

    const double residual = sampleInfo->d2cObs[ii] - pred;
#if GA_BIN_WEIGHTED
    if (useBinWeight) {
      score += sampleInfo->nBin[ii] * invMeanNBin * Square<double>(residual);
    } else {
      score += Square<double>(residual);
    }
#else
    score += Square<double>(residual);
#endif
    if (d2cPred == nullptr && score > cutoff) {
      return cutoff;
    }
  }
  return score;
}

// Score a candidate and also write its predicted d² per bin.
// Params: bicho — candidate Ne history; sampleInfo — observed d²; d2cPred —
//   receives the predicted d² per bin.
// Returns: the bin-weighted squared-residual score.
double CalculaSC(Bicho* bicho, GsampleInfo* sampleInfo, double* d2cPred) {
  return CalculaSCImpl(bicho, sampleInfo, d2cPred, MAX_DOUBLE);
}

// Score a candidate without producing the predicted-d² output.
// Params: bicho — candidate Ne history; sampleInfo — observed d².
// Returns: the bin-weighted squared-residual score.
double CalculaSCScoreOnly(Bicho* bicho, GsampleInfo* sampleInfo) {
  return CalculaSCImpl(bicho, sampleInfo, nullptr, MAX_DOUBLE);
}

// Score a candidate with an early-out: stop once the partial score exceeds
// cutoff (used to prune candidates that can't beat the current best).
// Params: bicho — candidate Ne history; sampleInfo — observed d²; cutoff —
//   score bound to bail out at.
// Returns: the score, or a value ≥ cutoff once the bound is crossed.
double CalculaSCScoreCutoff(Bicho* bicho, GsampleInfo* sampleInfo,
                            double cutoff) {
  return CalculaSCImpl(bicho, sampleInfo, nullptr, cutoff);
}
