
#include "./bicho.hpp"


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

inline void multiplePower(double base1, double base2, double base3,
                          double base4, int exponent, double *res1,
                          double *res2, double *res3, double *res4) {
  *res1 = *res2 = *res3 = *res4 = 1;
  while (exponent > 0) {
    if (exponent & 1) {
      *res1 *= base1;
      *res2 *= base2;
      *res3 *= base3;
      *res4 *= base4;
    }
    base1 *= base1;
    base2 *= base2;
    base3 *= base3;
    base4 *= base4;
    exponent >>= 1;
  }
}

// Core of computation:
// This subroutine calculates the values of SD2_c, SW_c and the predicted d2_c
// for a particular set of Ne_t values stored in Ne[]. HAPLOIDS
double CalculaSC(Bicho* bicho, GsampleInfo* sampleInfo) {
  int i, ii, exponente, hastasegmento;
  double Sd2, A;
  double Sw;

  double Ne12base, Ne122base, Ne102base;
  double Ne12ancho, Ne122ancho, Ne102ancho;
  double oneMinuscValSqAncho;

  double score, Necons, d2c;
  double p1a, p1b, r1a, s1;
  double a, b;
  double aPlus, bPlus;
  double acuOneMinusCvalSq;

  double Nec122, Nec102;

  double Ne12nSegMinusOne, Nec122nSegMinusOne;
  int nsegmentos = bicho->nSeg;

  Necons = bicho->NeBl[nsegmentos - 1];
  score = 0;
  hastasegmento = nsegmentos;
  if (hastasegmento > 1) {
    --hastasegmento;
  }
  for (ii = 0; ii < sampleInfo->nBins; ++ii) {
    // Removed the 0th case of the hastasegmento loop to be able to remove the
    // if that was inside it
    p1a = p1b = 1;
    s1 = 0;
    r1a = 1;
    Sd2 = 0;
    Sw = 0;
    acuOneMinusCvalSq = 1;
    Ne12nSegMinusOne = 1.0 - 2.0 / Necons;
    Nec122nSegMinusOne = sampleInfo->oneMinuscValSq[ii] *
             (1.0 - 2.2 / Necons);

    for (i = 0; i < hastasegmento; ++i) {
      exponente = bicho->segBl[i + 1] - bicho->segBl[i];
      Ne12base  = 1.0 - 2.0 / bicho->NeBl[i];
      Ne122base = 1.0 - 2.2 / bicho->NeBl[i];
      Ne102base = 1.0 - 0.2 / bicho->NeBl[i];
      // This could be made faster combining all powers into 1
      //Ne12ancho = powah(Ne12base, exponente);
      //Ne122ancho = powah(Ne122base, exponente);
      //Ne102ancho = powah(Ne102base, exponente);
      //oneMinuscValSqAncho = powah(sampleInfo->oneMinuscValSq[ii], exponente);
      multiplePower(Ne12base, Ne122base, Ne102base,
                  sampleInfo->oneMinuscValSq[ii], exponente, &Ne12ancho,
                  &Ne122ancho, &Ne102ancho, &oneMinuscValSqAncho);

      Nec122 = sampleInfo->oneMinuscValSq[ii] * Ne122base;
      Nec102 = sampleInfo->oneMinuscValSq[ii] * Ne102base;
      //A = (1 - Ne122ancho * cval12ancho) / (1 - (Nec122[i]));
      a = (1.0 - Ne122ancho * oneMinuscValSqAncho) / (1.0 - Nec122);
      //B = (1 - Ne12ancho) / (1 - Ne12[i]);
      b = (1.0 - Ne12ancho) / (1.0 - Ne12base);

      //SD2 += s1 * p1a * B + p1b / Neblock[i] * cval12aa *
                                //(Ne12[i] * B - Nec122[i] * A) /
                                //(Ne12[i] - Nec122[i]);
      Sd2 += s1 * p1a * b + p1b / bicho->NeBl[i] * acuOneMinusCvalSq *
                                (Ne12base * b - Nec122 * a) /
                                (Ne12base - Nec122);
      //SW += p1a * (1 - Ne12ancho) / (2 / Neblock[i]);
      Sw += p1a * (1.0 - Ne12ancho) / (2.0 / bicho->NeBl[i]);

      //s1 += r1a * cval12aa / Neblock[i] * (1 - Ne102ancho * cval12ancho) /
            //(1 - Nec102[i]);

      s1 += r1a * acuOneMinusCvalSq / bicho->NeBl[i] *
            (1.0 - Ne102ancho * oneMinuscValSqAncho) / (1.0 - Nec102);
      r1a *= Ne102ancho;
      p1a *= Ne12ancho;
      p1b *= Ne122ancho;
      acuOneMinusCvalSq *= oneMinuscValSqAncho;
    }

    // TODO(me): This could be precomputed every ii loop
    Ne12base = 1.0 - 2.0 / Necons;
    Nec122 = sampleInfo->oneMinuscValSq[ii] *
             (1.0 - 2.2 / Necons);
    aPlus = Nec122nSegMinusOne / (1.0 - Nec122nSegMinusOne);
    bPlus = Ne12nSegMinusOne * Necons / 2.0;

    //SD2 += s1 * p1a * Necons / 2 +
             //p1b / Necons * cval12aa * (Bplus - Aplus) /
                 //(Ne12[nsegmentos - 1] - Nec122[nsegmentos - 1]);

    Sd2 += s1 * p1a * Necons / 2.0 + p1b / Necons * acuOneMinusCvalSq *
                                         (bPlus - aPlus) / (Ne12nSegMinusOne - Nec122nSegMinusOne);

    //SW += p1a * Necons / 2;
    Sw += p1a * Necons / 2.0;

    d2c = Sd2 / Sw;

    // Correccion muestreo en generaciones sucesivas
    d2c*= sampleInfo->ngensamplingcorrec[ii];

    if (sampleInfo->haplotype !=1) {
      // Correccion acumulacion diploides *(1+c^2)
      d2c *= sampleInfo->onePluscValSq[ii];
    }

    if (sampleInfo->haplotype == 1) {  // CORRECCION MUESTREO HAPLOIDES
      bicho->d2cPred[ii] =
          (d2c * sampleInfo->basecallcorrec*sampleInfo->cValRep[ii] + sampleInfo->sampleZ4 ); 
    } 
    else if (sampleInfo->haplotype == 2) { // CORRECCION MUESTREO DIPLOIDES PHASED
        bicho->d2cPred[ii] =
            (d2c * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] + sampleInfo->sampleZ4 ); 
    } 
    else if (sampleInfo->haplotype == 3){ // CORRECCION MUESTREO LOW COVERAGE (PSEUDOHAPLOIDES)
        bicho->d2cPred[ii] =  (d2c / 4 * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] + sampleInfo->sampleZ3 );    
    }
    else if (sampleInfo->haplotype == 0) {  // CORRECCION MUESTREO DIPLOIDES UNPHASED
        bicho->d2cPred[ii] =
          (d2c * sampleInfo->basecallcorrec * sampleInfo->cValRep[ii] * sampleInfo->sampleX +sampleInfo->cValSq[ii] / Necons * sampleInfo->sampleX + sampleInfo->sampleY ) / sampleInfo->correccion;
    } 

    A = sampleInfo->d2cObs[ii] - bicho->d2cPred[ii];
    score += Square<double>(A);
  }
  //for (ii = 0; ii < sampleInfo->nBins; ++ii) {
  //  A = sampleInfo->d2cObs[ii] - bicho->d2cPred[ii];
  //  score += (A*=A);
  //}
  return (score);
}
// // Core of computation:
// // This subroutine calculates the values of SD2_c, SW_c and the predicted d2_c
// // for a particular set of Ne_t values stored in Ne[]. HAPLOIDS
// double CalculaSC_not_works(Bicho* bicho, GsampleInfo* sampleInfo) {
//   int i, ii, exponente, hastasegmento;
//   double Sd2, A;
//   double Sw;

//   double Ne12base, Ne122base, Ne102base;
//   double Ne12ancho, Ne122ancho, Ne102ancho;
//   double oneMinuscValSqAncho;

//   double score, Necons, d2c;
//   double p1a, p1b, r1a, s1;
//   double a, b;
//   double aPlus, bPlus;
//   double acuOneMinusCvalSq;

//   double Nec122, Nec102;

//   double Ne12nSegMinusOne, Nec122nSegMinusOne;
//   int nsegmentos = bicho->nSeg;

//   Necons = bicho->NeBl[nsegmentos - 1];
//   score = 0;
//   hastasegmento = nsegmentos;
//   if (hastasegmento > 1) {
//     --hastasegmento;
//   }
//   for (ii = 0; ii < sampleInfo->nBins; ++ii) {
//     // Removed the 0th case of the hastasegmento loop to be able to remove the
//     // if that was inside it
//     p1a = p1b = 1;
//     s1 = 0;
//     r1a = 1;
//     Sd2 = 0;
//     Sw = 0;
//     acuOneMinusCvalSq = 1;
//     Ne12nSegMinusOne = 1.0 - 2.0 / bicho->NeBl[nsegmentos - 1];
//     Nec122nSegMinusOne = sampleInfo->oneMinuscValSq[ii] *
//              (1.0 - 2.2 / bicho->NeBl[nsegmentos - 1]);

//     for (i = 0; i < hastasegmento; ++i) {
//       exponente = bicho->segBl[i + 1] - bicho->segBl[i];
//       Ne12base  = 1.0 - 2.0 / bicho->NeBl[i];
//       Ne122base = 1.0 - 2.2 / bicho->NeBl[i];
//       Ne102base = 1.0 - 0.2 / bicho->NeBl[i];
//       // This could be made faster combining all powers into 1
//       //Ne12ancho = powah(Ne12base, exponente);
//       //Ne122ancho = powah(Ne122base, exponente);
//       //Ne102ancho = powah(Ne102base, exponente);
//       //oneMinuscValSqAncho = powah(sampleInfo->oneMinuscValSq[ii], exponente);
//       multiplePower(Ne12base, Ne122base, Ne102base,
//                   sampleInfo->oneMinuscValSq[ii], exponente, &Ne12ancho,
//                   &Ne122ancho, &Ne102ancho, &oneMinuscValSqAncho);

//       Nec122 = sampleInfo->oneMinuscValSq[ii] * Ne122base;
//       Nec102 = sampleInfo->oneMinuscValSq[ii] * Ne102base;
//       //A = (1 - Ne122ancho * cval12ancho) / (1 - (Nec122[i]));
//       a = (1.0 - Ne122ancho * oneMinuscValSqAncho) / (1.0 - Nec122);
//       //B = (1 - Ne12ancho) / (1 - Ne12[i]);
//       b = (1.0 - Ne12ancho) / (1.0 - Ne12base);

//       //SD2 += s1 * p1a * B + p1b / Neblock[i] * cval12aa *
//                                 //(Ne12[i] * B - Nec122[i] * A) /
//                                 //(Ne12[i] - Nec122[i]);
//       Sd2 += s1 * p1a * b + p1b / bicho->NeBl[i] * acuOneMinusCvalSq *
//                                 (Ne12base * b - Nec122 * a) /
//                                 (Ne12base - Nec122);
//       //SW += p1a * (1 - Ne12ancho) / (2 / Neblock[i]);
//       Sw += p1a * (1.0 - Ne12ancho) / (2.0 / bicho->NeBl[i]);

//       //s1 += r1a * cval12aa / Neblock[i] * (1 - Ne102ancho * cval12ancho) /
//             //(1 - Nec102[i]);

//       s1 += r1a * acuOneMinusCvalSq / bicho->NeBl[i] *
//             (1.0 - Ne102ancho * oneMinuscValSqAncho) / (1.0 - Nec102);
//       r1a *= Ne102ancho;
//       p1a *= Ne12ancho;
//       p1b *= Ne122ancho;
//       acuOneMinusCvalSq *= oneMinuscValSqAncho;
//     }

//     // TODO(me): This could be precomputed every ii loop
//     Ne12base = 1.0 - 2.0 / bicho->NeBl[nsegmentos - 1];
//     Nec122 = sampleInfo->oneMinuscValSq[ii] *
//              (1.0 - 2.2 / bicho->NeBl[nsegmentos - 1]);
//     aPlus = Nec122nSegMinusOne / (1.0 - Nec122nSegMinusOne);
//     bPlus = Ne12nSegMinusOne * Necons / 2.0;

//     //SD2 += s1 * p1a * Necons / 2 +
//              //p1b / Necons * cval12aa * (Bplus - Aplus) /
//                  //(Ne12[nsegmentos - 1] - Nec122[nsegmentos - 1]);

//     Sd2 += s1 * p1a * Necons / 2.0 + p1b / Necons * acuOneMinusCvalSq *
//                                          (bPlus - aPlus) / (Ne12nSegMinusOne - Nec122nSegMinusOne);

//     //SW += p1a * Necons / 2;
//     Sw += p1a * Necons / 2.0;

//     d2c = Sd2 / Sw;
//     // Correccion acumulacion diploides *(1+c^2)
//     d2c *= sampleInfo->onePluscValSq[ii];

//     if (sampleInfo->haplotype == 2) {  // Correccion diploides unphased
//       bicho->d2cPred[ii] =
//           (d2c * sampleInfo->cValRep[ii] * sampleInfo->sampleX +
//            sampleInfo->cValSq[ii] / Necons * sampleInfo->sampleX +
//            sampleInfo->sampleY) /
//           sampleInfo->correccion;
//     } else if (sampleInfo->haplotype == 1) {  // Correccion diploides phased
//       bicho->d2cPred[ii] =
//           (d2c * sampleInfo->cValRep[ii] + sampleInfo->sampleZ4);
//     } else if (sampleInfo->haplotype == 0) {  // Correccion pseudohaploides
//       bicho->d2cPred[ii] =
//           (d2c / 4 * sampleInfo->cValRep[ii] * sampleInfo->sampleZ2 +
//            sampleInfo->sampleZ3);
//     }
//     A = sampleInfo->d2cObs[ii] - bicho->d2cPred[ii];
//     score += Square<double>(A);
//   }
//   return score;
// }
