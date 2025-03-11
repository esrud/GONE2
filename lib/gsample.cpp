
#include "gsample.hpp"
//std::random_device seed;
// mersenne_twister_engine 64bit (very good but very big)
//std::mt19937_64 rng(seed());
Xoshiro256plus xprng(0);
// uses the result of the engine to generate uniform dist
//std::uniform_real_distribution<> uniforme01(0.0,1.0);

void SortSampleBycVal(GsampleInfo* sampleInfo, int linesRead) {
  int maxcIdx;
  double maxc;
  int temp;
  for (int i = 0; i < (sampleInfo->sampleSize - 1); ++i) {
      maxcIdx = i;
      maxc = sampleInfo->cVal[sampleInfo->indx[maxcIdx]];
      for (int j = i + 1; j < linesRead; ++j) {
        if (maxc < sampleInfo->cVal[sampleInfo->indx[j]]) {
          maxcIdx = j;
          maxc = sampleInfo->cVal[sampleInfo->indx[maxcIdx]];
        }
      }
      temp = sampleInfo->indx[i];
      sampleInfo->indx[i] = sampleInfo->indx[maxcIdx];
      sampleInfo->indx[maxcIdx] = temp;
    }
}

int CompressSample(GsampleInfo* sampleInfo, int linesRead) {
  int minBinIdx;
  double minBinN;
  int prevBin, nextBin;
  double increm;
  while (linesRead > 18) {
    // Busca el mínimo
    minBinIdx = 0;
    minBinN = sampleInfo->nBin[sampleInfo->indx[minBinIdx]];
    for (int i = 1; i < linesRead; ++i) {
      if (minBinN > sampleInfo->nBin[sampleInfo->indx[i]]) {
        minBinIdx = i;
        minBinN = sampleInfo->nBin[sampleInfo->indx[i]];
      }
    }

    // This if is not used in the other iteration of this loop
    // The question is should I use a flag to avoid repeating code
    // or should I just repeat the code
    if (minBinN > sampleInfo->sizeBins) {
      break;
    }
    // Mira cual de los adyacentes es el menor
    prevBin = minBinIdx - 1;
    nextBin = minBinIdx + 1;
    if (prevBin < 0) {
      prevBin = nextBin;
    } else {
      if (nextBin >= linesRead) {
        // This is left here assuming some sort of speed-up
      } else {
        if (sampleInfo->nBin[sampleInfo->indx[prevBin]] >
            sampleInfo->nBin[sampleInfo->indx[nextBin]]) {
          prevBin = nextBin;
        }
      }
    }

    if (prevBin < minBinIdx) {
      minBinIdx = prevBin;
      prevBin += 1;
    }
    // los suma
    increm = (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] /
                  sampleInfo->cVal[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]] /
                  sampleInfo->cVal[sampleInfo->indx[prevBin]]) /
             (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]]);

    sampleInfo->cVal[sampleInfo->indx[minBinIdx]] = 1.0 / increm;

    increm = (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] *
                  sampleInfo->d2cObs[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]] *
                  sampleInfo->d2cObs[sampleInfo->indx[prevBin]]) /
             (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]]);

    sampleInfo->d2cObs[sampleInfo->indx[minBinIdx]] = increm;

    sampleInfo->nBin[sampleInfo->indx[minBinIdx]] =
        sampleInfo->nBin[sampleInfo->indx[minBinIdx]] +
        sampleInfo->nBin[sampleInfo->indx[prevBin]];
    // comprime la lista
    --linesRead;
    for (int i = prevBin; i < linesRead; ++i) {
      sampleInfo->cVal[sampleInfo->indx[i]] =
          sampleInfo->cVal[sampleInfo->indx[i + 1]];
      sampleInfo->d2cObs[sampleInfo->indx[i]] =
          sampleInfo->d2cObs[sampleInfo->indx[i + 1]];
      sampleInfo->nBin[sampleInfo->indx[i]] =
          sampleInfo->nBin[sampleInfo->indx[i + 1]];
    }
  }

  if ((sampleInfo->flags & FLAG_NBINS) == 0) {
    CalculateNBins(sampleInfo, linesRead);
  }
  while (linesRead > sampleInfo->nBins) {
    // Busca el mínimo
    minBinIdx = 0;
    minBinN = sampleInfo->nBin[sampleInfo->indx[minBinIdx]];
    for (int i = 1; i < linesRead; ++i) {
      if (minBinN > sampleInfo->nBin[sampleInfo->indx[i]]) {
        minBinIdx = i;
        minBinN = sampleInfo->nBin[sampleInfo->indx[i]];
      }
    }
    // Mira cual de los adyacentes es el menor
    prevBin = minBinIdx - 1;
    nextBin = minBinIdx + 1;
    if (prevBin < 0) {
      prevBin = nextBin;
    } else {
      if (nextBin >= linesRead) {
      } else {
        if (sampleInfo->nBin[sampleInfo->indx[prevBin]] >
            sampleInfo->nBin[sampleInfo->indx[nextBin]]) {
          prevBin = nextBin;
        }
      }
    }
    if (prevBin < minBinIdx) {
      minBinIdx = prevBin;
      prevBin += 1;
    }
    // los suma
    increm = (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] /
                  sampleInfo->cVal[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]] /
                  sampleInfo->cVal[sampleInfo->indx[prevBin]]) /
             (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]]);
    sampleInfo->cVal[sampleInfo->indx[minBinIdx]] = 1.0 / increm;
    increm = (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] *
                  sampleInfo->d2cObs[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]] *
                  sampleInfo->d2cObs[sampleInfo->indx[prevBin]]) /
             (sampleInfo->nBin[sampleInfo->indx[minBinIdx]] +
              sampleInfo->nBin[sampleInfo->indx[prevBin]]);
    sampleInfo->d2cObs[sampleInfo->indx[minBinIdx]] = increm;
    sampleInfo->nBin[sampleInfo->indx[minBinIdx]] =
        sampleInfo->nBin[sampleInfo->indx[minBinIdx]] +
        sampleInfo->nBin[sampleInfo->indx[prevBin]];
    // comprime la lista
    --linesRead;
    for (int i = prevBin; i < linesRead; ++i) {
      sampleInfo->cVal[sampleInfo->indx[i]] =
          sampleInfo->cVal[sampleInfo->indx[i + 1]];
      sampleInfo->d2cObs[sampleInfo->indx[i]] =
          sampleInfo->d2cObs[sampleInfo->indx[i + 1]];
      sampleInfo->nBin[sampleInfo->indx[i]] =
          sampleInfo->nBin[sampleInfo->indx[i + 1]];
    }
  }
  return linesRead;
}

void CalculateSumNBins(GsampleInfo* sampleInfo) {
  sampleInfo->sumNBin = 0;
  for (int i = 0; i < sampleInfo->nBins; ++i) {
    int indx = sampleInfo->indx[i];
    sampleInfo->sumNBin += sampleInfo->nBin[indx];

    if (sampleInfo->cVal[indx] < sampleInfo->cValMin) {
      sampleInfo->cValMin = sampleInfo->cVal[indx];
    }
    if (sampleInfo->cVal[indx] > sampleInfo->cValMax) {
      sampleInfo->cValMax = sampleInfo->cVal[indx];
    }
    sampleInfo->oneMinuscValSq[indx] =
        Square<double>(1.0 - sampleInfo->cVal[indx]);
    if ((sampleInfo->flags & FLAG_REP) > 0) {
      sampleInfo->cValRep[indx] = 1.0;
    } else {
      sampleInfo->cValRep[indx] = sampleInfo->oneMinuscValSq[indx];
    }
    sampleInfo->onePluscValSq[indx] =
        (1 + Square<double>(sampleInfo->cVal[indx]));
    sampleInfo->cValSq[indx] = Square<double>(sampleInfo->cVal[indx]);
  }
}

void CalculateNBins(GsampleInfo* sampleInfo, const int linesRead) {
  // CALCULA EL No DE LINEAS FINALES A PARTIR DEL NUMERO INICIAL
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
  if ((sampleInfo->nBins < 30) && (linesRead>30))  {
    sampleInfo->nBins = 30;
  }
}

void CalculateAverageNe(GsampleInfo* sampleInfo) {
  int conta = 0;
  double dpob, Ne0;
  sampleInfo->NeMed = 0;

  for (int i = 2; i < 8; ++i) {
    int indx = sampleInfo->indx[i];
    if (sampleInfo->d2cObs[indx] > 0) {
      if ((sampleInfo->haplotype == 1) || (sampleInfo->haplotype == 2)) {// haploides o diploides phased
        dpob = (sampleInfo->d2cObs[indx] - sampleInfo->sampleZ3) \
         / (sampleInfo->sampleZ2*sampleInfo->cValRep[indx] * sampleInfo->basecallcorrec);
      } 
      else { // diploides unphased
        if (sampleInfo->mix){
          dpob=sampleInfo->d2cObs[indx];
        }
        else{
        dpob = 4.0 * (sampleInfo->d2cObs[indx] - sampleInfo->sampleZ1 ) /
               (sampleInfo->cValRep[indx] * sampleInfo->sampleZ2 * sampleInfo->basecallcorrec);
        } 
      }
      Ne0 = 2 * (sampleInfo->onePluscValSq[indx] / (dpob * 2 * (1 - sampleInfo->oneMinuscValSq[indx])) -
                 1.1 * sampleInfo->oneMinuscValSq[indx] / (1 - sampleInfo->oneMinuscValSq[indx]));
      if (Ne0 > 10) {
        sampleInfo->NeMed += Ne0;
        conta++;
      }
    }
  }
  // POBLACIONES SUBDIVIDIDAS:
  sampleInfo->NeMed = 1000;
  // FIN POBLACIONES SUBDIVIDIDAS
  if (conta > 0) {
    sampleInfo->NeMed = sampleInfo->NeMed / conta;
  } else {
    sampleInfo->NeMed = 2000;
  }
  if (sampleInfo->NeMed < 10) {
    sampleInfo->NeMed = 10;
  }
}

// TODO(me): Check if closing the file is really needed when using
//           ifstream
int ProcessFile(std::string fichero, double clow, double chigh,
                GsampleInfo* sampleInfo) {
  /*
   * Carga el archivo d2 y prepara la estructura sampleInfo con los
   * parámetros de la muestra
   */

  int nLinea = -4;
  // Start of input file processing
  std::ifstream entrada;
  entrada.open(fichero, std::ios::in);
  if (!(entrada.is_open())) {
    std::cerr << " Error opening file " << fichero << "\n";
    exit(1);
  }

  double col0, col1, col2;
  std::string line;
  while (getline(entrada, line)) {
    if (line[0] != '#' && line[0] != '*' && line[0] != '>' && line[0] != '/') {
      std::istringstream iss(line);
      if (nLinea == -4) {
        // We are at the start of the file, get the important stuff
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
        nLinea++;
      } 
      else if (nLinea == -3) {
        // On the second line extract more important stuff
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
        //        std::cout<< sampleInfo->sampleSize<< std::endl;

        // Si popInfo->haplotype = 0 (dips fase desconocida) ó = 3 (pseudohaploides) :
        //      sampleInfo->sampleSize indica el numero de individuos diploides
        // Si popInfo->haplotype = 1 (haploides) :
        //      sampleInfo->sampleSize indica el numero de individuos haploides
        // Si popInfo->haplotype = 2 (diploides phased) :
        //      sampleInfo->sampleSize indica el DOBLE del numero de individuos diploides
        double sample_size_h = 2 * sampleInfo->sampleSize;

        // sampleInfo->sampleX y sampleInfo->sampleY solo se usan en diploides con fase desconocida
        sampleInfo->sampleX =
            (Cube<double>(sample_size_h - 2.0) +
             8.0 / 5.0 * Square<double>(sample_size_h - 2.0) +
             4 * (sample_size_h - 2.0)) /
            (Cube<double>(sample_size_h - 1.0) + Square<double>(sample_size_h - 1.0));
        sampleInfo->sampleY = (2.0 * sample_size_h - 4.0) / Square<double>(sample_size_h - 1.0);

        // sampleInfo->sampleZ1 solo se usa en pseudohaploides.
        sampleInfo->sampleZ1 = (1.0 + sampleInfo->sampleSize) / (sampleInfo->sampleSize * (sampleInfo->sampleSize - 1.0));

        // sampleInfo->sampleZ2 solo se usa en pseudohaploides, dips phased y haploides, Aunque sampleSize está
        // multiplicado por dos en dips phased y haploides.
        sampleInfo->sampleZ2 = 1.0 - 0.2 / (sampleInfo->sampleSize - 1.0); // Corregido error. Era 0.8

        // sampleInfo->sampleZ4 solo se usan en dips phased y haploides, sampleSize está
        // multiplicado por dos en dips phased y haploides.
        // sampleInfo->sampleZ3 solo se usa con coverage>0
        sampleInfo->sampleZ3 = 1.0 / (sampleInfo->sampleSize - 1.0);
        if ((sampleInfo->haplotype==2)){
          sampleInfo->sampleZ4 = (1.0+sampleInfo->fVal) / (sampleInfo->sampleSize - 1.0);
        }
        else{
          sampleInfo->sampleZ4 = (1.0) / (sampleInfo->sampleSize - 1.0);
        }
        // sampleZ5 es para low coverage o pseudohaploides
        sampleInfo->sampleZ5 = (sampleInfo->sampleSize + 1.0) / (sampleInfo->sampleSize - 1.0)/sampleInfo->sampleSize;

        nLinea++;
      } 
      else if (nLinea == -2) {
        if (!(iss >> col0)) {
          std::cerr << " Format error in file " << fichero << std::endl;
          entrada.close();
          exit(1);
        }
        if (col0 < -1 || col0 > 1) {
          std::cerr << " Specify a f value in the range -1 and +1."
                    << std::endl;
          entrada.close();
          exit(1);
        }
        sampleInfo->fValSample = col0;
        nLinea++;
      } 
      else if (nLinea == -1) {
        if (!(iss >> col0)) {
          std::cerr << " Format error in file " << fichero << std::endl;
          entrada.close();
          exit(1);
        }
        sampleInfo->binExtra = col0;
        nLinea++;
      } 
      else {
        if (!(iss >> col0 >> col1 >> col2)) {
          std::cerr << " Format error in file " << fichero << std::endl;
          entrada.close();
          exit(1);
        }

        if ((col0 < 0) || (col1 < 0) || (col1 > 0.5)) {
          std::cerr << " Wrong data type. Probably some input values are out "
                       "of range."
                    << std::endl;
          entrada.close();
          exit(1);
        }

        if ((col1 >= clow) && (col1 <= chigh) && (col0 > 0)) {
          sampleInfo->nBin[nLinea] = col0;
          sampleInfo->cVal[nLinea] = col1;
          sampleInfo->d2cObs[nLinea] = col2;
          sampleInfo->indx[nLinea] = nLinea;

          nLinea++;
        }
      }
      if (nLinea == kNumLinMax) {
        break;
      }
    }
  }
  // std::cout << "nLinea:"<<nLinea<<"\n"; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  entrada.close();

  if (sampleInfo->haplotype==0){
    sampleInfo->fVal =
      (1.0 + sampleInfo->fValSample * (2.0 * sampleInfo->sampleSize- 1.0)) /
      (2.0 * sampleInfo->sampleSize - 1.0 + sampleInfo->fValSample);
  }
  else if(sampleInfo->haplotype==2){
    sampleInfo->fVal =
     (1.0 + sampleInfo->fValSample * (sampleInfo->sampleSize - 1.0)) /
      (sampleInfo->sampleSize - 1.0 + sampleInfo->fValSample);
  }
  else{
    sampleInfo->fVal = 0;
  }

  // sampleInfo->fVal =
  //     (1.0 + sampleInfo->fValSample * (2 * sampleInfo->sampleSize - 1.0)) /
  //     (2 * sampleInfo->sampleSize - 1 + sampleInfo->fValSample);
  SortSampleBycVal(sampleInfo, nLinea);
  nLinea = CompressSample(sampleInfo, nLinea);
  return nLinea;
}
