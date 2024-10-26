
#include "sample.hpp"
#include "./pool.hpp"
#include <string.h>

void CalculateFrequencies(PopulationInfo* popInfo, SampleInfo* sampleInfo) {
  // Calculo de frecuencia del alelo noref y de homo noref para todos los loci
  // variables en la muestra
  int i, j;
  double contaIndX = 0;
  int locus;
  int indi;
  double counter = 0;
  int mincontaIndX;
  uint8_t ff,gg;


  // Looks like we need a reference to the shuffled array
  sampleInfo->hetEsp=0;
  mincontaIndX=int(float(sampleInfo->sampleSizeIndi)*2.0*(0.8));// 20% missing data
  for (j = 0; j < sampleInfo->numLoci; ++j) {
    locus = sampleInfo->shuffledLoci[j];
    sampleInfo->frec[locus] = 0;
    sampleInfo->homo[locus] = 0;
    sampleInfo->segrega[locus] = false;
    contaIndX = 0;
    if (popInfo->haplotype!=2){
      for (i = 0; i < sampleInfo->sampleSizeIndi; i++) {
        indi = sampleInfo->shuffledIndi[i];
        ff = popInfo->indi[indi][locus];
        if (ff < 9) {
          // acumulador de frecuencia del alelo noref
          sampleInfo->frec[locus] += static_cast<double>(ff);
          if (ff == 2) {
            ++sampleInfo->homo[locus];
          }  // acumulador de homo noref
          ++(++contaIndX);
        }
      }
    }
    else{
      for (i = 0; i < sampleInfo->sampleSizeIndi; ++(++i)) {
        indi = sampleInfo->shuffledIndi[i];
        ff = popInfo->indi[indi][locus];
        indi = sampleInfo->shuffledIndi[i+1];
        gg = popInfo->indi[indi][locus];
        if ((ff < 9) and (gg < 9)) {
          // acumulador de frecuencia del alelo noref
          sampleInfo->frec[locus] += static_cast<double>(ff);
          ++(++contaIndX);
          // acumulador de frecuencia del alelo noref
          sampleInfo->frec[locus] += static_cast<double>(gg);
          ++(++contaIndX);
          if ((ff == 2) && (gg ==2)) {
            ++(++sampleInfo->homo[locus]);
          }  // acumulador de homo noref
        }
      }
    }

    if (contaIndX == 0) {
      std::cerr << "There is no genotyping data for at least one SNP\n";
      exit(EXIT_FAILURE);
    }
//    if ((sampleInfo->frec[locus] > 0) && (sampleInfo->frec[locus] < contaIndX)) {
    if ((sampleInfo->frec[locus] > 0) && (contaIndX>=mincontaIndX) && (sampleInfo->frec[locus] < contaIndX)) { // excluded SNPs with > 20% of missing data
      ++sampleInfo->numSegLoci;
      sampleInfo->segrega[locus] = true;
      sampleInfo->frec[locus] /= contaIndX;
      sampleInfo->homo[locus] /= (contaIndX / 2.0);
      sampleInfo->hetEsp+=2*sampleInfo->frec[locus]*(1-sampleInfo->frec[locus]);
    }  // contador de segregantes
    popInfo->avgNumIndiAnalyzed += static_cast<double>(contaIndX);
    counter++;

    if (sampleInfo->numSegLoci >= sampleInfo->numLoci) {
      break;
    }
  }
  popInfo->avgNumIndiAnalyzed /= counter * 2.0;
  sampleInfo->hetEsp /= static_cast<double>(sampleInfo->numSegLoci);
  sampleInfo->hetEspAll = sampleInfo->hetEsp * static_cast<double>(sampleInfo->numSegLoci) / static_cast<double>(sampleInfo->numLoci);
  popInfo->hetEsp=sampleInfo->hetEsp;
  popInfo->hetEspAll = sampleInfo->hetEspAll;

}

void CalculateF(SampleInfo* sampleInfo, PopulationInfo* popInfo) {
  // Calculo de f sample y control del exceso del límite de loci:
  double acupq = 0;
  double acuPp2 = 0;
  int locus;
  for (int j = 0; j < sampleInfo->numLoci; ++j) {
    locus = sampleInfo->shuffledLoci[j];
    if (sampleInfo->segrega[locus]) {
      acuPp2 += (sampleInfo->homo[locus] -
                 Square<double>(sampleInfo->frec[locus]));
      acupq +=
          (sampleInfo->frec[locus] * (1.0 - sampleInfo->frec[locus]));
    }
  }
  sampleInfo->f = acuPp2 / acupq;
  sampleInfo->hetpq = acupq / sampleInfo->numSegLoci;
}


bool CalculateHeterozigosity(PopulationInfo* popInfo, SampleInfo* sampleInfo) {
    // Distribucion de heterocigosidades de los individuos (ver
    // Medidas_de_Parent_y_Het.docx)
    int* pIndi = sampleInfo->shuffledIndi;
    int* pLocus;

    int contaLocX;
    int i, j;
    bool isSegLocus;
    uint8_t ff,gg;

    // double a, sx2 = 0, sx3 = 0;
    // sampleInfo->hetAvg = 0;
    // for (i = 0; i < sampleInfo->sampleSizeIndi; ++i) {
    //   sampleInfo->het[*pIndi] = 0;
    //   pLocus = sampleInfo->shuffledLoci;
    //   contaLocX = 0;
    //   for (j = 0; j < sampleInfo->numSegLoci; ++j) {
    //     isSegLocus = false;
    //     while (!isSegLocus) {
    //       if (sampleInfo->segrega[*pLocus]) {
    //         ff = popInfo->indi[*pIndi][*pLocus];
    //         if (ff == 1) {
    //           // acumulador de heterocigosidad de cada indiv
    //           ++sampleInfo->het[*pIndi];
    //           ++contaLocX;
    //         } else if (ff < 9) {
    //           ++contaLocX;
    //         }
    //         isSegLocus = true;
    //       }
    //       ++pLocus;
    //     }
    //   }
    //   if (contaLocX == 0) {
    //     std::cerr << "There is no genotyping data for at least one individual"
    //               << std::endl;
    //     delete popInfo;
    //     delete sampleInfo;
    //     return false;
    //   }
    //   sampleInfo->het[*pIndi] /= contaLocX;  // het de cada individuo
    //   sampleInfo->hetAvg += sampleInfo->het[*pIndi];
    //   ++pIndi;
    // }

    if (popInfo->haplotype ==0){
      sampleInfo->hetAvg = 0;
      for (i = 0; i < sampleInfo->sampleSizeIndi; ++i) {
        double hetIndi=0;
        pLocus = sampleInfo->shuffledLoci;
        contaLocX = 0;
        for (j = 0; j < sampleInfo->numSegLoci; ++j) {
          isSegLocus = false;
          while (!isSegLocus) {
            if (sampleInfo->segrega[*pLocus]) {
              ff = popInfo->indi[*pIndi][*pLocus];
              if (ff == 1) {
                // acumulador de heterocigosidad de cada indiv
                ++hetIndi;
                ++contaLocX;
              } else if (ff < 9) {
                ++contaLocX;
              }
              isSegLocus = true;
            }
            ++pLocus;
          }
        }
        if (contaLocX == 0) {
          std::cerr << "There is no genotyping data for at least one individual"
                    << std::endl;
          delete popInfo;
          delete sampleInfo;
          return false;
        }
        hetIndi /= contaLocX;  // het de cada individuo
        sampleInfo->hetAvg += hetIndi;
        ++pIndi;
      }
      sampleInfo->hetAvg /= static_cast<double>(sampleInfo->sampleSizeIndi);
      sampleInfo->hetAvgAll = sampleInfo->hetAvg * static_cast<double>(sampleInfo->numSegLoci) / static_cast<double>(sampleInfo->numLoci);
      popInfo->hetAvg=sampleInfo->hetAvg;
      popInfo->hetAvgAll = sampleInfo->hetAvgAll;
    }
    else{ // PHASED
      sampleInfo->hetAvg = 0;
      for (i = 0; i < sampleInfo->sampleSizeIndi; ++(++i)) {
        double hetIndi=0;
        pLocus = sampleInfo->shuffledLoci;
        contaLocX = 0;
        for (j = 0; j < sampleInfo->numSegLoci; ++j) {
          isSegLocus = false;
          while (!isSegLocus) {
            if (sampleInfo->segrega[*pLocus]) {
              ff = popInfo->indi[*pIndi][*pLocus];
              gg = popInfo->indi[*pIndi+1][*pLocus];
              if (ff+gg == 2) {
                // acumulador de heterocigosidad de cada indiv
                ++hetIndi;
                ++contaLocX;
              } else if ((ff < 9) && (gg < 9)) {
                ++contaLocX;
              }
              isSegLocus = true;
            }
            ++pLocus;
          }
        }
        if (contaLocX == 0) {
          std::cerr << "There is no genotyping data for at least one individual"
                    << std::endl;
          delete popInfo;
          delete sampleInfo;
          return false;
        }
        hetIndi /= contaLocX;  // het de cada individuo
        sampleInfo->hetAvg += hetIndi;
        ++(++pIndi);
      }
      sampleInfo->hetAvg /= (static_cast<double>(sampleInfo->sampleSizeIndi) / 2.0);
      sampleInfo->hetAvgAll = sampleInfo->hetAvg * static_cast<double>(sampleInfo->numSegLoci) / static_cast<double>(sampleInfo->numLoci);
      popInfo->hetAvg=sampleInfo->hetAvg;
      popInfo->hetAvgAll = sampleInfo->hetAvgAll;
    }


    // pIndi = sampleInfo->shuffledIndi;
    // for (i = 0; i < sampleInfo->sampleSizeIndi; ++i) {
    //   // het individual dividida por het media total
    //   sampleInfo->het[*pIndi] /= sampleInfo->hetAvg;

    //   a = sampleInfo->het[*pIndi] - 1.0;
    //   sx2 += Square<double>(a);
    //   sx3 += Cube<double>(a);

    //   ++pIndi;
    // }
    // sampleInfo->hetVar = sx2 / (sampleInfo->sampleSizeIndi - 1);
    // sampleInfo->hetSesg  = sx3 / sampleInfo->sampleSizeIndi;
    // double hetDT = sqrt(sampleInfo->hetVar);
    // sampleInfo->hetSesg /= Cube<double>(hetDT);
    return true;
}

bool CalculateKinshipDistribution(PopulationInfo* popInfo,
                                  SampleInfo* sampleInfo) {
  // TODO(me): this could all problaby done inside the heterozygosity
  // calculation but the code would be much uglier. On the other side,
  // performance
  // Distribucion de parentescos de los individuos (ver
  // Medidas_de_Parent_y_Het.docx)
  double b = static_cast<double>(2 * sampleInfo->sampleSizeIndi - 1) /
      static_cast<double>(2 * (sampleInfo->sampleSizeIndi - 1));
  int* pIndi = sampleInfo->shuffledIndi;
  int* pLocus;
  int i, j;
  double contaLocX;
  double a, ff;
  bool isSegLocus;
  sampleInfo->parentAvg = 0;
  for (i = 0; i < sampleInfo->sampleSizeIndi; ++i) {
    sampleInfo->parent[*pIndi] = 0;
    pLocus = sampleInfo->shuffledLoci;
    contaLocX = 0;
    // TODO(me): all this seems too convoluted. Since we know
    // the number of segregating loci, we should be able to loop
    // only through them and not have to do that while loop
    // We could store the segregating loci in a separate array that
    // holds the index to the locus in popInfo->indi[x][locus]
    // We'll see how much of the performance is wasted in this loop
    // Also, similarly, this should be possible with the heterozygosity loop
    for (j = 0; j < sampleInfo->numSegLoci; ++j) {
      isSegLocus = false;
      while (!isSegLocus) {
        if (sampleInfo->segrega[*pLocus]) {
          ff = popInfo->indi[*pIndi][*pLocus];
          if (ff < 9) {
            ++contaLocX;
            a = static_cast<double>(ff) / 2.0;
            a -= sampleInfo->frec[*pLocus];
            sampleInfo->parent[*pIndi] += Square<double>(a);
          }
          isSegLocus = true;
        }
        ++pLocus;
      }
    }
    if (contaLocX < 100) {
      std::cerr << "There are too few genotyping data for individual " << i
                << std::endl;
      delete sampleInfo;
      delete popInfo;
      return false;
    }
    sampleInfo->parent[*pIndi] /= contaLocX;
    sampleInfo->parent[*pIndi] *= (2 * b / sampleInfo->hetpq);
    sampleInfo->parentAvg += sampleInfo->parent[*pIndi];
    ++pIndi;
  }
  sampleInfo->parentAvg /= sampleInfo->sampleSizeIndi;

  return true;
}

void CalculateD2Parallel(PopulationInfo* popInfo, SampleInfo* sampleInfo,
                 AppParams* params) {
  // Calculo de D2:

  int* valid_idx = new int[popInfo->numLoci];
  int counter = 0;
  int i, j, k;
  int _i, _irepe;
  double gBase;
  double xlc = 0, xhc = 0;
  long int refextra;
  double infx, supx;
  int bin, extrabin;
  double epsilon=0;
  double iota=0;
  double kappa=0;
  double lambda=1;
  double desde,hasta;
  int genot1,genot2;
  bool okchromsize;

  params->progress.InitCurrentTask(sampleInfo->numSegLoci - 1);

  // Inicializamos la tabla de índices válidos (en los que segrega[x] es true)
  for (int idx = 0; idx < popInfo->numLoci; idx++) {
    if (sampleInfo->segrega[sampleInfo->shuffledLoci[idx]] &&
        (sampleInfo->frec[sampleInfo->shuffledLoci[idx]] > params->MAF) &&
        ((1.0 - sampleInfo->frec[sampleInfo->shuffledLoci[idx]]) >
         params->MAF)) {
      valid_idx[counter] = sampleInfo->shuffledLoci[idx];
      counter++;
    }
  }

  gBase = 1.0 / (2.0 * params->hc);
// gBase = std::max((1.0 / (2.0 * params->hc)), 1.0);
//  binBase = trunc((gBase - 10.0) / 5.0) + 6.0;

  if (params->cMMb > 0) {
    for (i = 0; i < popInfo->numLoci; ++i) {
      popInfo->posiCM[i] = params->cMMb * popInfo->posiBP[i] / 1000000;
    }
    popInfo->Mtot = params->cMMb * popInfo->Mbtot / 100;
  }

  // Checking chomosome sizes 
  j=0;
  k=0;
  okchromsize=true;
  desde=popInfo->posiCM[j];
  for (i = 1; i < popInfo->numLoci; ++i) {
    if ((popInfo->cromo[i]!=popInfo->cromo[j]) || (i==(popInfo->numLoci-1))){
      if((popInfo->posiCM[i-1]-desde) < 20.0){
        std::cerr << "Chromosome "<< popInfo->cromo[j]<<" is too small. All chromosomes must be larger than 20cM" << std::endl;
        okchromsize=false;        
      }
      popInfo->chrsize[k] = popInfo->posiCM[i-1]-desde; // Se miden los cromosomas
      ++k;
      j=i;
      desde=popInfo->posiCM[j];
    }
  }
  if (!okchromsize){
    exit(EXIT_FAILURE);
  }

  if (params->distance == 0) {  // SIN TRANSFORMAR
    xlc = params->lc;
    xhc = params->hc;
  }
  if (params->distance == 1) {  // Haldane
    if (params->lc >= 0.5) {
      xlc = 100.0;
    } else {
      xlc = -log(1.0 - 2.0 * params->lc) / 2.0;
    }
    if (params->hc >= 0.5) {
      xhc = 200.0;
    } else {
      xhc = -log(1.0 - 2.0 * params->hc) / 2.0;
    }
  }
  if (params->distance == 2) {  // Kosambi
    if (params->lc >= 0.5) {
      xlc = 100.0;
    } else {
      xlc = -log((1.0 - 2.0 * params->lc) / (1.0 + 2.0 * params->lc)) / 4.0;
    }
    if (params->hc >= 0.5) {
      xhc = 200.0;
    } else {
      xhc = -log((1.0 - 2.0 * params->hc) / (1.0 + 2.0 * params->hc)) / 4.0;
    }
  }
  infx = xlc * 100.0;  // en cM
  supx = xhc * 100.0;  // en cM

  int locus1;
  int locus2;
  int* pIndi;

  double cenmor = 0;
  double contaIndX = 0;
  double tacuHoHo, tacuHoHetHetHo, tacuHetHet;
  uint8_t ss;

  double D, W,D2,W2;
  double recrate, generacion;

  int numThreads = params->numThreads;
  if (numThreads == 0) {
    numThreads = omp_get_max_threads();
  }

  double xW[MAXBINS] = {0}, xD2[MAXBINS] = {0};
  long int nxc[MAXBINS] = {0};
  double xc[MAXBINS] = {0};

  // Per thread accumulators
  double** txW  = new double*[numThreads];
  double** txD2 = new double*[numThreads];
  long int** tnxc = new long int*[numThreads];
  double** txc = new double*[numThreads];

  double *texpNData = new double[numThreads]{};
  double *teffNData = new double[numThreads]{};

  // Allocate each thread's array
  for (int z = 0; z < numThreads; z++) {
    txW[z] = new double[MAXBINS]{};
    txD2[z] = new double[MAXBINS]{};
    tnxc[z] = new long int[MAXBINS]{};
    txc[z] = new double[MAXBINS]{};
  }

  int* binMaxes = new int[numThreads]{};
  double expNData = 0;
  double effNData = 0;
  int increbinBase;

  if (sampleInfo->hayrecentbins){
    increbinBase=3;
  }
  else{
    increbinBase=0;
  }

  omp_set_num_threads(numThreads);
  double frecI, frecJ;
  std::sort(valid_idx, valid_idx + sampleInfo->numSegLoci);
  sampleInfo->numSNPs = 0;
  int repes=1;
  if (sampleInfo->haplotype==3){
    repes=10;
  }
//  std::cout<<"Inds: "<<sampleInfo->sampleSizeIndi<<std::endl;
#pragma omp parallel for private(D, W, D2, W2, recrate, generacion, bin, extrabin,i, j, _i, genot1, genot2, locus1, locus2, tacuHoHo, tacuHoHetHetHo, tacuHetHet, pIndi, cenmor, contaIndX, frecI, frecJ, ss,_irepe) schedule(dynamic, 1000)
  for (i = 0; i < sampleInfo->numSegLoci - 1; ++i) {
    int tid = omp_get_thread_num();
    locus1 = valid_idx[i];
    contaIndX = 0;
    double threadEffNCounter = 0;
    double threadExpNCounter = 0;
    double posiCM1 = popInfo->posiCM[locus1];
    int cromo1 = popInfo->cromo[locus1];
    bool cond;
    for (j = i + 1; j < sampleInfo->numSegLoci; ++j) {
      locus2 = valid_idx[j];
      // Are we in the same cromosome?
      if (popInfo->cromo[locus2] != cromo1) {
        continue;
      }
      cenmor = fabs(popInfo->posiCM[locus2] - posiCM1);
      
      if (sampleInfo->hayrecentbins){
        if ((cenmor > 56)) {continue;}
        if ((cenmor > 30) && (cenmor < 39)) {continue;}
        if ((cenmor > 15) && (cenmor < 20)) {continue;}
        if ((cenmor > supx) && (cenmor < 10)) {continue;}
        if ((cenmor < infx)) {continue;}
        extrabin=-1;
        if (cenmor>=39){extrabin=0;}
        else if(cenmor>=20){extrabin=1;}
        else if(cenmor>=10){extrabin=2;}
      }
      else{
        if ((cenmor > supx) || (cenmor < infx)) {continue;}
        extrabin=-1;
      }

      // acumuladores de genotipos
      D2=W2=0;
      for (_irepe=0;_irepe<repes;++_irepe){
        tacuHoHo = tacuHoHetHetHo = tacuHetHet = frecI = frecJ = 0;
        pIndi = sampleInfo->shuffledIndi;
        contaIndX = 0;
        for (_i = 0; _i < sampleInfo->sampleSizeIndi; ++_i) {
          genot1=popInfo->indi[*pIndi][locus1];
          genot2=popInfo->indi[*pIndi][locus2];
          if (sampleInfo->haplotype==3){  // Alternative for haplotype=3
            if (genot1==1){
              if (xprng.uniforme01() > 0.5){
                genot1=0;
              }
              else{
                genot1=2;
              }
            }
            if (genot2==1){
              if (xprng.uniforme01() > 0.5){
                genot2=0;
              }
              else{
                genot2=2;
              }
            }
          }
          ss = genot1 + genot2;
          if (_irepe==0){threadExpNCounter++;}
          //++texpNData[tid];
          cond = (ss < 9);
          contaIndX += cond;
          frecI += genot1 * cond;
          frecJ += genot2 * cond;
          tacuHetHet += (ss == 2) * (genot1 == 1);
          tacuHoHetHetHo += (ss == 3);
          tacuHoHo += (ss == 4); // Homo_no_ref and Homo_no_ref
          ++pIndi;
        }
        if (contaIndX == 0) {
          continue;
        }
        frecI /= 2 * contaIndX;
        frecJ /= 2 * contaIndX;
        W = frecI * frecJ;
        if (params->haplotype == 0){
          D = -2.0 * W +
            (2.0 * tacuHoHo + tacuHoHetHetHo + tacuHetHet / 2.0) / contaIndX;
        }
        else{
          D= tacuHoHo / contaIndX - frecI * frecJ;
        }
        D *= D;
        W *= (1.0 - frecI) * (1.0 - frecJ);
        D2+=D;
        W2+=W;
      }
      D2/=repes;
      W2/=repes;
      if (params->distance == 1) {
        recrate = (1.0 - exp(-0.02 * cenmor)) / 2.0;
      } else if (params->distance == 2) {
        recrate = exp(0.04 * cenmor);
        recrate = 0.5 * (recrate - 1) / (recrate + 1);
      } else {
        recrate = cenmor / 100.0;
        if (recrate > 0.5) {
          recrate = 0.5;
        }
      }
      if (extrabin < 0){
        generacion = 0.5 / ( recrate);
        bin = trunc((generacion - gBase) / 5);
        bin = bin+increbinBase;
      }
      else{
        bin=extrabin;
      }
      // Maybe store binMax in a per-thread array and then loop in the end
      // through all the local binMaxes
      if (bin > binMaxes[tid]) {
        binMaxes[tid] = bin;
      }
      threadEffNCounter += contaIndX;
      //teffNData[tid] += contaIndX;
      tnxc[tid][bin] += 1;
      txc[tid][bin] += 1/recrate; // Para calcular media harmonica
//      txc[tid][bin] += recrate; // Para calcular media Aritmetica
      txD2[tid][bin] += D2;
      txW[tid][bin] += W2;
    }
    teffNData[tid] += threadEffNCounter;
    texpNData[tid] += threadExpNCounter;
    if (tid == 0 && i % 100 == 0) {
      params->progress.SetTaskProgress(i+1);
      //params->progress.PrintProgress();
    }
  }
  // Reduction of data
  for (int nThread = 0; nThread < numThreads; ++nThread) {
    for (int nBin = 0; nBin <= binMaxes[nThread]; ++nBin) {
      nxc[nBin] += tnxc[nThread][nBin];
      xc[nBin]  += txc[nThread][nBin];
      xD2[nBin] += txD2[nThread][nBin];
      xW[nBin]  += txW[nThread][nBin];
    }

    if (binMaxes[nThread] > sampleInfo->binMax) {
      sampleInfo->binMax = binMaxes[nThread];
    }
    effNData += teffNData[nThread];
    expNData += texpNData[nThread];
  }

  sampleInfo->binExtra= increbinBase;
  refextra=nxc[increbinBase]; // Elimina los bines extras con pocos datos
  // std::cout <<"refextra "<<refextra<<std::endl;
  // for (j=0;j<sampleInfo->binMax;++j){
  //   if (nxc[j] > 0) {
  //     std::cout <<  std::setprecision(5)<<"j: "<<j<<"  nxc: "<<nxc[j]<<"  xc: "<<1/(xc[j]/nxc[j])<<"  d2: "<<xD2[j]/xW[j]<<std::endl;
  //   }
  // }

  for (i=increbinBase-1; i>=0; --i){
      if (nxc[i]<refextra*(3-i) * 0.6){
          for (j=i;j<sampleInfo->binMax;++j){
              nxc[j] = nxc[j+1];
              xc[j]  = xc[j+1];
              xD2[j] = xD2[j+1];
              xW[j]  = xW[j+1];
          }
        --sampleInfo->binMax;
        --sampleInfo->binExtra;
      }
  }

  popInfo->expNData = expNData;
  popInfo->effNData = effNData;
  double XXW=0;
  long int XXn=0;
  for (j = 0; j <= sampleInfo->binMax; ++j) {
    if (nxc[j] > 0) {
      sampleInfo->nxc[j] = nxc[j];
      sampleInfo->xc[j] = xc[j];
      sampleInfo->xc[j] /= sampleInfo->nxc[j];
      sampleInfo->xc[j] = 1 / sampleInfo->xc[j]; // media harmonica
      sampleInfo->numSNPs += sampleInfo->nxc[j];  // numero de parejas de SNPs
      sampleInfo->d2[j] = xD2[j] / xW[j];
        XXW+=xW[j];
        XXn+=nxc[j];
    }
  }

  popInfo->propMiss = 1.0 - static_cast<double>(popInfo->effNData) /
                                static_cast<double>(popInfo->expNData);
  // Numero efectivo de individuos
  sampleInfo->numIndiEff = sampleInfo->sampleSizeIndi * popInfo->effNData / popInfo->expNData;

  sampleInfo->increW=0;
  if (params->basecallerror>0){
    XXW/=XXn;
    epsilon = params->basecallerror;
    iota = (1-2/(2* sampleInfo->numIndiEff +1));
    kappa = 1/(sampleInfo->numIndiEff +1);
    sampleInfo->increW=XXW-pow(((pow(XXW,.5)- epsilon)/(1-4 * epsilon)),2);
    for (j = 0; j <= sampleInfo->binMax; ++j) {
      if (nxc[j] > 0) {
        sampleInfo->d2[j] = (xD2[j]/nxc[j]-sampleInfo->increW*kappa) / \
                            (xW[j]/nxc[j]-sampleInfo->increW*iota);
      }
    }
  }
  // if (params->coverage>0){
  //   XXW/=XXn;
  //   lambda = pow((1+pow(0.5,params->coverage-1)),2);
  //   sampleInfo->increW=XXW*(1-1/lambda);
  //   iota = (1-2/(2* sampleInfo->numIndiEff +1));
  //   kappa = 1/(sampleInfo->numIndiEff +1);
  //   for (j = 0; j <= sampleInfo->binMax; ++j) {
  //     if (nxc[j] > 0) {
  //       sampleInfo->d2[j] = (xD2[j]/nxc[j]-sampleInfo->increW*kappa) / \
  //                           (xW[j]/nxc[j]-sampleInfo->increW*iota);
  //     }
  //   }
  // }
  if (popInfo->haplotype==0){
    popInfo->f =
      (1.0 + sampleInfo->f * (2.0 * popInfo->avgNumIndiAnalyzed - 1.0)) /
      (2.0 * popInfo->avgNumIndiAnalyzed - 1.0 + sampleInfo->f);
  }
  else if(popInfo->haplotype==2){
    popInfo->f =
     (1.0 + sampleInfo->f * (popInfo->avgNumIndiAnalyzed - 1.0)) /
      (popInfo->avgNumIndiAnalyzed - 1.0 + sampleInfo->f);
  }

  // std::cout<<sampleInfo->f<<std::endl;
  // std::cout<<popInfo->avgNumIndiAnalyzed<<std::endl;
  // std::cout<<popInfo->f<<std::endl;

  // Free the rest of the memory
  delete[] tnxc;
  delete[] txc;
  delete[] txD2;
  delete[] txW;
  delete[] teffNData;
  delete[] texpNData;
  delete[] binMaxes;
  delete[] valid_idx;
}

void CalculateD2ParallelFst(std::string fichero,PopulationInfoMix* popInfoMix, PopulationInfo* popInfo, SampleInfo* sampleInfo,
                 AppParams* params) {
  // Calculo de D2:

  int* valid_idx = new int[popInfo->numLoci];
  int counter = 0;
  int i, j;
  int _i, _irepe;
  double gBase;
  // double xlc = 0, xhc = 0;
  long int refextra;
  // double infx, supx;
  int bin, extrabin;
  double epsilon=0;
  double iota=0;
  double kappa=0;
  double lambda=1;
  double desde,hasta;
  int genot1,genot2;
  bool okchromsize;
  double Dw2, Db2, DbDw;
  double Fst;

  params->progress.InitCurrentTask(sampleInfo->numSegLoci - 1);

  // Inicializamos la tabla de índices válidos (en los que segrega[x] es true)
  for (int idx = 0; idx < popInfo->numLoci; idx++) {
    if (sampleInfo->segrega[sampleInfo->shuffledLoci[idx]]) {
      valid_idx[counter] = sampleInfo->shuffledLoci[idx];
      counter++;
    }
  }

  if (params->cMMb > 0) {
    for (i = 0; i < popInfo->numLoci; ++i) {
      popInfo->posiCM[i] = params->cMMb * popInfo->posiBP[i] / 1000000;
    }
    popInfo->Mtot = params->cMMb * popInfo->Mbtot / 100;
  }

  // Checking chomosome sizes 
  j=0;
  okchromsize=true;
  desde=popInfo->posiCM[j];
  for (i = 1; i < popInfo->numLoci; ++i) {
    if ((popInfo->cromo[i]!=popInfo->cromo[j]) || (i==(popInfo->numLoci-1))){
      if((popInfo->posiCM[i-1]-desde) < 20.0){
        std::cerr << "Chromosome "<< popInfo->cromo[j]<<" is too small. All chromosomes must be larger than 20cM" << std::endl;
        okchromsize=false;        
      }
      if((popInfo->posiCM[i-1]-desde) > 1000.0){
        std::cerr << "Chromosome "<< popInfo->cromo[j]<<" is too large. All chromosomes must be smaller than 1000cM" << std::endl;
        okchromsize=false;        
      }
      j=i;
      desde=popInfo->posiCM[j];
    }
  }
  if (!okchromsize){
    exit(EXIT_FAILURE);
  }
  // xlc = -log(1.0 - 2.0 * 0.20) / 2.0; //  desde c=0.20
  // xhc = -log(1.0 - 2.0 * 0.25) / 2.0; // hasta c=0.25
  // infx = xlc * 100.0;  // en cM
  // supx = xhc * 100.0;  // en cM

  double mapdist[MAXDIST] = {0};
  double increcM = 0.01;
  
  int locus1;
  int locus2;
  int* pIndi;

  double cenmor = 0;
  double contaIndX = 0;
  double tacuHoHo, tacuHoHetHetHo, tacuHetHet;
  uint8_t ss;

  double D, W,D2,W2;
  double recrate, generacion;

  int numThreads = params->numThreads;
  if (numThreads == 0) {
    numThreads = omp_get_max_threads();
  }

  double xW[2] = {0}, xD2[2] = {0};
  long int nxc[2] = {0};
  double xc[2] = {0};

  // Per thread accumulators
  double** txW  = new double*[numThreads];
  double** txD2 = new double*[numThreads];
  long int** tnxc = new long int*[numThreads];
  double** txc = new double*[numThreads];

  double *texpNData = new double[numThreads]{};
  double *teffNData = new double[numThreads]{};

  // Allocate each thread's array
  for (int z = 0; z < numThreads; z++) {
    txW[z] = new double[2]{};
    txD2[z] = new double[2]{};
    tnxc[z] = new long int[2]{};
    txc[z] = new double[2]{};
  }

  double expNData = 0;
  double effNData = 0;
  int increbinBase;

  omp_set_num_threads(numThreads);
  double frecI, frecJ;
  std::sort(valid_idx, valid_idx + sampleInfo->numSegLoci);
  // sampleInfo->numSNPs = 0;
  int repes=1;
  //  std::cout<<"Inds: "<<sampleInfo->sampleSizeIndi<<std::endl;
  #pragma omp parallel for private(D, W, D2, W2, recrate, generacion, bin, extrabin,i, j, _i, genot1, genot2, locus1, locus2, tacuHoHo, tacuHoHetHetHo, tacuHetHet, pIndi, cenmor, contaIndX, frecI, frecJ, ss,_irepe) schedule(dynamic, 1000)
  for (i = 0; i < sampleInfo->numSegLoci - 1; ++i) {
    int tid = omp_get_thread_num();
    locus1 = valid_idx[i];
    contaIndX = 0;
    double threadEffNCounter = 0;
    double threadExpNCounter = 0;
    double posiCM1 = popInfo->posiCM[locus1];
    int cromo1 = popInfo->cromo[locus1];
    bool cond;
    for (j = i + 1; j < sampleInfo->numSegLoci; ++j) {
      locus2 = valid_idx[j];
      // Are we in the same cromosome?
      bin=-1;
      if (popInfo->cromo[locus2] != cromo1) {
        bin=0;
      }
      else{
        cenmor = fabs(popInfo->posiCM[locus2] - posiCM1);
        if (cenmor>5){ // No se consideran las parejas más próximas de 5 cM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          cenmor/=increcM;
          #pragma omp atomic
          ++mapdist[int(cenmor)];
          bin=1;
        }
      }
      if (bin>=0){
        // acumuladores de genotipos
        D2=W2=0;
        for (_irepe=0;_irepe<repes;++_irepe){
          tacuHoHo = tacuHoHetHetHo = tacuHetHet = frecI = frecJ = 0;
          pIndi = sampleInfo->shuffledIndi;
          contaIndX = 0;
          for (_i = 0; _i < sampleInfo->sampleSizeIndi; ++_i) {
            genot1=popInfo->indi[*pIndi][locus1];
            genot2=popInfo->indi[*pIndi][locus2];
            if (sampleInfo->haplotype==3){  // Alternative for haplotype=3
              if (genot1==1){
                if (xprng.uniforme01() > 0.5){
                  genot1=0;
                }
                else{
                  genot1=2;
                }
              }
              if (genot2==1){
                if (xprng.uniforme01() > 0.5){
                  genot2=0;
                }
                else{
                  genot2=2;
                }
              }
            }
            ss = genot1 + genot2;
            if (_irepe==0){threadExpNCounter++;}
            //++texpNData[tid];
            cond = (ss < 9);
            contaIndX += cond;
            frecI += genot1 * cond;
            frecJ += genot2 * cond;
            tacuHetHet += (ss == 2) * (genot1 == 1);
            tacuHoHetHetHo += (ss == 3);
            tacuHoHo += (ss == 4); // Homo_no_ref and Homo_no_ref
            ++pIndi;
          }
          if (contaIndX == 0) {
            continue;
          }
          frecI /= 2 * contaIndX;
          frecJ /= 2 * contaIndX;
          W = frecI * frecJ;
          if (params->haplotype == 0){
            D = -2.0 * W +
              (2.0 * tacuHoHo + tacuHoHetHetHo + tacuHetHet / 2.0) / contaIndX;
          }
          else{
            D= tacuHoHo / contaIndX - frecI * frecJ;
          }
          D *= D;
          W *= (1.0 - frecI) * (1.0 - frecJ);
          D2+=D;
          W2+=W;
        }
      }
      D2/=repes;
      W2/=repes;
      // if (bin == 1) {
      //   recrate = (1.0 - exp(-0.02 * cenmor)) / 2.0;
      // } 
      // else{
      //   recrate=0.5;
      // }

      threadEffNCounter += contaIndX;
      //teffNData[tid] += contaIndX;
      tnxc[tid][bin] += 1;
      // txc[tid][bin] += 1/recrate; // Para calcular media harmonica
      //      txc[tid][bin] += recrate; // Para calcular media aritmetica
      txD2[tid][bin] += D2;
      txW[tid][bin] += W2;
    }
    teffNData[tid] += threadEffNCounter;
    texpNData[tid] += threadExpNCounter;
    if (tid == 0 && i % 100 == 0) {
      params->progress.SetTaskProgress(i+1);
      //params->progress.PrintProgress();
    }
  }
  // Reduction of data
  for (int nThread = 0; nThread < numThreads; ++nThread) {
    for (int nBin = 0; nBin <= 1; ++nBin) {
      nxc[nBin] += tnxc[nThread][nBin];
      // xc[nBin]  += txc[nThread][nBin];
      xD2[nBin] += txD2[nThread][nBin];
      xW[nBin]  += txW[nThread][nBin];
    }
    effNData += teffNData[nThread];
    expNData += texpNData[nThread];
  }
  double numIndi = sampleInfo->sampleSizeIndi * effNData / expNData;
  

  double xWt = (xW[0]+xW[1])/(nxc[0]+nxc[1]);
  for (j = 0; j <= 1; ++j) {
    if (nxc[j] > 0) {
      // xc[j] /= nxc[j];
      // xc[j] = 1 / xc[j]; // media harmonica
      xW[j] /= nxc[j];
      xD2[j] /= nxc[j];
    }
  }

  // Calculo de la distribución de distancias para la integral:
  double mindistance=5/100; // No se consideran las parejas más próximas de 5 cM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  double contadist=0;
  double maxdistance=0;
  int maxdistanceindx=0;
  for (j = 0; j < MAXDIST; ++j){
    if (mapdist[j]>0){
      contadist+=mapdist[j];
      if(mapdist[j]>maxdistance){
        maxdistance=float(j)/10000; // en Morgans
        maxdistanceindx=j;
      }
    }
  }
  for (j = 0; j <=maxdistanceindx; ++j){
      mapdist[j]/=contadist;
  }

  if (popInfo->haplotype==0){
    popInfo->f =
      (1.0 + sampleInfo->f * (2.0 * numIndi - 1.0)) /
      (2.0 * numIndi - 1.0 + sampleInfo->f);
      xWt /= (1-2/(2*numIndi+1));
  }
  else if(popInfo->haplotype==2){
    popInfo->f =
     (1.0 + sampleInfo->f * (numIndi - 1.0)) /
      (numIndi - 1.0 + sampleInfo->f);
      xWt /= (1-2/(numIndi+1));
  }

  double d2s05 = xD2[0]/xW[0];
  double d2slink = xD2[1]/xW[1];
  double fp =popInfo->f;

  double fp12 = (1 + fp) * (1 + fp);
  double sample_size_h = numIndi * 2;
  double samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
  samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
  samplex += 4 * (sample_size_h - 2.0);
  samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
  double sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

  double d2;
  switch (sampleInfo->haplotype){
    case 0:
      d2= (d2s05 - sampley * fp12) / samplex;
      break;
    case 1:
      d2=(d2s05 - 1/(sample_size_h-1));
      break;
    case 2:
      d2=(d2s05 - fp12/(sample_size_h-1));
      break;
    case 3:
      d2=(d2s05 - 1/(numIndi-1));
      break;
  }

  // std::cout<< std::setprecision(11);
  // std::cout <<d2s05<<"  "<< d2slink<< std::endl;

  if (d2>0){
    // first, find the maximum possible value of Fst
    Fst = std::sqrt(d2);
    double Ne=10E9;
    double Ne05=10E9, Nelink=10E9;
    double Ne05ant=10E9;
    double Nelinkant=10E9;
    double KK=0;
    double increFst = Fst;

    for (i=0;i<5;++i){
        increFst /= 10;
        while ((Ne05-Nelink)<=0){
            if (Fst<increFst){break;}
            Ne05ant=Ne05;
            Nelinkant=Nelink;
            Fst-=increFst;
            Ne05 = MixEcuacion05(Ne05, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2s05, KK);
            Nelink = MixEcuacionlink(maxdistance,mindistance,maxdistanceindx,&mapdist[0],Nelink, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2slink, KK);
            if ((Ne05-Nelink)>0){
                Fst+=increFst;
                Ne05=Ne05ant;
                Nelink=Nelinkant;
                break;
            }
        }
    }
    Fst-=increFst;
    Ne05 = MixEcuacion05(Ne05, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2s05, KK);
    Nelink = MixEcuacionlink(maxdistance,mindistance,maxdistanceindx,&mapdist[0],Nelink, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2slink, KK);
    double difmax=Ne05-Nelink;
    double maxFst=Fst;

  // std::cout<< std::setprecision(11);
  // std::cout <<popInfo->f<<"  "<<Fst<<"  "<< Ne05<<"  "<< Nelink<< std::endl;

    // then, search for the minimum possible value of Fst
    double minFst=0;
    double difmin=0;
    Fst = 0;
    minFst=Fst;
    Ne05 = MixEcuacion05(Ne05, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2s05, KK);
    Nelink = MixEcuacionlink(maxdistance,mindistance,maxdistanceindx,&mapdist[0],Nelink, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2slink, KK);
    difmin=Ne05-Nelink;

  // std::cout<< std::setprecision(11);
  //  std::cout <<popInfo->f<<"  "<<Fst<<"  "<< Ne05<<"  "<< Nelink<< std::endl;

    //  if ((difmin*difmax)<0){
    if (difmax>0){
      if (difmin<0){
        // Next, search for the confluence of the two estimates of Ne
        increFst =(maxFst-minFst)/10;
        Fst=maxFst;
        Ne05 = MixEcuacion05(Ne05, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2s05, KK);
        Nelink = MixEcuacionlink(maxdistance,mindistance,maxdistanceindx,&mapdist[0],Nelink, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2slink, KK);
        int signo=-1;
        double dif=(Ne05-Nelink);
        int countchanges=0;
        while (countchanges<10){
            // Check for valid Fst
            while ((dif*signo<0)){
                Fst +=increFst*signo; // reset the previous value
                if (Fst>maxFst){Fst=maxFst;}
                if (Fst<minFst){Fst=minFst;}
                Ne05 = MixEcuacion05(Ne05, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2s05, KK);
                Nelink = MixEcuacionlink(maxdistance,mindistance,maxdistanceindx,&mapdist[0],Nelink, numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst, &popInfo->m, &xWt,&Dw2,&Db2,&DbDw,d2slink, KK);
                dif=Ne05-Nelink;
                // std::cout << signo<<"  "<< Fst<<"  "<<minFst<<"  "<<maxFst<<"  "<<Ne05<<"  "<< Nelink<< std::endl;
            }
            signo=-signo;
            ++countchanges;
            increFst/=10;
        }
        popInfo->Fst = Fst;
        popInfo->m = (1-Fst)/(16*Nelink*Fst);

  // std::cout<< std::setprecision(11);
  // std::cout <<popInfo->f<<"  "<<Fst<<"  "<< Ne05<<"  "<< Nelink<< std::endl;

      }
      else{
        popInfo->Fst = Fst =0;
        popInfo->m=0.5;
      }
      double d2spred=0;
      Ne = FixMixEcuacion05link(numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst,&popInfo->m,&xWt,&Dw2,&Db2,&DbDw,0.5,&d2spred,d2s05, KK);
      popInfoMix->xc[0]=0.5;
      popInfoMix->Nelink[0]=2*Ne;
      popInfoMix->mig[0]=popInfo->m;
      popInfoMix->Fst[0]=popInfo->Fst;
      popInfoMix->Dw2[0]=Dw2;
      popInfoMix->Db2[0]=Db2;
      popInfoMix->DbDw[0]=DbDw;
      popInfoMix->Dt2[0]=Dw2+Db2+2*DbDw;
      popInfoMix->Wt[0]=xWt;
      popInfoMix->d2spred[0]=d2spred;
      popInfoMix->d2sobs[0]=d2s05;

      std::string fichero_sal_NeH = fichero + "_GONE2_Ne_mix";
      std::ofstream salida;
      salida.open(fichero_sal_NeH, std::ios::out);
      salida <<std::fixed<< std::setprecision(7);
      salida <<"# Estimates of Ne using observed LD values (weighted d² correlation coefficients) for different\n";
      salida <<"# recombination rate bins. Extending the approximation of Hayes et al.(2003), which reffers\n";
      salida <<"# to a single population and Ne changing linearly with time, to a model of two populations:\n";
      salida <<"# each Ne estimate from a recombination rate is assumed that corresponds to a particular \n";
      salida <<"# generation in the past. A metapopulation of two subpopulations of equal size is assumed.\n";
      salida <<"# The Fst (subpopulation differentiation index) and the reciprocal migration rate m were\n";
      salida <<"# estimated from the observed d² values of unlinked and weakly linked markers, and by using\n";
      salida <<"# these estimates (Fst="<< popInfo->Fst <<", m="<< popInfo->m <<") the particular Ne and the expected LD\n";
      salida <<"# were calculated for each bin of recrecombination rate.\n";
      // salida <<"# Rec_rate_bin\tgeneration\tNe_metapopulation\td²_wi_subpops\td²_bw_subpops\t2_x_d²_wi_x_bw\td²_predicted_metapop\td²_predicted_sample\td²_observed_sample\n";
      salida <<"# Rec_rate_bin\tgeneration\tNe_metapopulation\td²_predicted_metapop\td²_predicted_sample\td²_observed_sample\n";

      int contador=0;
      for (int nBin=-1;nBin<sampleInfo->binMax;++nBin){
          if (nBin>=0){
            Ne = FixMixEcuacion05link(numIndi, popInfo->f, popInfo->haplotype, popInfo->basecallcorrec, Fst,&popInfo->m,&xWt,&Dw2,&Db2,&DbDw,sampleInfo->xc[nBin],&d2spred,sampleInfo->d2[nBin], KK);
              popInfoMix->xc[contador]=sampleInfo->xc[nBin];
              popInfoMix->Nelink[contador]=2*Ne;
              popInfoMix->mig[contador]=popInfo->m;
              popInfoMix->Fst[contador]=popInfo->Fst;
              popInfoMix->Dw2[contador]=Dw2;
              popInfoMix->Db2[contador]=Db2;
              popInfoMix->DbDw[contador]=DbDw;
              popInfoMix->Dt2[contador]=Dw2+Db2+2*DbDw;
              popInfoMix->Wt[contador]=xWt;
              popInfoMix->d2spred[contador]=d2spred;
              popInfoMix->d2sobs[contador]=sampleInfo->d2[nBin];
          }
          double generacion=1/(2*popInfoMix->xc[contador]);
          if (generacion<155){
            if (popInfoMix->Nelink[contador]<100000000){
              // salida <<std::fixed<< std::setprecision(5)<< popInfoMix->xc[contador]<<"\t"<<std::fixed<< std::setprecision(0)<< generacion<<"\t"<<std::fixed<< std::setprecision(2)<< popInfoMix->Nelink[contador]<< "\t";
              // salida <<std::fixed<< std::setprecision(8)<< Dw2/xWt <<"\t"<< Db2/xWt <<"\t"<< 2*DbDw/xWt <<"\t"<< popInfoMix->Dt2[contador]/xWt <<"\t"<<popInfoMix->d2spred[contador]<<"\t"<<popInfoMix->d2sobs[contador]<<"\n";
              salida <<std::fixed<< std::setprecision(5)<< popInfoMix->xc[contador]<<"\t"<<std::fixed<< std::setprecision(0)<< generacion<<"\t"<<std::fixed<< std::setprecision(2)<< popInfoMix->Nelink[contador]<< "\t";
              salida <<std::fixed<< std::setprecision(8)<< popInfoMix->Dt2[contador]/xWt <<"\t"<<popInfoMix->d2spred[contador]<<"\t"<<popInfoMix->d2sobs[contador]<<"\n";
            }
          }
          ++contador;
      }
      salida.close();
    }
    else{
      std::cerr << "#\n# The Ne estimate does not converge to a solution within the range of Fst from 0 to 1."<<std::endl;
      // exit(EXIT_FAILURE);
    }
  }
  else{
    std::cerr << "# Ne cannot be estimated because the LD is smaller for linked than for unlinked markers.\n"<<std::endl;
  }
  // Free the rest of the memory
  delete[] tnxc;
  delete[] txc;
  delete[] txD2;
  delete[] txW;
  delete[] teffNData;
  delete[] texpNData;
  delete[] valid_idx;
}

//  MIGRACION entre dos subpoblaciones CON LIGAMIENTO:
double MixEcuacionlink(double maxdistance,double mindistance, int indmaxdistanceindx,double* mapdist,double Neini,  double effeneind, double fp, int haplotype, double basecallcorrec, double Fst, double* m, double* Wt,double* Dw2,double* Db2, double* DbDw, double d2sobs, double KK)
{
    int DOSMIL, i, ii, j, k, conta,minj;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increL, increNe, sumad2spred, sumafrec;
    double distancia,fp12, distmin, K12;
    bool flagbreak = false;
    double AA,BB,CC,DD,EE,c,c2,c12,m12,m22;
    double Db2p,Dw2p,DbDwp;

    K12 = (1.0 + KK / 4.0);
    increL=0.0001; // El cromosoma de 10 Morgans (max) se parte en incres de 1/100.000 Morgans
    distmin=mindistance+increL/2;
    minj=int(distmin/increL);

    DOSMIL = 2000;
    fp12 = (1 + fp) * (1 + fp);
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

    double Ne = Neini;
    if (Ne<1000){Ne=1000;}
    if (Ne>1000000){Ne=1000000;}

    *m= (1-Fst)/(16*Ne*Fst);
    if (*m>0.5){*m=0.5;}
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 10; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                sumad2spred = 0;
                distancia = distmin;
                sumafrec = 0;
                j = minj;
                while (distancia < maxdistance)
                {
                  c = (1 - exp(-2 * distancia)) / 2;
                  c2 = c * c;
                  c12 = (1 - c) * (1 - c);
                  m12=(1-2*(*m))*(1-2*(*m));

                  // AA=(1-m12)*m12*(1-1/(2*Ne));
                  // BB=2*(1-m12)*(1-2.2/(4*Ne))*(1-c);
                  // CC=1-(1-2.2/(4*Ne))*c12;
                  // DD=1-m12*(1-1/(2*Ne))*(1-c);
                  // EE=(1-m12)*(1-m12)*(1-c)*(1-2.2/(4*Ne));

                  // *Db2=(*Wt)*Fst*Fst;
                  // *DbDw=(*Db2)*(AA/DD+(EE/CC+(AA*BB)/(CC*DD))/(4*Ne*DD))+(*Wt)*(1-Fst)*(1-Fst)/(16*Ne*Ne*CC*DD);
                  // *Dw2=(*Db2)*(EE/CC+(AA*BB)/(CC*DD))+(*Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*CC);

                  *Db2=(*Wt)*Fst*Fst;
                  *Dw2=(*Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*(1-c12)+2.2*c12);
                  *DbDw=(4*(*Wt)*Fst*Fst*m12*(*m)+(*Dw2)/(4*Ne))/(1-m12*(1-c)*(1-1/2/Ne));

                  Dw2p = (*Dw2)*c12 + (*Db2)*(1-m12)*(1-m12) + (*DbDw)*2*(1-m12)*(1-c);
                  DbDwp = (*DbDw)*m12*(1-c) + (*Db2)*m12*(1-m12);
                  Db2p = (*Db2)*m12*m12;

                  *Dw2=Dw2p;
                  *DbDw=DbDwp;
                  *Db2=Db2p;

                  switch (haplotype){
                    case 0:
                      d2spred=((*Dw2)+4*(*Db2)+4*(*DbDw))/(*Wt)*samplex*basecallcorrec+sampley*fp12;
                      break;
                    case 1:
                      d2spred=((*Dw2)+(*Db2)+2*(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size_h-1);
                      break;
                    case 2:
                      d2spred=((*Dw2)+(*Db2)+2*(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size_h-1);
                      break;
                    case 3:
                      d2spred=((*Dw2)/4+(*Db2)+(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size-1);
                      break;
                  }
                  sumad2spred += *(mapdist+j) * d2spred; // La integral
                  distancia += increL;
                  sumafrec += *(mapdist+j);
                  ++j;
                }
                if (sumafrec > 0)
                {
                    d2spred = (sumad2spred / sumafrec);
                }
                else
                {
                    std::cerr << "Within chromosome integral cannot be computed." << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (d2spred > d2sobs)
                {
                    if (increNe < 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if (d2spred < d2sobs)
                {
                    if (increNe > 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe /= 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                    *m= (1-Fst)/(16*Ne*Fst);
                    if (*m>0.5){*m=0.5;}
                }
                if (abs(increNe) < 0.01)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
    return Ne;
}
//  MIGRACION entre dos subpoblaciones CON c=0.5:
double MixEcuacion05(double Neini,  double effeneind, double fp, int haplotype, double basecallcorrec, double Fst, double* m, double* Wt,double* Dw2,double* Db2, double* DbDw, double d2sobs, double KK)
{
    int MIL, DOSMIL, i, ii, j, k, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double d2spred, increL, ele, increNe, sumad2spred, sumafrec;
    double fp12, K12;
    bool flagbreak = false;
    double AA,BB,CC,DD,EE,c,c2,c12,m12,m22;
    double Db2p,Dw2p,DbDwp;

    K12 = (1.0 + KK / 4.0);
    MIL = 1000;
    DOSMIL = 2000;
    fp12 = (1 + fp) * (1 + fp);
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

    double Ne = Neini;
    if (Ne<1000){Ne=1000;}
    if (Ne>1000000){Ne=1000000;}

    c = 0.5;
    c2 = c * c;
    c12 = (1 - c) * (1 - c);

    *m= (1-Fst)/(16*Ne*Fst);
    if (*m>0.5){*m=0.5;}
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (4 * (ii + 1));
        for (i = 0; i < 20; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {
                m12=(1-2*(*m))*(1-2*(*m));

                // AA=(1-m12)*m12*(1-1/(2*Ne));
                // BB=2*(1-m12)*(1-2.2/(4*Ne))*(1-c);
                // CC=1-(1-2.2/(4*Ne))*c12;
                // DD=1-m12*(1-1/(2*Ne))*(1-c);
                // EE=(1-m12)*(1-m12)*(1-c)*(1-2.2/(4*Ne));

                // *Db2=(*Wt)*Fst*Fst;
                // *DbDw=(*Db2)*(AA/DD+(EE/CC+(AA*BB)/(CC*DD))/(4*Ne*DD))+(*Wt)*(1-Fst)*(1-Fst)/(16*Ne*Ne*CC*DD);
                // *Dw2=(*Db2)*(EE/CC+(AA*BB)/(CC*DD))+(*Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*CC);

                *Db2=(*Wt)*Fst*Fst;
                *Dw2=(*Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*(1-c12)+2.2*c12);
                *DbDw=(4*(*Wt)*Fst*Fst*m12*(*m)+(*Dw2)/(4*Ne))/(1-m12*(1-c)*(1-1/2/Ne));

                Dw2p = (*Dw2)*c12 + (*Db2)*(1-m12)*(1-m12) + (*DbDw)*2*(1-m12)*(1-c);
                DbDwp = (*DbDw)*m12*(1-c) + (*Db2)*m12*(1-m12);
                Db2p = (*Db2)*m12*m12;

                *Dw2=Dw2p;
                *DbDw=DbDwp;
                *Db2=Db2p;

                switch (haplotype){
                  case 0:
                    d2spred=((*Dw2)+4*(*Db2)+4*(*DbDw))/(*Wt)*samplex*basecallcorrec+sampley*fp12;
                    break;
                  case 1:
                    d2spred=((*Dw2)+(*Db2)+2*(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size_h-1);
                    break;
                  case 2:
                    d2spred=((*Dw2)+(*Db2)+2*(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size_h-1);
                    break;
                  case 3:
                    d2spred=((*Dw2)/4+(*Db2)+(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size-1);
                    break;
                }

                if (d2spred > d2sobs)
                {
                    if (increNe < 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if (d2spred < d2sobs)
                {
                    if (increNe > 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe /= 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                    *m= (1-Fst)/(16*Ne*Fst);
                    if (*m>0.5){*m=0.5;}
                }
                if (abs(increNe) < 0.01)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.1)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
    return Ne;
}

//  MIGRACION entre dos subpoblaciones CON c=0.5:
double FixMixEcuacion05link(double effeneind, double fp, int haplotype, double basecallcorrec, double Fst, double* m, double* Wt,double* Dw2,double* Db2, double* DbDw, double c,double* d2spred,double d2sobs, double KK)
{
    int MIL, DOSMIL, i, ii, j, k, conta;
    double sample_size, sample_size_h, samplex, sampley;
    double increL, ele, increNe, sumad2spred, sumafrec;
    double fp12, K12;
    bool flagbreak = false;
    double AA,BB,CC,DD,EE,c2,c12,m12,m22;
    double Db2p,Dw2p,DbDwp;

    K12 = (1.0 + KK / 4.0);
    MIL = 1000;
    DOSMIL = 2000;
    fp12 = (1 + fp) * (1 + fp);
    sample_size = effeneind;
    sample_size_h = sample_size * 2;
    samplex = (sample_size_h - 2.0) * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 8.0 / 5.0 * (sample_size_h - 2.0) * (sample_size_h - 2.0);
    samplex += 4 * (sample_size_h - 2.0);
    samplex /= ((sample_size_h - 1.0) * (sample_size_h - 1.0) * (sample_size_h - 1.0) + (sample_size_h - 1.0) * (sample_size_h - 1.0));
    sampley = (2.0 * sample_size_h - 4.0) / ((sample_size_h - 1.0) * (sample_size_h - 1.0));

    c2 = c * c;
    c12 = (1 - c) * (1 - c);
    m12=(1-2*(*m))*(1-2*(*m));

    double Ne = 5000;
    for (ii = 0; ii < 2; ++ii)
    {
        increNe = Ne / (20 * (ii + 1));
        for (i = 0; i < 20; ++i)
        {
            for (conta = 0; conta < DOSMIL; ++conta)
            {

                // AA=(1-m12)*m12*(1-1/(2*Ne));
                // BB=2*(1-m12)*(1-2.2/(4*Ne))*(1-c);
                // CC=1-(1-2.2/(4*Ne))*c12;
                // DD=1-m12*(1-1/(2*Ne))*(1-c);
                // EE=(1-m12)*(1-m12)*(1-c)*(1-2.2/(4*Ne));

                // *Db2=(*Wt)*Fst*Fst;
                // *DbDw=(*Db2)*(AA/DD+(EE/CC+(AA*BB)/(CC*DD))/(4*Ne*DD))+(*Wt)*(1-Fst)*(1-Fst)/(16*Ne*Ne*CC*DD);
                // *Dw2=(*Db2)*(EE/CC+(AA*BB)/(CC*DD))+(*Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*CC);

                *Db2=(*Wt)*Fst*Fst;
                *Dw2=(*Wt)*(1-Fst)*(1-Fst)*(1+c2)/(4*Ne*(1-c12)+2.2*c12);
                *DbDw=(4*(*Wt)*Fst*Fst*m12*(*m)+(*Dw2)/(4*Ne))/(1-m12*(1-c)*(1-1/2/Ne));

                Dw2p = (*Dw2)*c12 + (*Db2)*(1-m12)*(1-m12) + (*DbDw)*2*(1-m12)*(1-c);
                DbDwp = (*DbDw)*m12*(1-c) + (*Db2)*m12*(1-m12);
                Db2p = (*Db2)*m12*m12;

                *Dw2=Dw2p;
                *DbDw=DbDwp;
                *Db2=Db2p;

                switch (haplotype){
                  case 0:
                    *d2spred=((*Dw2)+4*(*Db2)+4*(*DbDw))/(*Wt)*samplex*basecallcorrec+sampley*fp12;
                    break;
                  case 1:
                    *d2spred=((*Dw2)+(*Db2)+2*(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size_h-1);
                    break;
                  case 2:
                    *d2spred=((*Dw2)+(*Db2)+2*(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size_h-1);
                    break;
                  case 3:
                    *d2spred=((*Dw2)/4+(*Db2)+(*DbDw))/(*Wt)*basecallcorrec+1/(sample_size-1);
                    break;
                }

                if (*d2spred > d2sobs)
                {
                    if (increNe < 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if (*d2spred < d2sobs)
                {
                    if (increNe > 0)
                    {
                        increNe /= -5;
                        break;
                    }
                }
                if ((Ne + increNe) < 3)
                {
                    increNe /= 5;
                    break;
                }
                else
                {
                    Ne += increNe;
                }
                if (abs(increNe) < 0.01)
                {
                    break;
                }
                if (Ne > 100000000)
                {
                    flagbreak = true;
                    break;
                }
                // if (flagbreak){
                //     break;
                // }
            }
            if (abs(increNe) < 0.01)
            {
                break;
            }
            if (flagbreak)
            {
                break;
            }
        }
        if (flagbreak)
        {
            break;
        }
    }
    return Ne;
}