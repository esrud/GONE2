//
//  main.cpp
//  GONE (Genetic Optimization for Ne Estimation)
//
//  Created by Enrique Santiago on 12/11/19.
//  Copyright © 2019 Enrique Santiago. All rights reserved.

#include "libgone.hpp"

void gone(AppParams* params, std::string fichero, int argc, char* argv[], PopulationInfo *popInfo, SampleInfo *sInfo) {
  /*
   * Main genetic algorithm loop
  */
  double clow = 0, chigh = 0.5;
  int j;

  if (params->hc <= params->lc) {
    std::cerr << " Invalid range of recombination frequencies." << std::endl;
    exit(1);
  }

  xprng.setSeed(params->semilla);

  std::string fichsal = fichero.substr(0, fichero.length() - 7);

  std::string fichero_sal_NeH = fichsal + "_GONE2_Ne";
  std::string fichero_sal_d2 = fichsal + "_GONE2_d2";
  std::string fichero_sal_evol = fichsal + "_GONE2_evol";  // --debug
  std::string fichero_dbg = fichsal + "_GONE2_dbg";      // --debug
  std::ofstream salida;

  GsampleInfo* sampleInfo = new GsampleInfo();
  sampleInfo->flags = params->flags;
  sampleInfo->hetEspAll = popInfo->hetEspAll;
  sampleInfo->hetEsp = popInfo->hetEsp;
  sampleInfo->hetAvgAll = popInfo->hetAvgAll;
  sampleInfo->hetAvg = popInfo->hetAvg;
  
  if (params->basecallerror == 0){
      sampleInfo->basecallcorrec = 1;
  }
  else{
      sampleInfo->basecallcorrec = pow((1 - 4* params->basecallerror * (1- params->basecallerror)),2);
  }
  popInfo->basecallcorrec=sampleInfo->basecallcorrec;
  
  sampleInfo->mix = sInfo->mix;
  // Set args based on flags
  if ((sampleInfo->flags & FLAG_RESIZE_BINS) > 0) {
    sampleInfo->sizeBins = params->sizeBins;
  } else {
    sampleInfo->sizeBins = DEFAULT_BIN_SIZE;
  }
  sampleInfo->nBins = kNumBins;

  int linesRead = ProcessFile(fichero, clow, chigh, sampleInfo);
  if (linesRead < 15) {
    std::cerr << "There are not enough recombination bins to perform the analysis.\n";
    exit(1);
  }
  // Copy values needed for mix calculation
  memcpy(&(sInfo->xc[0]), &(sampleInfo->cVal[0]), linesRead * sizeof(double));
  for (int i = 0; i < linesRead; ++i) {
    sInfo->nxc[i] = static_cast<long int>(sampleInfo->nBin[i]);
  }
  memcpy(&(sInfo->d2[0]), &(sampleInfo->d2cObs[0]), linesRead * sizeof(double));
  sInfo->binMax = linesRead;

  sampleInfo->cValMin = MAX_DOUBLE;

  if (sampleInfo->haplotype == 0) { // OJO AQUI: PROVISIONAL
      sampleInfo->correccion = 1.0 / Square<double>((1.0 + sampleInfo->fVal));
  }

  sampleInfo->muestraSalida = params->muestraSalida;

  CalculateSumNBins(sampleInfo);
  CalculateAverageNe(sampleInfo);

  // Correccion muestreo generaciones sucesivas:
  double g=params->ngensampling;
  for (int i = 0; i < sampleInfo->nBins; ++i) {
    int indx = sampleInfo->indx[i];
    double ac=1-sampleInfo->cVal[indx];
    if (g==1){
      sampleInfo->ngensamplingcorrec[indx]=1;
    }
    else{
      sampleInfo->ngensamplingcorrec[indx]=(g-2*ac+2*pow(ac,g+1)-ac*ac*g)/(g*g*(ac-1)*(ac-1));
    }
  }  

  if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
    salida.open(fichero_sal_evol, std::ios::out);
    salida << "Gener\tSCbest\tSCmed1\tnsegbest\tnsegmed\n";
    salida.close();
  }

  // Lo que sigue es para una sola Poblacion panmíctica
  if ((!params->mix)){

    // COMIENZAN LOS CICLOS
    Pool* pool = new Pool();
    SetInitialPoolParameters(pool, sampleInfo->cValMin);
    PrePopulatePool(pool, sampleInfo);

    if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
      salida.open(fichero_dbg, std::ios::app);
      for (j = 0; j < pool->parents[0].nSeg; ++j) {
        salida << j << "\t" << pool->parents[0].segBl[j] << "\t"
              << pool->parents[0].NeBl[j] / 2.0 << "\n";
      }
      salida << "\n";
      salida.close();
    }
    double avgD2Pred[kNumLinMax] = {};
    double avgSCVal = 0.0;
    double avgNe[MAXBINS] = {};

    static const Pool emptyPool = Pool();
    int numThreads = params->numThreads;
    if (numThreads == 0) {
      numThreads = omp_get_max_threads();
    }
    double **tavgNe = new double*[numThreads];
    double **tavgD2Pred = new double*[numThreads];
    double *tavgSCVal = new double[numThreads]{};
    int *maxNeConta = new int[numThreads]{};
    // Allocate mamory for avgNe and avgD2Pred
    for (int z = 0; z < numThreads; ++z) {
      tavgD2Pred[z] = new double[kNumLinMax]{};
      tavgNe[z] = new double[kNumLinMax]{};
      std::fill_n(tavgNe[z], MAXBINS, 1.0);
    }
    params->progress.InitCurrentTask((float)GONE_ROUNDS);
    #pragma omp parallel
    {
    int tid = omp_get_thread_num();
    Pool* privpool = new Pool();
    Bicho* bestBicho;

    for (int _i = tid; _i < GONE_ROUNDS; _i+=numThreads) {
      *privpool = emptyPool;
      SetInitialPoolParameters(privpool, sampleInfo->cValMin);
      PrePopulatePool(privpool, sampleInfo);
      if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
        RunDbg(privpool, sampleInfo, fichsal);
      } else {
        Run(privpool, sampleInfo, fichsal);
      }

      bestBicho = &privpool->parents[0];
      for (int i = 0; i < sampleInfo->nBins; ++i) {
        tavgD2Pred[tid][i] += bestBicho->d2cPred[i] / GONE_ROUNDS;
      }
      tavgSCVal[tid] += bestBicho->SCval / GONE_ROUNDS;

      int conta = 0;
      for (int i = 0; i < bestBicho->nSeg; ++i) {
        for (int j = bestBicho->segBl[i]; j < bestBicho->segBl[i + 1]; ++j) {
          conta = conta + 1;
        }
      }
      if (conta > kNumLinMax) {
        conta = kNumLinMax;
      }
      int conta2 = conta;
      double sumNe[kNumGenMax] = {0};

      // MEDIA GEOMETRICA DE LOS muestrasalida PRIMEROS
      for (int i = 0; i < conta; ++i) {
        sumNe[i] = 1;
      }  //
      for (int ii = 0; ii < sampleInfo->muestraSalida; ++ii) {
        conta = 0;
        for (int i = 0; i < privpool->parents[ii].nSeg; ++i) {
          for (int j = privpool->parents[ii].segBl[i];
              j < privpool->parents[ii].segBl[i + 1]; ++j) {
            sumNe[conta] *= pow(privpool->parents[ii].NeBl[i],
                                1.0 / sampleInfo->muestraSalida);
            conta = conta + 1;
            if (conta > privpool->poolParams.gmax[0]) {
              break;
            }
          }
          if (conta > privpool->poolParams.gmax[0]) {
            break;
          }
        }
      }

      for (int i = 0; i < conta2; ++i) {
        sumNe[i] /= 2;
        tavgNe[tid][i] *= pow(sumNe[i], 1.0 / GONE_ROUNDS);
        //tavgNe[tid][i] += sumNe[i] / GONE_ROUNDS;
      }
      if (conta2 > maxNeConta[tid]) {
        maxNeConta[tid] = conta2;
      }
      if (tid == 0) {
        params->progress.SetTaskProgress(_i+1);
        //params->progress.PrintProgress();
      }
    }
    delete privpool;
    }
    // Reduction
    std::fill_n(avgNe, MAXBINS, 1.0);
    for (int nt = 0; nt < numThreads; ++nt) {
      // Reduce avgNe
      for (int i = 0; i < maxNeConta[nt]; ++i) {
        avgNe[i] *= tavgNe[nt][i];
      }
      // Reduce avgD2Pred
      for (int i = 0; i < sampleInfo->nBins; ++i) {
        avgD2Pred[i] += tavgD2Pred[nt][i];
      }
      avgSCVal += tavgSCVal[nt];
      delete[] tavgNe[nt];
      delete[] tavgD2Pred[nt];
    }
    // Free the memory
    delete[] tavgNe;
    delete[] tavgD2Pred;
    delete[] tavgSCVal;
    delete[] maxNeConta;
    // Write output
    salida.open(fichero_sal_NeH, std::ios::out);
    int generacion;
    if (params->haplotype == 1){
      salida <<"Generation\tNe_haploids\n";
      for (int j = 0; j < pool->poolParams.gmax[2] + 1; ++j) {
        generacion = j+1;
        if (generacion<151){
          avgNe[j] *= 2;
          salida << generacion << "\t" << std::max(avgNe[j], MIN_NE_SIZE) << "\n";
        }
      }
    }
    else{
      salida <<"Generation\tNe_diploids\n";
      for (int j = 0; j < pool->poolParams.gmax[2] + 1; ++j) {
        generacion = j+1;
        if (generacion<151){
          salida << generacion << "\t" << std::max(avgNe[j], MIN_NE_SIZE) << "\n";
        }
      }
    }
    salida.close();

    salida.open(fichero_sal_d2, std::ios::out);
    salida <<"c_bin\tnumber_of_SNP_pairs\tObserved_d2\tPredicted_d2\n";

    for (int i = 0; i < linesRead; ++i) {
      if (sampleInfo->cVal[i] != 0) {
        salida << std::fixed << std::setprecision(8) << sampleInfo->cVal[i] << "\t"
          << std::fixed << std::setprecision(0) << sampleInfo->nBin[i]
          << "\t" << std::fixed << std::setprecision(8) << sampleInfo->d2cObs[i] << "\t"
          << avgD2Pred[i] << "\n";
      }
    }
    salida.close();

    delete pool;
    delete sampleInfo;
  }
}
