//
//  main.cpp
//  GONE (Genetic Optimization for Ne Estimation)
//
//  Created by Enrique Santiago on 12/11/19.
//  Copyright © 2019 Enrique Santiago. All rights reserved.

#include "libgone.hpp"
#include <iomanip>

void gone(AppParams* params, std::string fichero, int argc, char* argv[], PopulationInfo *popInfo, SampleInfo *sInfo) {
  double clow = 0, chigh = 0.5, kk;
  int j;

  if (params->hc <= params->lc) {
    std::cerr << " Invalid range of recombination frequencies." << std::endl;
    exit(1);
  }

  

  xprng.setSeed(params->semilla);

  std::string fichsal = fichero.substr(0, fichero.length() - 7);

  std::string fichero_sal_NeH = fichsal + "_GONE2_Ne";
  std::string fichero_sal_d2 = fichsal + "_GONE2_d2";
//  std::string fichero_sal_log = fichsal + "_GONE_log";
//  std::string fichero_sal_rearrange = fichsal + "_GONE_input";
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
  
  sampleInfo->mix = params->mix;
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
  for (size_t i = 0; i < linesRead; ++i) {
    sInfo->nxc[i] = static_cast<long int>(sampleInfo->nBin[i]);
  }
  memcpy(&(sInfo->d2[0]), &(sampleInfo->d2cObs[0]), linesRead * sizeof(double));
  sInfo->binMax = linesRead;
  //  std::remove(fichero.c_str());
  //  std::cout << "LinesRead:"<<linesRead<<"\n"; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  double sample_size_h = popInfo->avgNumIndiAnalyzed;
  sampleInfo->cValMin = MAX_DOUBLE;

  if ((sampleInfo->haplotype == 0)) { // OJO AQUI: PROVISIONAL
      double sample_size_h = popInfo->avgNumIndiAnalyzed * 2.0;
      // sampleInfo->fVal = (1.0 + sampleInfo->fValSample * (sample_size_h - 1.0)) /
      //               (sample_size_h - 1 + sampleInfo->fValSample);
      // std::cout<<sampleInfo->fVal<<std::endl;
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

  //  salida.open(fichero_sal_log, std::ios::out);
  //  salida << "(GONE v2)\n\n";
  //  salida << "Input File: " << fichero << "\n\n";
  //  salida << "Command:";
  //  for (i = 0; i < argc; ++i) {
  //    salida << " " << argv[i];
  //  }
  //  salida << "\n\n";
  //  salida.close();

  if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
    salida.open(fichero_sal_evol, std::ios::out);
    salida << "Gener\tSCbest\tSCmed1\tnsegbest\tnsegmed\n";
    salida.close();
  }

  // Lo que sigue es para una sola Poblacion panmíctica
  if (!params->mix){

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

    #ifdef LATEST_OMP
    // Main genetic loop
    #pragma omp parallel for simd reduction(+ : avgSCVal, avgNe, avgD2Pred)
    for (int _i = 0; _i < GONE_ROUNDS; ++_i) {
      Pool* privpool = new Pool();
      Bicho* bestBicho;
      SetInitialPoolParameters(privpool, sampleInfo->cValMin);
      PrePopulatePool(privpool, sampleInfo);
      if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
        RunDbg(privpool, sampleInfo, fichsal);
      } else {
        Run(privpool, sampleInfo, fichsal);
      }
      bestBicho = &pool->parents[0];
      for (int i = 0; i < sampleInfo->nBins; ++i) {
        avgD2Pred[i] += bestBicho->d2cPred[i] / GONE_ROUNDS;
      }
      avgSCVal += bestBicho->SCval / GONE_ROUNDS;


      int conta = 0;
      for (int i = 0; i < bestBicho->nSeg; ++i) {
        for (int j = bestBicho->segBl[i]; j < bestBicho->segBl[i + 1]; ++j) {
          conta = conta + 1;
        }
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
        avgNe[i] += sumNe[i] / GONE_ROUNDS;
      }
      delete privpool;
    }
    #else
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
    Bicho* Bicho0;
    Bicho* Bicho1;
    Bicho* Bicho2;
    double SCsube, SCbaja, SCant, incresube, increbaja, incre, testNe1, testNe2;

    for (int _i = tid; _i < GONE_ROUNDS; _i+=numThreads) {
      *privpool = emptyPool;
      int ii=0;
      int jj=0;
      int pp=0;
      SetInitialPoolParameters(privpool, sampleInfo->cValMin);
      PrePopulatePool(privpool, sampleInfo);
      if ((sampleInfo->flags & FLAG_DEBUG) > 0) {
        RunDbg(privpool, sampleInfo, fichsal);
      } else {
        Run(privpool, sampleInfo, fichsal);
      }

      // // POST-AJUSTE DE LOS params->muestraSalida BICHOS MEJORES:
      // int basa = sampleInfo->muestraSalida;
      //   int ancho;
      //   int ancho1;
      //   int ancho2;
      //   double max1,max2,min1,min2;
      //   int increii=0;// ***
      //   int increlateral=0;// ***
      //   int desdeii=0;// ***
      //   int hastaii=0;// ***
      //   double sumHar=0;
      //   double sumGeo=1;
      // for (jj=0; jj<sampleInfo->muestraSalida; ++jj){
      //   Bicho0 = &privpool->parents[jj]; // Apunta al original
      //   Bicho1 = &privpool->parents[basa+0]; // Extendida en número de segmentos
      //   Bicho2 = &privpool->parents[basa+1]; // Extendida en número de segmentos (para operaciones)
      //   Bicho0->SCval = CalculaSC(Bicho0, sampleInfo);
      //   // Extiende Bicho1:
      //   int conta3=0;
      //   for (ii = 0; ii < Bicho0->nSeg; ++ii) {
      //     Bicho1->NeBl[conta3]=Bicho0->NeBl[ii];
      //     Bicho1->segBl[conta3]=Bicho0->segBl[ii];
      //     ++conta3;
      //     int ancho=Bicho0->segBl[ii+1]-Bicho0->segBl[ii];
      //     if (((ancho>3) && (Bicho0->segBl[ii]<10)) || ((ancho>5) && (Bicho0->segBl[ii]<50)) || ((ancho>8) && (Bicho0->segBl[ii]<80))) {
      //       ancho1=int(ancho/2);
      //       ancho2=ancho-ancho1;
      //       Bicho1->NeBl[conta3]=Bicho0->NeBl[ii];
      //       Bicho1->segBl[conta3]=Bicho0->segBl[ii]+ancho1;
      //       ++conta3;
      //     }
      //   }
      //   Bicho1->nSeg=conta3;
      //   Bicho1->segBl[conta3]=Bicho0->segBl[ii];
      //   Bicho1->SCval = CalculaSC(Bicho1, sampleInfo);
      //   Bicho1->efval=Bicho0->efval;

      //   // Copia Bicho1 en Bicho2
      //   for (ii = 0; ii <= Bicho1->nSeg;++ii){
      //     Bicho2->NeBl[ii]=Bicho1->NeBl[ii];
      //     Bicho2->segBl[ii]=Bicho1->segBl[ii];
      //   }
      //   Bicho2->nSeg=Bicho1->nSeg;
      //   Bicho2->efval=Bicho1->efval;
      //   Bicho2->SCval = Bicho1->SCval;

      //   // SUBE Y BAJA TODOS LOS SEGMENTOS POR PARES ADYACENTES:
      //   int alternaii=1000;// ***
      //   for (int kk=0; kk<30000;++kk){ // Repeticiones de afinamiento creciente
      //     ++alternaii;
      //     if (alternaii > 0){alternaii=0;}
      //     if (alternaii == 0){// *** Aleatorio + 1
      //       increlateral = 1;
      //       increii = 1;
      //       desdeii = static_cast<int>(xprng.uniforme01() * (Bicho1->nSeg - 1));
      //       // (--desdeii);
      //       // if (desdeii<0){desdeii=0;}
      //       hastaii = desdeii + 1; // Este ya no lo hace
      //     }
      //     else if (alternaii == 1){ // *** Aleatorio + 2
      //       increlateral = 2;
      //       increii = 1;
      //       desdeii = static_cast<int>(xprng.uniforme01() * (Bicho1->nSeg - 2));
      //       (--desdeii);if (desdeii<0){desdeii=0;}
      //       hastaii = desdeii + 1; // Este ya no lo hace
      //     }
      //     else{// *** Aleatorio + 3
      //       increlateral = 3;
      //       increii = 1;
      //       desdeii = static_cast<int>(xprng.uniforme01() * (Bicho1->nSeg - 3));
      //       (--desdeii);if (desdeii<0){desdeii=0;}
      //       hastaii = desdeii + 1; // Este ya no lo hace
      //     }

      //     // else if (alternaii <= 10){// *** El primero y segundo con los 10 siguientes
      //     //   increlateral = alternaii+1;
      //     //   increii = 1;
      //     //   desdeii = 0;
      //     //   hastaii = 20; // Este ya no lo hace
      //     // }
      //     // else if (alternaii == 11){// *** Hacia alante
      //     //   increlateral = 1;
      //     //   increii = 1;
      //     //   desdeii = 0;
      //     //   hastaii = Bicho1->nSeg - 1; // Este ya no lo hace
      //     // }
      //     // else if (alternaii <= 21){// *** El primero y segundo con los 10 siguientes
      //     //   increlateral = alternaii-10;
      //     //   increii = 1;
      //     //   desdeii = 0;
      //     //   hastaii = 20; // Este ya no lo hace
      //     // }
      //     // else {// *** Hacia atras
      //     //   increlateral = -1;
      //     //   increii = -1;
      //     //   desdeii = std::min(20,Bicho1->nSeg - 1);
      //     //   hastaii = 0; // Este ya no lo hace
      //     // }
      //     // SIGUEN POSIBILIDADES DE FOR: ADELANTE, ATRAS o ALTERNADS
      //     for (ii = desdeii; ii != hastaii; ii += increii){
      //       if (((ii+increlateral)<(Bicho1->nSeg - 1)) && ((ii+increlateral)>=0)){ 
      //         // Limites para el de referencia
      //         max1 = Bicho1->NeBl[ii] * 2; // Doble o mitad
      //         min1 = Bicho1->NeBl[ii] / 2;
      //         if (min1<5){min1=5;}

      //         // Limites para el otro
      //         max2 = Bicho1->NeBl[ii+increlateral] * 2; // Doble o mitad
      //         min2 = Bicho1->NeBl[ii+increlateral] / 2;  
      //         if (min2<5){min2=5;}

      //         // sumHar = (1/Bicho1->NeBl[ii]+1/Bicho1->NeBl[ii+increlateral]); // Suma de inversas para harmonica
      //         sumGeo = Bicho1->NeBl[ii] * Bicho1->NeBl[ii+increlateral]; // Suma de inversas para harmonica
              
      //         //double prodGeo = Bicho1->NeBl[ii]*Bicho1->NeBl[ii+increlateral]; // Producto para media geometrica

      //         // maneja la de referencia para subir o bajar
      //         // primero prueba el efecto de subir un solo incremento
      //         incresube= (max1 - Bicho1->NeBl[ii]) / 10; // Incrementos del 5%
      //         testNe1 = Bicho1->NeBl[ii]+incresube;
      //         // testNe2=1/(sumHar-1/testNe1); // Harmonica
      //         testNe2=sumGeo/testNe1; // Geometrica
      //         if (testNe2>min2){
      //           Bicho2->NeBl[ii]=testNe1;
      //           Bicho2->NeBl[ii+increlateral]=testNe2;
      //           SCsube = CalculaSC(Bicho2, sampleInfo);
      //         }
      //         else{
      //           SCsube = std::numeric_limits<double>::max();
      //         }
      //         // luego prueba el efecto de bajar un solo incremento
      //         increbaja= (min1-Bicho1->NeBl[ii]) / 20; // Incrementos del 5% (negativos)
      //         testNe1 = Bicho1->NeBl[ii]+increbaja;
      //         // testNe2=1/(sumHar-1/testNe1); // Harmonica
      //         testNe2=sumGeo/testNe1; // Geometrica
      //         while (testNe2<5){
      //           increbaja/=2;
      //           testNe1 = Bicho1->NeBl[ii]+increbaja;
      //           // testNe2=1/(sumHar-1/testNe1); // Harmonica
      //           testNe2=sumGeo/testNe1; // Geometrica
      //         }
      //         if (testNe2<max2){
      //           Bicho2->NeBl[ii]=testNe1;
      //           Bicho2->NeBl[ii+increlateral]=testNe2;
      //           SCbaja = CalculaSC(Bicho2, sampleInfo);
      //         }
      //         else{
      //           SCbaja = std::numeric_limits<double>::max();
      //         }
      //         if (SCsube<SCbaja){
      //           if (SCsube<Bicho1->SCval){
      //             incre=incresube;
      //           }
      //           else{
      //             incre=0;
      //           }
      //         }
      //         else if (SCsube>SCbaja){
      //           if (SCbaja<Bicho1->SCval){
      //             incre=increbaja;
      //           }
      //           else{
      //             incre=0;
      //           }
      //         }
      //         else{
      //           incre=0;
      //         }
      //         Bicho2->NeBl[ii] = Bicho1->NeBl[ii];
      //         Bicho2->NeBl[ii+increlateral] = Bicho1->NeBl[ii+increlateral];
      //         if (incre!=0){
      //           if (incre>0){ // A subir
      //             Bicho2->NeBl[ii] += incre;
      //             // Bicho2->NeBl[ii+increlateral] = 1/(sumHar-1/Bicho2->NeBl[ii]); // Harmonica
      //             Bicho2->NeBl[ii+increlateral] = sumGeo/Bicho2->NeBl[ii]; // Geometrica
      //             SCant = CalculaSC(Bicho2, sampleInfo);
      //             while ((Bicho2->NeBl[ii]<max1) && (Bicho2->NeBl[ii+increlateral]>min2)){
      //                 Bicho2->NeBl[ii]+=incre;
      //                 // Bicho2->NeBl[ii+increlateral] = 1/(sumHar-1/Bicho2->NeBl[ii]); // Harmonica
      //                 Bicho2->NeBl[ii+increlateral] = sumGeo/Bicho2->NeBl[ii]; // Geometrica
      //                 SCsube = CalculaSC(Bicho2, sampleInfo);
      //                 if (SCsube > SCant){
      //                   break;
      //                 }
      //                 SCant=SCsube;
      //             }
      //             testNe1 = (Bicho2->NeBl[ii] - incre + 1 * Bicho1->NeBl[ii]) / 2; // AJUSTE SUAVE (1/5 del mejor)
      //             // testNe2 = 1/(sumHar-1/testNe1); // Harmonica
      //             testNe2 = sumGeo/testNe1; // Geometrica
      //             Bicho1->NeBl[ii] = testNe1;
      //             Bicho1->NeBl[ii+increlateral] = testNe2;
      //             Bicho1->SCval = CalculaSC(Bicho1, sampleInfo);
      //           } 
      //           else{ // A bajar
      //             Bicho2->NeBl[ii] += incre;
      //             // Bicho2->NeBl[ii+increlateral] = 1/(sumHar-1/Bicho2->NeBl[ii]); // Harmonica
      //             Bicho2->NeBl[ii+increlateral] = sumGeo/Bicho2->NeBl[ii]; // Geometrica
      //             while (Bicho2->NeBl[ii+increlateral]<5){
      //               incre/=2;
      //               Bicho2->NeBl[ii] += incre;
      //               // Bicho2->NeBl[ii+increlateral] = 1/(sumHar-1/Bicho2->NeBl[ii]); // Harmonica
      //               Bicho2->NeBl[ii+increlateral] = sumGeo/Bicho2->NeBl[ii]; // Geometrica
      //             }
      //             SCant = CalculaSC(Bicho2, sampleInfo);
      //             while ((Bicho2->NeBl[ii]>min1) && (Bicho2->NeBl[ii+increlateral]<max2) && (Bicho2->NeBl[ii+increlateral]>5)){
      //                 Bicho2->NeBl[ii]+=incre;
      //                 // Bicho2->NeBl[ii+increlateral] = 1/(sumHar-1/Bicho2->NeBl[ii]); // Harmonica
      //                 Bicho2->NeBl[ii+increlateral] = sumGeo/Bicho2->NeBl[ii]; // Geometrica
      //                 SCbaja = CalculaSC(Bicho2, sampleInfo);
      //                 if (SCbaja > SCant){
      //                   break;
      //                 }
      //                 SCant=SCbaja;
      //             }
      //             testNe1 = (Bicho2->NeBl[ii] - incre + 1 * Bicho1->NeBl[ii]) / 2; // AJUSTE SUAVE (1/5 del mejor)
      //             // testNe2 = 1/(sumHar-1/testNe1); // Harmonica
      //             testNe2 = sumGeo/testNe1; // Geometrica
      //             Bicho1->NeBl[ii] = testNe1;
      //             Bicho1->NeBl[ii+increlateral] = testNe2;
      //             Bicho1->SCval = CalculaSC(Bicho1, sampleInfo);
      //           }
      //         }
      //         // Copia completa de Bicho1 en Bicho2
      //         for (pp = 0; pp <= Bicho1->nSeg;++pp){
      //           Bicho2->NeBl[pp]=Bicho1->NeBl[pp];
      //           Bicho2->segBl[pp]=Bicho1->segBl[pp];
      //         }
      //         Bicho2->nSeg=Bicho1->nSeg;
      //         Bicho2->efval=Bicho1->efval;
      //         Bicho2->SCval = Bicho1->SCval;

      //         // Copia de Bicho1 a Bicho2 solo lo que cambia
      //         // Bicho2->NeBl[ii] = Bicho1->NeBl[ii];
      //         // Bicho2->NeBl[ii+increlateral] = Bicho1->NeBl[ii+increlateral];
      //         // Bicho2->SCval = Bicho1->SCval;

      //       }
      //     }
      //   }

      //   privpool->parents[jj] = *Bicho1; // ACTUALIZA EL BICHO.SI SE ANULA ESTA LINEA, SE ANULA EL POST-TRATAMIENTO

      // }
      // // FIN DEL POSTAJUSTE


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
    #endif
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
        salida << sampleInfo->cVal[i] << "\t" << sampleInfo->nBin[i] << "\t" << sampleInfo->d2cObs[i] << "\t"
            << avgD2Pred[i] << "\n";
      }
    }
    salida.close();

    delete pool;
    delete sampleInfo;

    // if (!params->quiet) {
    //   std::cout << " Ne estimation saved to " << fichero_sal_NeH << std::endl<< std::endl;
    // }
  }
}
