#pragma once

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>

#include "population.hpp"
#include "constants.hpp"
#include "simplemath.hpp"
#include "console.hpp"

typedef struct {
  int haplotype;
  bool mix;
  int sampleSizeIndi;
  int numLoci;
  int numSegLoci;
  long int numSNPs;
  int binMax;
  int binExtra;
  double numIndiEff;
  double increW;
  double f;
  double hetpq;
  double hetAvg;
  double hetEspAll;
  double hetEsp;
  double hetAvgAll;
  double hetVar;
  double hetSesg;
  double parentAvg;
  double frec[MAXLOCI];
  double homo[MAXLOCI];
  double het[MAXIND];
  double parent[MAXIND];
  bool segrega[MAXLOCI];
  long int nxc[MAXBINS];
  double d2[MAXBINS];
  double xc[MAXBINS];
  bool hayrecentbins;
  int shuffledIndi[MAXIND];
  int shuffledLoci[MAXLOCI];
} SampleInfo;

typedef struct{
  double xc[MAXBINS];
  double ene[MAXBINS];
  double gen[MAXBINS];
  double Nelink[MAXBINS];
  double mig[MAXBINS];
  double Fst[MAXBINS];
  double Dw2[MAXBINS];
  double Db2[MAXBINS];
  double DbDw[MAXBINS];
  double Dt2[MAXBINS];
  double Wt[MAXBINS];
  double d2sobs[MAXBINS];
  double d2spred[MAXBINS];
}PopulationInfoMix;

void CalculateFrequencies(PopulationInfo* popInfo, SampleInfo* sampleInfo);
void CalculateF(SampleInfo* sampleInfo, PopulationInfo* popInfo);
bool CalculateHeterozigosity(PopulationInfo* popInfo, SampleInfo* sampleInfo);
bool CalculateKinshipDistribution(PopulationInfo* popInfo,
                                  SampleInfo* sampleInfo);

void CalculateD2Parallel(PopulationInfo* popInfo, SampleInfo* sampleInfo,
                 AppParams* params);
void CalculateD2ParallelFst(std::string fichero,PopulationInfoMix* popInfoMix, PopulationInfo* popInfo, SampleInfo* sampleInfo,
                 AppParams* params);
void CalculateD2GPU(PopulationInfo* popInfo, SampleInfo* sampleInfo,
                 AppParams* params);
double MixEcuacion05(double Neini,double numIndi, double efe, int haplotype, double basecallcorrec, double Fst, double* m,double* xWt,double* Dw2,double* Db2,double* DbDw,double d2s05, double KK);                 
double MixEcuacionlink(double maxdistance, double mindistance, int indmaxdistanceindx,double* mapdist,double Neini,double numIndi, double efe, int haplotype, double basecallcorrec, double Fst, double* m,double* xWt,double* Dw2,double* Db2,double* DbDw,double d2s05, double KK);                 
double FixMixEcuacion05link(double numIndi, double efe, int haplotype, double basecallcorrec, double Fst, double* m, double* Wt,double* Dw2,double* Db2, double* DbDw, double cc,double* d2spred,double d2sobs, double KK);
