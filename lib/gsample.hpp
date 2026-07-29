#pragma once

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <random>

#include "./constants.hpp"
#include "./ga_config.hpp"
#include "./simplemath.hpp"
#include "rng/Xoshiro256plus.h"

typedef struct {
  int haplotype;
  bool mix;
  int binExtra;
  double coveragecorrec;
  double ngensamplingcorrec[kNumLinMax];
  double basecallcorrec;
  double hetEsp;
  double hetEspAll;
  double hetAvg;
  double hetAvgAll;
  double sampleSize;
  double fValSample;
  double fVal;
  double nBin[kNumLinMax];
  double cVal[kNumLinMax];
  double oneMinuscValSq[kNumLinMax];
  double oneMinuscValSqPow[kNumLinMax][kNumGenMax];
  double onePluscValSq[kNumLinMax];
  double cValSq[kNumLinMax];
  double cValRep[kNumLinMax];
  double d2cObs[kNumLinMax];
  int indx[kNumLinMax];
  long int nBins;
  double sumNBin;
  double cValMin;
  double cValMax;
  double sampleX;
  double sampleY;
  double sampleZ1;
  double sampleZ2;
  double sampleZ3;
  double sampleZ4;
  double sampleZ5;
  double correccion;
  bool structuredMix;
  double structScale[kNumLinMax];
  double structOffset[kNumLinMax];
  long int sizeBins;
  int muestraSalida;
  bool hayrecentbins;
  double NeMed;
  uint8_t flags;
  GAConfig gaConfig;
} GsampleInfo;

int ProcessFile(std::string fichero, double clow, double chigh,
                GsampleInfo* sampleInfo);

void CalculateNBins(GsampleInfo* sampleInfo, const int linesRead);
int CompressSample(GsampleInfo* sampleInfo, int linesRead);
void CalculateSumNBins(GsampleInfo* sampleInfo);
void CalculateAverageNe(GsampleInfo* sampleInfo);

extern thread_local Xoshiro256plus xprng;
