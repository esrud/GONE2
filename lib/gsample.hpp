#pragma once

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <random>

#include "./constants.hpp"
#include "./simplemath.hpp"
#include "rng/Xoshiro256plus.h"

typedef struct {
  int haplotype;
  bool mix;
  int binExtra;
  double coveragecorrec;
  double ngensamplingcorrec[MAXIND];
  double basecallcorrec;
  double hetEsp;
  double hetEspAll;
  double hetAvg;
  double hetAvgAll;
  double sampleSize;
  double fValSample;
  double fVal;
  double nBin[MAXIND];
  double cVal[MAXIND];
  double oneMinuscValSq[MAXIND];
  double onePluscValSq[MAXIND];
  double cValSq[MAXIND];
  double cValRep[MAXIND];
  double d2cObs[MAXIND];
  int indx[MAXIND];
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
  long int sizeBins;
  int muestraSalida;
  bool hayrecentbins;
  double NeMed;
  uint8_t flags;
} GsampleInfo;

int ProcessFile(std::string fichero, double clow, double chigh,
                GsampleInfo* sampleInfo);

void CalculateNBins(GsampleInfo* sampleInfo, const int linesRead);
int CompressSample(GsampleInfo* sampleInfo, int linesRead);
void CalculateSumNBins(GsampleInfo* sampleInfo);
void CalculateAverageNe(GsampleInfo* sampleInfo);

//std::random_device seed;
// mersenne_twister_engine 64bit (very good but very big)
//extern std::mt19937_64 rng;
// uses the result of the engine to generate uniform dist
//extern std::uniform_real_distribution<> uniforme01;
extern Xoshiro256plus xprng;
