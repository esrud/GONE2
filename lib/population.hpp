#pragma once

#include <stdint.h>
#include <iostream>

#include "constants.hpp"

typedef struct {
  int haplotype;
  bool mix;
  int numIndi;
  int numLoci;
  int numCromo;
  double avgNumIndiAnalyzed;
  double hetEspAll;
  double hetEsp;
  double hetAvgAll;
  double hetAvg;
  double expNData;
  double effNData;
  double f;
  double Fst;
  double basecallcorrec;
  double currentNe;
  double m;
  double propMiss;
  double Mtot;
  double Mbtot;
  double chrsize[MAXCROMO];
  double posiCM[MAXLOCI];
  int posiBP[MAXLOCI];
  int cromo[MAXLOCI];
  bool hayrecentbins;
  uint8_t indi[MAXIND][MAXLOCI];
} PopulationInfo;

