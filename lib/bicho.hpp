#pragma once

#include "./gsample.hpp"
#include "./constants.hpp"

typedef struct {
  double efval;
  double SCval;
  int nSeg;
  int segBl[kNumLinMax];
  double NeBl[kNumLinMax];
} Bicho;


double CalculaSC(Bicho* bicho, GsampleInfo* sampleInfo, double* d2cPred);
double CalculaSCScoreOnly(Bicho* bicho, GsampleInfo* sampleInfo);
double CalculaSCScoreCutoff(Bicho* bicho, GsampleInfo* sampleInfo,
                            double cutoff);
