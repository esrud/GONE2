#pragma once

#include "./gsample.hpp"
#include "./constants.hpp"

typedef struct {
  double efval;
  double SCval;
  int nSeg;
  int segBl[kNumLinMax];
  double NeBl[kNumLinMax];
  double d2cPred[kNumLinMax];
} Bicho;


double CalculaSC(Bicho* bicho, GsampleInfo* sampleInfo);
