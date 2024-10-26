#pragma once

#include <algorithm>
#include <string>
#include "./bicho.hpp"

typedef struct {
  int gmax[3];
  int topePosInicio[4];
  double fval;
} PoolParams;

typedef struct {
  double frecMut;
  double efectoMut;
  double efectoMutLateral;
  double frecRec;
  double frecNoRec;
  double frecMutLateral;
  int maxSegmentos;
  bool solape;
} MutParams;

typedef struct {
  Bicho parents[MAX_POPULATION];
  Bicho children[MAX_POPULATION];
  PoolParams poolParams;
  MutParams mutParams;
  int numGenerations;
  int numChildren;
} Pool;

void SetInitialPoolParameters(Pool* pool, double mincval);
void PrePopulatePool(Pool* pool, GsampleInfo* sampleInfo);
void MutateNeRnd(Bicho* bicho, Pool* pool);
void MutateLateral(Bicho* bicho, Pool* pool);
void Run(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal);
void RunDbg(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal);
