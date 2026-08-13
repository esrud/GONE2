#pragma once

#include <getopt.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "constants.hpp"
#include "progress.hpp"

typedef struct {
    int haplotype;
    bool mix;
    double basecallerror;
    double miss;
    double coverage;
    double ngensampling;
    int numThreads;
    int numSample;
    long int numSNPs;
    double hc;
    double lc;
    double cMMb;
    double MAF;
    int distance;
    bool quiet;
    bool printToStdOut;
    std::string fich;
    std::string fileOut;
    std::string ftype;
    std::string realNeFile;
    int flags;
    int muestraSalida;
    int semilla;
    int sizeBins;
    int nbins;
    bool hayrecentbins;
    bool twoPassGA;
    bool mixRun;
    bool showNeBand;
    bool mixRoundProbe;
    bool mixConverged;
    bool mixStructureEstimated;
    bool mixFatal;
    std::string mixFailReason;
    std::string analysisFailReason;
    std::vector<double> chartNTMix;
    ProgressStatus progress;
} AppParams;

// Infix prepended before each _GONE2_* output-file suffix: "_GONE2_mix_"
// for an -x (metapopulation) run, otherwise "_GONE2_".
// Params: mixRun — whether this is an -x run.
// Returns: the filename prefix.
inline std::string Gone2OutPrefix(bool mixRun) {
  return mixRun ? "_GONE2_mix_" : "_GONE2_";
}

void HandleInput(int argc, char * argv[], AppParams* params);
void SetDefaultParameters(AppParams* params);
bool GetFileType(std::string fname, std::string *ftype);
