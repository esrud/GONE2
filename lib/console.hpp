#pragma once

#include <getopt.h>
#include <omp.h>
#include <iostream>
#include <string>
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
    int flags;
    int muestraSalida;
    int semilla;
    int sizeBins;
    int nbins;
    bool hayrecentbins;
    ProgressStatus progress;
} AppParams;


void HandleInput(int argc, char * argv[], AppParams* params);
void SetDefaultParameters(AppParams* params);
bool GetFileType(std::string fname, std::string *ftype);
