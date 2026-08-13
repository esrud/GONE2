#pragma once

#include <climits>

#ifndef MAXLOCI
#define MAXLOCI 2000000
#endif
#ifndef MAXIND
#define MAXIND 2000
#endif
static_assert(MAXLOCI >= 1 && MAXLOCI <= INT_MAX,
              "MAXLOCI must be in [1, INT_MAX]");
static_assert(MAXIND >= 1 && MAXIND <= INT_MAX,
              "MAXIND must be in [1, INT_MAX]");
#define MAXCROMO 1000
#define MAXBINS 2000
#ifndef GONE_ROUNDS
#define GONE_ROUNDS 50
#endif
#ifndef GA_BIN_WEIGHTED
#define GA_BIN_WEIGHTED 1
#endif
#define MAXDIST 100000

#define MAX_POPULATION 128
#define kNumLinMax 256
#define kMaxD2PredBins 64
static_assert(kNumLinMax >= 132, "kNumLinMax must cover transient GA splices");
static_assert(kMaxD2PredBins >= 60, "kMaxD2PredBins must cover d2 bins");
#define kNumGenMax 2000
#define kResolucion 3
#define kMuestraSalida 10
#define kNumBins 50
#define kMinBinsForGA 15

#define kTopePosInicio1 60
#define kTopePosInicio2 120
#define kTopePosInicio3 240
#define kTopePosInicio4 480

#define kNumDesInicio 3000
#define kNumDes 1000
#ifndef kNumGenerations
#define kNumGenerations 750
#endif
#define kTopeSalto 9.0
#define kInvTopeSalto (1.0 / kTopeSalto)
#define kTopeSalto2 50.0
#define kInvTopeSalto2 (1.0 / kTopeSalto2)
#define kEfectoMutSuave 0.02
#define kFrecInversion 0.3

#define kTercio1 10
#define kTercio2 90
#define kTercio12 (kTercio1 + kTercio2)
#define kNHijos (kTercio1 + kTercio2)

#define FLAG_RESIZE_BINS (1 << 6)
#define FLAG_DEBUG (1 << 5)
#define FLAG_REP (1 << 4)
#define FLAG_NE (1 << 3)
#define FLAG_LC (1 << 2)
#define FLAG_HC (1 << 1)
#define FLAG_NBINS 1
#define FLAG_POR_BLOQUES true

#define DEFAULT_BIN_SIZE 10000
#define DEFAULT_BIN_SIZE_MIX 1000
#define MAX_NE_SIZE 10000000.0
#define MIN_NE_SIZE 1.0

#define MAX_DOUBLE std::numeric_limits<double>::max()
