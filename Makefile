MAXLOCI = 2000000
MAXIND  = 2000

CC      = g++
CLUSTERC = /DATA/APPS/gcc/7.2.0/bin/g++

COMMON_FLAGS = -Wall -fopenmp
# GONE_COMBO is on for every target so the GA is identical across builds:
# gone2 and gone2_tui must produce the same estimates (they differ only
# in the TUI).
CFLAGS       = -O3 -DMAXLOCI=$(MAXLOCI) -DMAXIND=$(MAXIND) -DGONE_COMBO

OFNAME   = gone2
NCURSES_OFNAME = gone2_tui
SOURCES  = gone2.cpp lib/*.cpp
NCURSES_LIBS = -lncursesw -lpthread

all: gone

info:
	@echo
	@echo "Compiling GONE2 with MAXLOCI=$(MAXLOCI) and MAXIND=$(MAXIND)"
	@echo "  approx RAM: ~$(shell awk 'BEGIN{printf "%.1f", $(MAXLOCI)*$(MAXIND)/1e9}') GB"
	@echo "  (genotype matrix is MAXLOCI x MAXIND bytes, 1 byte/genotype; it dominates the footprint)"
	@echo
	@echo "Raise either limit to process more loci/individuals; RAM grows as their product:"
	@echo "    make MAXLOCI=10000000 MAXIND=5000 gone   # ~50.0 GB - make sure it fits in memory"
	@echo
	@echo "Available targets: gone, gone-tui"
	@echo

# Combo runtime (all targets, via -DGONE_COMBO in CFLAGS): reads the input
# and computes d² exactly once, then runs the Ne-estimation GA three times —
# once per algorithm (trunc05_kick, L2, L1+kicks) — switching the smoothness
# penalty and kick schedule at runtime between passes, and keeps whichever
# pass produced the lowest bin-weighted d² residual for the _GONE2_Ne /
# _GONE2_d2 output. gone and gone-tui run this identical GA; the gone-tui
# build only adds the live TUI (current pass label, best-fit curve).
#
# Mutation rates use the original empirical defaults
# (frecMut=0.3, efectoMut=0.3). The 3×3 sweep on data/simus (May
# 2026) found these near-optimal at full strength; a static
# frecMut=0.5 and a Rechenberg 1/5 self-adapting efectoMut both
# improved fast-preset RMSE by 5–7 % but the gap closed at full
# strength, so neither replaced the defaults.
# (GA_FREC_MUT / GA_EFECTO_MUT / GA_FREC_MUT_LATERAL remain
# overridable from the command line for future sweeps.)
gone: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -o $(OFNAME) $(SOURCES)

gone-tui: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -DGONE_NCURSES_TUI \
	    -o $(NCURSES_OFNAME) $(SOURCES) $(NCURSES_LIBS)

cluster: info
	$(CLUSTERC) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) $(SOURCES)

clean:
	rm -f $(OFNAME) $(NCURSES_OFNAME)

.PHONY: all info gone gone-tui cluster clean
