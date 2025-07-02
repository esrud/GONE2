MAXLOCI=2000000
MAXIND=2000
CC=g++
CLUSTERC=/DATA/APPS/gcc/7.2.0/bin/g++
COMMON_FLAGS=-Wall -fopenmp
MAC_PATH_OPENMP=/opt/homebrew/Cellar/libomp/19.1.5
MAC_FLAGS=-Wall -Xpreprocessor -fopenmp -lomp -I$(MAC_PATH_OPENMP)/include -L$(MAC_PATH_OPENMP)/lib -std=c++11
CFLAGS=-O3 -DMAXLOCI=$(MAXLOCI) -DMAXIND=$(MAXIND)
FASTFLAGS =-Ofast
OFNAME=gone2
LIBFILES=lib/*.cpp

all: gone

info:
	@echo
	@echo Compiling GONE with MAXLOCI = $(MAXLOCI) and MAXIND = $(MAXIND)
	@echo
	@echo To change the parameters you can run the following command
	@echo
	@echo "        make MAXLOCI=YOUR_MAX_LOCI MAXIND=YOUR_MAX_IND rule"
	@echo
	@echo where:
	@echo "    YOUR_MAX_LOCI is the maximum number of loci to process"
	@echo "    YOUR_MAX_IND is the maximum number of individuals to process"
	@echo "    RULE is one of gone, enrique, static"
	@echo
	@echo e.g\:
	@echo "        make MAXLOCI=10000000 MAXIND=5000 gone"
	@echo
	@echo Please note that increasing the numbers will need more RAM
	@echo to be able to run it
	@echo
#	@read -p "Press any key to continue..." pp
	@echo
gone: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -o $(OFNAME) gone2.cpp $(LIBFILES)
clean:
	rm $(OFNAME) 
fast: info
	$(CC) $(COMMON_FLAGS) $(FASTFLAGS) -o $(OFNAME) gone2.cpp $(LIBFILES)
static: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) gone2.cpp $(LIBFILES)
enrique: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -mcmodel=medium -o $(OFNAME) gone2.cpp $(LIBFILES)
arch: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -o $(OFNAME) -mavx2 -march=native gone2.cpp $(LIBFILES)
unoptimized: info
	$(CC) $(COMMON_FLAGS) -mcmodel=large -o $(OFNAME) gone2.cpp $(LIBFILES)
debug: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -g -o $(OFNAME) gone2.cpp $(LIBFILES)
cluster: info
	$(CLUSTERC) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) gone2.cpp $(LIBFILES)
profile: info
	clang++ -fopenmp -g -O3 -mavx2 -Wall -pedantic -o gone2_perf gone2.cpp $(LIBFILES)
dbg: info
	$(CC) -fopenmp -g -O3 -Wall -pedantic -o gone2_perf gone2.cpp $(LIBFILES)
macos: info
	$(CC) $(MAC_FLAGS) $(CFLAGS) -o $(OFNAME) gone2.cpp $(LIBFILES)
