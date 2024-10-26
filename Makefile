MAXLOCI=2000000
MAXIND=1000
CXX=g++
CLUSTERC=/DATA/APPS/gcc/7.2.0/bin/g++
COMMON_FLAGS=-Wall -fopenmp
CFLAGS=-O3 -DMAXLOCI=$(MAXLOCI) -DMAXIND=$(MAXIND)
FASTFLAGS =-Ofast
OFNAME=gone2
LIBFILES=lib/*.cpp

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
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -o $(OFNAME) gone2.cpp $(LIBFILES)
clean:
	rm $(OFNAME) 
fast: info
	$(CXX) $(COMMON_FLAGS) $(FASTFLAGS) -o $(OFNAME) gone2.cpp $(LIBFILES)
static: info
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) gone2.cpp $(LIBFILES)
enrique: info
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -mcmodel=medium -o $(OFNAME) gone2.cpp $(LIBFILES)
arch: info
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -o $(OFNAME) -mavx2 -march=native gone2.cpp $(LIBFILES)
unoptimized: info
	$(CXX) $(COMMON_FLAGS) -mcmodel=large -o $(OFNAME) gone2.cpp $(LIBFILES)
debug: info
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -g -o $(OFNAME) gone2.cpp $(LIBFILES)
cluster: info
	$(CLUSTERC) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) gone2.cpp $(LIBFILES)
profile: info
	clang++ -fopenmp -g -O3 -mavx2 -Wall -pedantic -o gone2_perf gone2.cpp $(LIBFILES)
static: info
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) gone2.cpp $(LIBFILES)
dbg: info
	$(CXX) -fopenmp -g -O3 -Wall -pedantic -o gone2_perf gone2.cpp $(LIBFILES)
