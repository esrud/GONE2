
#include "console.hpp"

// Print the usage/help text (options and examples) to stderr.
// Params: appName — program name (argv[0]) used in the usage examples.
void printHelp(char* appName) {
  fprintf(
    stderr,
    "GONE2 - Genetic Optimization for Ne Estimation (v3.0 - Jul 2026)\n"
    "Authors: Enrique Santiago - Carlos Köpke\n"
    "       This software estimates past demography from the distribution of LD\n"
    "       between pairs of SNP located at different distances on a genetic map.\n"
    "\n"
    "USAGE: %s [OPTIONS] <file_name_with_extension>\n"
    "       Where file_name is the name of the data file in vcf, ped or tped\n"
    "           format. The filename must include the .vcf, .ped or tped\n"
    "           extension, depending on its format.\n"
    "       Estimates are made using the genetic map available in the .tped\n"
    "           file, when using the tped format, or in an accompanying .map file,\n"
    "           when using the ped or vcf formats. This .map file has the same \n"
    "           name as the .ped or .vcf data file and must be located in the\n"
    "           same directory as the data file. When a constant recombination\n"
    "           rate per Mb is assumed (option -r), the physical map is used to\n"
    "           infer an approximate genetic map, and the detailed genetic map\n"
    "           information in the data file is ignored, if available.\n"
    "       Three files are created with the following extensions:\n"
    "           _GONE_Ne: Estimates of Ne backward in time.\n"
    "           _GONE_d2: Number of SNP pairs used in the analysis, observed LD \n"
    "                 (weighted squared correlation d²) and predicted LD per bin \n"
    "                 of recombination frequency.\n"
    "           _GONE_STATS: Summary statistics.\n\n"
    "OPTIONS:\n"
    "    -h     Print out this help\n"
    "    -g     Type of genotyping data. 0:unphased diploids; 1:haploids; \n"
    "           2:phased diploids; 3:low-coverage (0 by default). Low-coverage\n"
    "           assumes diploid unphased genotypes and can be used with any\n"
    "           distribution of coverage within and between individuals.\n"
    "    -x     The sample is considered to be a random set of individuals from\n"
    "           a metapopulation with subpopulations of equal size. After the\n"
    "           After structure search converges, the GA is run against the\n"
    "           Fst-corrected within-subpop d² to produce a continuous\n"
    "           per-generation N_T curve.\n"
    "    -b     Average base calling error rate per site (0 by default). \n"
    "    -i     Number of individuals to use in the analysis (all by default)\n"
    "    -s     Number of SNPs to use in the analysis (all by default)\n"
    "    -t     Number of threads to be used in parallel computation (default: %d)\n"
    "    -l     Lower bound of recombination rates to be considered (default: 0.001)\n"
    "    -u     Upper bound of recombination rates to be considered (default: 0.05)\n"
    "    -e     Reinforcement estimates of recent generations\n"
    "    -r     If specified, constant rec rate in cM/Mb across the genome\n"
    "    -M     If specified, minor allele frequency cut-off\n"
    "    -o     Specifies the output filename. If not specified, the output \n"
    "           filename is built from the name of the input file.\n"
    "    -S     Integer to seed random number generator. Taken \n"
    "           from the system, if not given.\n"
    "    -E     Overlay on the TUI chart the min/max Ne bounding curves\n"
    "           across the GA rounds plus a density cloud of the per-round\n"
    "           Ne values, alongside the best-fit curve.\n"
    "    -f     Path to a reference Ne-history file.\n"
    "\n"
    "EXAMPLES:\n"
    "    - Analysis of high quality diploid unphased data in \"file.ped\" (PLINK\n"
    "      format) assumes a constant recombination rate of 1.1 cM per Mb across\n"
    "      the genome (no need for a detailed genetic map within the .map file).\n"
    "      16 threads will be used:\n"
    "          %s -r 1.1 -t 16 file.ped\n"
    "    - A subsample of 10000 SNPs of the individuals in \"file.ped\" assuming\n"
    "      assuming that they were randomly sampled from a metapopulation composed\n"
    "      of two subpopulations:\n"
    "          %s -x -s 100000 file.ped\n"
    "    - Analysis of diploid high quality phased data in \"file.vcf\" (format vcf).\n"
    "      assumes that the genetic locations of the SNPs are given in the\n"
    "      \"file.map\" file (PLINK format) available in the same directory:\n"
    "          %s -g 2 file.vcf\n"
    "    - Analysis of diploid high quality phased data in \"file.vcf\" (format vcf).\n"
    "      assumes a constant recombination rate of 1.1 cM per Mb across the genome:\n"
    "          %s -g 2 -r 1.1 file.vcf\n"
    "    - Analysis of diploid high quality unphased data in a .tped file,\n"
    "      performed on a radom subset of 50 individuals and 100,000 SNPs.\n"
    "          %s -i 50 -s 10000 file.tped\n"
    "    - Analysis of low quality unphased data in a .tped file containing the\n"
    "      locations on a genetic map. Low-coverage (no need to specify depth) and\n"
    "      a genotyping error rate of 0.001 across genomes are assumed.\n"
    "          %s -g 3 -b 0.001 file.tped\n"
    "\n",
    appName, omp_get_max_threads(), appName, appName, appName, appName, appName, appName);
}

// Parse the command-line options into params, starting from the defaults.
// Applies each flag, validates the input file and its type, and exits the
// program on -h or on an invalid/missing argument.
// Params: argc/argv — command-line arguments; params — AppParams filled in
//   place with the parsed run configuration.
void HandleInput(int argc, char * argv[], AppParams* params) {
  SetDefaultParameters(params);
  for (;;) {
    switch (getopt(argc, argv, "hxg:b:n:i:s:l:u:er:M:t:qpvS:b:zvo:f:EB:")) {
      case '?':
      case 'h':
      default:
        printHelp(argv[0]);
        exit(EXIT_FAILURE);
      case 'g':
        params->haplotype = std::atoi(optarg);
        if (params->haplotype > 3 || params->haplotype < 0) {
          std::cerr << "Invalid code for type of genotyping data" << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 'x':
        params->mix = true;
        continue;
      case 'b':
        params->basecallerror = std::atof(optarg);
        if (params->basecallerror > 0.2 || params->basecallerror < 0) {
          std::cerr << "Invalid base call error rate. Must be < 0.2"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 'm':
        params->miss = std::atof(optarg);
        if (params->miss < 0 || params->miss > 0.90) params->miss = 0.2;
        continue;
      case 'i':
        params->numSample = std::atoi(optarg);
        if (params->numSample < 3) {
          std::cerr << "Invalid sample number" << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 's':
        params->numSNPs = std::atoi(optarg);
        if (params->numSNPs < 10) {
          std::cerr << "Invalid number of SNPs" << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 'l':
        params->flags |= FLAG_LC;
        params->lc = std::atof(optarg);
        if (params->lc <= 0 || params->lc > 0.01) {
          std::cerr << "Invalid minimum recombination rate c. Must be: "
                       "0< c <=0.01." << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 'u':
        params->flags |= FLAG_HC;
        params->hc = std::atof(optarg);
        if (params->hc <= 0.01 || params->hc > 0.15) {
          std::cerr << "Invalid maximum recombination rate c. Must be: "
                       "0.01< c <=0.15" << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 'e':
        params->hayrecentbins = true;
        continue;
      case 'r':
        params->cMMb = std::atof(optarg);
        if (params->cMMb <= 0) {
          std::cerr << "Invalid ratio cM/Mb" << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 'M':
        params->MAF = std::atof(optarg);
        if (params->MAF < 0 || params->MAF > 0.4) {
          std::cerr << "Invalid MAF value" << std::endl;
          exit(EXIT_FAILURE);
        }
        continue;
      case 't':
        params->numThreads = std::atoi(optarg);
        continue;
      case 'B':
        params->sizeBins = std::atoi(optarg);
        if (params->sizeBins < 1) {
          std::cerr << "Invalid minimum number of SNP pairs per bin. Must be"
                       " at least 1." << std::endl;
          exit(EXIT_FAILURE);
        }
        params->flags |= FLAG_RESIZE_BINS;
        continue;
      case 'q':
        params->quiet = true;
        continue;
      case 'p':
        params->printToStdOut = true;
        continue;
      case 'o':
        params->fileOut = optarg;
        continue;
      case 'f':
        params->realNeFile = optarg;
        continue;
      case 'S':
        params->semilla = std::atoi(optarg);
        continue;
      case 'E':
        params->showNeBand = true;
        continue;
      case 'z':
        params->flags |= FLAG_REP;
        continue;
      case 'v':
        params->flags |= FLAG_DEBUG;
        continue;
      case -1:
        break;
    }
    break;
  }

  if (optind < argc) {
    params->fich = argv[optind];
  }

  if (params->fich == "") {
    std::cerr << "Missing data file name" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (params->fileOut == "") {
    params->fileOut = params->fich;
  }
  bool name_ok = GetFileType(params->fich, &params->ftype);
  if ((!name_ok) || ((params->ftype!="tped") && (params->ftype!="ped") &&(params->ftype!="vcf"))) {
    std::cerr << "Invalid file name extension. Must be .vcf, .ped or .tped. " << std::endl;
    exit(EXIT_FAILURE);
  }
  if (params->haplotype == 2) {
    params->numSample*=2;
  }

  if (params->haplotype == 3){
      params->hayrecentbins=false;
      if (params->hc > 0.05){
        params->hc = 0.05;
      }
  }
  if (params->mix){
    if ((params->haplotype !=0)){
      std::cerr << "There is no possibility of analyzing phased-diploid, haploid and low-coverage data to infer Ne of metapopulations." << std::endl;
      exit(EXIT_FAILURE);
    }
    if (params->ngensampling != 1){
      std::cerr << "There is no possibility of analyzing more than one generation of sampling." << std::endl;
      exit(EXIT_FAILURE);
    }
    params->hayrecentbins = false;
    params->flags |= FLAG_LC;
    params->flags |= FLAG_HC;
    params->lc = 0.001;
    params->hc = 0.15;
    params->twoPassGA = true;
    params->mixRun = true;
    if (params->numSNPs == 0) {
      params->numSNPs = 20000;
    }
  }
  if (params->hayrecentbins){
    params->flags |= FLAG_HC;
    params->hc = 0.09;
  }
}

// Initialise an AppParams to its default run configuration (before any
// command-line flags are applied).
// Params: params — AppParams filled in place with defaults.
void SetDefaultParameters(AppParams* params) {
  params->haplotype = 0;
  params->mix = false;
  params->basecallerror = 0;
  params->miss = 0.2;
  params->coverage = 0;
  params->ngensampling = 1; 
  params->numThreads = 0;
  params->numSample = 0;
  params->numSNPs = 0;
  params->hc = 0.05;
  params->lc = 0.001;
  params->cMMb = 0;
  params->MAF = 0;
  params->distance = 1;
  params->quiet = false;
  params->printToStdOut = false;
  params->flags = 0;
  params->semilla = 0;
  params->muestraSalida = kMuestraSalida;
  params->nbins = 50;
  params->sizeBins = 1;
  params ->hayrecentbins = false;
  params->twoPassGA = false;
  params->mixRun = false;
  params->showNeBand = false;
  params->mixRoundProbe = false;
  params->mixConverged = true;
  params->mixStructureEstimated = false;
  params->mixFatal = false;
  params->mixFailReason.clear();
  params->analysisFailReason.clear();
}

// Determine the input file type from its extension (ped/vcf/tped).
// Params: fname — input file name; ftype — set to the lowercased extension
//   on success.
// Returns: true if the extension is a recognised type, false otherwise.
bool GetFileType(std::string fname, std::string *ftype) {
  std::string delimiter = ".";
  size_t del_start = fname.rfind(delimiter);
  if (del_start == std::string::npos) {
    return false;
  }
  del_start++;
  *ftype = fname.substr(del_start, fname.size() - del_start);
  return true;
}
