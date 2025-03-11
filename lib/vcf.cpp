#include "vcf.hpp"

bool ReadVcf(std::string fichVcf, PopulationInfo *popInfo) {
  /*
   * Takes as input the file name and a pointer to the population
   * matrix and returns if the file was processed correcty
   */
  char base1, base2;
  int contaIndBase = 0, nline=0, i, j;
  int conta = 0, posi = 0, posi2 = 0, longi = 0;
  std::string line;
  std::string cromocod;
  std::string cromocodback = "laksjhbqne";
  std::vector<int> rangocromo(MAXCROMO);
  std::vector<char> base(MAXLOCI);
  std::vector<double> chrnum(MAXCROMO);
  std::vector<double> chrprop(MAXCROMO);
  std::string str_vcf = "AGCT,agct";
  std::string str_base = "";
  bool rightletter,hayalelos;

  // READING .vcf DATA:
  std::ifstream entrada;
  // Bucle de lectura del fichero tped
  entrada.open(fichVcf, std::ios::in);
  if (!entrada.good()) {
    std::cerr << "Could not open \"" << fichVcf << "\". Does the file exist?"
              << std::endl;
    return false;
  }
  int ncromos = 0;
  int contalines = 0;
  int nLoci = 0;
  int nIndi = 0;
  while (std::getline(entrada, line)) {
    ++nline;
    if (line.at(0) == '#') {
      continue;
    }
    longi = static_cast<int>(line.length());
    if (longi < 12) {
      std::cerr << "Line too short in vcf file" << std::endl;
      return false;
    }
    conta = 0;
    posi = 0;
    // chr name
    posi2 = posi;
    posi = static_cast<int>(line.find_first_of("\t", posi2));
    if (posi < 0) {
      std::cerr << "Line too short in vcf file" << std::endl;
      return false;
    }
    cromocod = line.substr(0, posi);
    ++posi;
    ++conta;

    posi2 = posi;
    posi = static_cast<int>(line.find_first_of("\t", posi2));
    if (posi < 0) {
      std::cerr << "Line too short in vcf file" << std::endl;
      return false;
    }
    int position_basepair = std::stoi(line.substr(posi2, posi-posi2));
    popInfo->posiBP[nLoci] = position_basepair;
    ++posi;
    ++conta;

    while ((posi < longi) && (conta < 3)){
        posi2=posi;
        posi=int(line.find_first_of("\t",posi2));
        if (posi < 0) {
            std::cerr << "Line too short in vcf file" << std::endl;
            exit(EXIT_FAILURE);
        }
        ++posi;
        ++conta;
    }

    // Mira el alelo REF
    posi2=posi;
    posi=int(line.find_first_of("\t",posi2));
    j=posi-posi2;
    rightletter=true;
    hayalelos=true;
    for (i=0;i<j;++i){
        str_base=line.substr(posi2+i,1);
        if (str_base=="."){
            hayalelos=false;
        }
        else if (str_vcf.find(str_base) == std::string::npos){
            rightletter=false;
            break;
        }
    }
    if (!rightletter) {
        std::cerr << "Wrong reference allele in line "<<nline<<" of vcf file. Only bases A, G, C and T (case-insensitive) are allowed." << std::endl;
        exit(EXIT_FAILURE);
    }
    ++posi;
    ++conta;

    // Mira el alelo ALT
    posi2=posi;
    posi=int(line.find_first_of("\t",posi2));
    j=posi-posi2;
    for (i=0;i<j;++i){
        str_base=line.substr(posi2+i,1);
        if (str_base=="."){
            hayalelos=false;
        }
        else if (str_vcf.find(str_base) == std::string::npos){
            rightletter=false;
            break;
        }
    }
    if (!rightletter) {
        std::cerr << "Wrong alternative allele in line "<<nline<<" of vcf file (maybe also in other lines). Only bases A, G, C and T (case-insensitive) are allowed." << std::endl;
        exit(EXIT_FAILURE);
    }
    ++posi;
    ++conta;

    while ((posi < longi) && (conta < 8)) {
      posi2 = posi;
      posi = static_cast<int>(line.find_first_of("\t", posi2));
      if (posi < 0) {
        std::cerr << "Line too short in vcf file" << std::endl;
        return false;
      }
      ++posi;
      ++conta;
    }
    // Mira si la posicion 8 es GT (genotipo)
    bool haygenotipos = line.substr(posi, 2) == "GT";

    // Lee los  genotipos:
    if (haygenotipos && hayalelos) {
      // contabiliza el cromosoma que habia leido al principio de la linea
      if (cromocod != cromocodback) {
        cromocodback = cromocod;
        rangocromo[ncromos] = contalines;
        ++ncromos;
      }
      popInfo->cromo[contalines] = ncromos;
      ++contalines;  // para crom
      // Avanza hasta el primer genotipo:
      posi2 = posi;
      posi = static_cast<int>(line.find_first_of("\t", posi2));
      if (posi < 0) {
        std::cerr << "Line too short in vcf file" << std::endl;
        return false;
      }
      ++posi;
      ++conta;
      // empieza a leer genotipos
      nIndi = 0;
      // asigna genot.
      while (posi < longi) {
        base1 = line.at(posi);
        // Haplotype check
        if (popInfo->haplotype != 1) {
          posi2 = posi;
          posi = line.find_first_of("/|", posi2);
          if (posi < 0) {
            break;
          }
          ++posi;
          base2 = line.at(posi);
        } else {
          base2 = base1;
        }
        if ((base1 != '.') && (base2 != '.')) {
          if (base[nLoci] == '\0') {
            base[nLoci] = base1;
          }
          base1 = base1 != base[nLoci] ? 'X' : base[nLoci];
          base2 = base2 != base[nLoci] ? 'X' : base[nLoci];

          if ((popInfo->haplotype == 0) || (popInfo->haplotype == 3)) {
            // Diploides Fase desconocida y Pseudohaploides
            // 0:homo ref, 1:het, 2:homo noref
            if (base1 == base2) {
              popInfo->indi[nIndi][nLoci] =
                base1 == base[nLoci] ? 0 : 2;
            } else {
              popInfo->indi[nIndi][nLoci] = 1;
            }
          } 
          else if ((popInfo->haplotype == 1)) {
            // Haploides 
            popInfo->indi[nIndi][nLoci] = (base1 == base[nLoci] ? 0 : 2);
          } 
          else {
            // Fase conocida (2)
            popInfo->indi[nIndi][nLoci] =
              base1 == base[nLoci] ? 0 : 2;
            popInfo->indi[nIndi + 1][nLoci] =
              base1 == base[nLoci] ? 0 : 2;
          }
        } 
        else {
          // '9' = Genotipo sin asignar
          // Fase desconocida, haploides y pseudohaploides
          popInfo->indi[nIndi][nLoci] = 9;
          if (popInfo->haplotype == 2) {
            // Fase conocida
            popInfo->indi[nIndi + 1][nLoci] = 9;
          }
        }
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of("\t", posi2));
        if (posi < 0) {
          posi = longi;
        }
        ++posi;

        nIndi++;
        if (popInfo->haplotype == 2){// Fase Conocida son el doble de individuos
          nIndi++;
        }
        if (nIndi > MAXIND) {
          std::cerr << "Reached limit of sample size (" << MAXIND << ")"
                    << std::endl;
          return false;
        }
      }

      if (nLoci == 0) {
        contaIndBase = nIndi;
      }

      if (nIndi != contaIndBase) {
        std::cerr << "Some SNP in the sample are not represented in all the "
                     "individuals"
                  << std::endl;
        return false;
      }

      ++nLoci;
      if (nLoci > MAXLOCI) {
        std::cerr << "Reached max number of loci (" << MAXLOCI << ")"
                  << std::endl;
        return false;
      }
    }
  }
  entrada.close();
  popInfo->numLoci = nLoci;
  popInfo->numIndi = nIndi;
  popInfo->numCromo = ncromos;
  rangocromo[popInfo->numCromo] = nLoci;
  popInfo->Mbtot = 0;
  for (int i = 0; i < popInfo->numCromo; ++i) {
    int idx = rangocromo[i];
    int next_idx = rangocromo[i+1] - 1;
    popInfo->Mbtot += popInfo->posiBP[next_idx] - popInfo->posiBP[idx];
  }
  popInfo->Mbtot /= 1000000.0;
  return true;
}
