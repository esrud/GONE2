
#include "tped.hpp"
#include <ostream>

bool ReadTped(std::string fichTped, PopulationInfo *popInfo) {
  /*
   * Takes as input the file name and a pointer to the population
   * matrix and returns the number ofHR	Chromosome code
   * These are the columns of a tped file
      CHR       Chromosome code
      SNP	Variant identifier
      (C)M	Position in morgans or centimorgans
      POS	Base-pair coordinate
      COUNTED	Counted allele (defaults to A1)
      ALT	Other allele(s), comma-separated
      <FID>_<IID>...	Allelic dosages (0/1/2/'NA' for diploid variants, 0/2/'NA' for haploid) individuals in the file
   */
  char base1, base2;
  int contaIndBase = 0, nline =0;
  int conta = 0, posi = 0, posi2 = 0, longi = 0;
  std::string line;
  std::string cromocod;
  std::string cromocodback = "laksjhbqne";
  std::vector<int> rangocromo(MAXCROMO);
  std::vector<char> base(MAXLOCI);
  std::vector<double> chrnum(MAXCROMO);
  std::vector<double> chrprop(MAXCROMO);
  std::string str_plink = "AGCTNagctn0123456789";
  std::string str_base = "";
  bool rightletter;

  // READING .tped DATA:
  std::ifstream entrada;
  // Bucle de lectura del fichero tped
  entrada.open(fichTped, std::ios::in);
  if (!entrada.good()) {
    std::cerr << "Could not open \"" << fichTped << "\". Does the file exist?"
              << std::endl;
    return false;
  }
  int ncromos = 0;
  int contalines = 0;
  int nLoci = 0;
  int nIndi = 0;
  while (std::getline(entrada, line)) {
    longi = static_cast<int>(line.length());
    if (longi < 12) {
      std::cerr << "Line too short in tped file" << std::endl;
      return false;
    }
    conta = 0;
    posi = 0;
    // chr name
    posi2 = posi;
    posi = static_cast<int>(line.find_first_of(" \t", posi2));
    if (posi < 0) {
      std::cerr << "Line too short in tped file" << std::endl;
      return false;
    }
    cromocod = line.substr(0, posi);
    if (cromocod != cromocodback) {
      cromocodback = cromocod;
      rangocromo[ncromos] = contalines;
      ++ncromos;
    }
    popInfo->cromo[contalines] = ncromos;

    ++contalines;  // para crom

    ++posi;
    ++conta;

    // Ignora las primeras columnas
    while ((posi < longi) && (conta < 2)) {
      posi2 = posi;
      posi = static_cast<int>(line.find_first_of(" \t", posi2));
      if (posi < 0) {
        std::cerr << "Line too short in tped file" << std::endl;
        return false;
      }
      ++posi;
      ++conta;
    }
    // Position in morgans or centimorgans
    posi2 = posi;
    posi = static_cast<int>(line.find_first_of(" \t", posi2));
    double position_centimorgans = std::stod(line.substr(posi2, posi-posi2));
    popInfo->posiCM[nLoci] = position_centimorgans;
    ++posi;
    ++conta;
    // Base-pair coordinate
    posi2 = posi;
    posi = static_cast<int>(line.find_first_of(" \t", posi2));
    int position_basepair = std::stoi(line.substr(posi2, posi-posi2));
    popInfo->posiBP[nLoci] = position_basepair;

    ++posi;
    ++conta;
    if (conta == 4) {
      nIndi = 0;
      // asigna genot.
      while (posi < longi) {
        base1 = line.at(posi);
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of(" \t", posi2));
        if (posi < 0) {
          break;
        }
        ++posi;
        base2 = line.at(posi);
        rightletter=true;
        str_base=base1;
        if (str_plink.find(str_base) == std::string::npos){rightletter=false;}
        str_base=base2;
        if (str_plink.find(str_base) == std::string::npos){rightletter=false;}
        if (!rightletter) {
            std::cerr << "Wrong allele in line "<<nline<<" of tped file (maybe also in other lines). Only bases A, G, C, T and N (case-insensitive) and numbers from 0 to 9 are allowed." << std::endl;
            exit(EXIT_FAILURE);
        }

        // TODO(me): What about the haplotype?
        if ((base1!='0') && (base2!='0') && (base1!='N') && (base2!='N') && (base1!='n') && (base2!='n')) {
          if (base[nLoci] == '\0') {
            base[nLoci] = base1;
          }
          if (base1 != base[nLoci]) {
            base1 = 'X';
          }
          if (base2 != base[nLoci]) {
            base2 = 'X';
          }

          if ((popInfo->haplotype == 0) || (popInfo->haplotype == 3)){
            // Fase desconocida y pseudohaploides
            // 0:homo ref, 1:het, 2:homo noref
            if (base1 == base2) {
              popInfo->indi[nIndi][nLoci] = (base1 == base[nLoci] ? 0 : 2);
            } else {
              popInfo->indi[nIndi][nLoci] = 1;
            }
          } 
          else if (popInfo->haplotype ==1){ // HAPLOIDES
            if (base1 == base[nLoci]) {
              popInfo->indi[nIndi][nLoci] = 0;
            } else {
              popInfo->indi[nIndi][nLoci] = 2;
            }
          }
          else {
            // Fase conocida
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
            // FASE CONOCIDA (2) genera dos individuos haploides
            popInfo->indi[nIndi + 1][nLoci] = 9;
          }
        }
        posi2 = posi;
        posi = line.find_first_of(" \t", posi2);
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
  popInfo->Mtot = 0;
  popInfo->Mbtot = 0;
  for (int i = 0; i < popInfo->numCromo; ++i) {
    int idx = rangocromo[i];
    int next_idx = rangocromo[i+1] - 1;
    popInfo->Mtot  += popInfo->posiCM[next_idx] - popInfo->posiCM[idx];
    popInfo->Mbtot += popInfo->posiBP[next_idx] - popInfo->posiBP[idx];
  }
  popInfo->Mtot /= 100.0;       // en Morgans
  popInfo->Mbtot /= 1000000.0;
  return true;
}
