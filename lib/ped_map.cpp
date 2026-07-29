
#include "ped_map.hpp"
#include "./pool.hpp"

namespace {
struct PedAlleleTable {
  bool ok[256];
  // std::string's constructor is non-constexpr until C++20, so walk
  // a string literal by pointer instead — that's constexpr from C++14.
  constexpr PedAlleleTable() : ok{} {
    for (const char* p = "AGCTNagctn0123456789"; *p; ++p) {
      ok[static_cast<unsigned char>(*p)] = true;
    }
  }
};
const PedAlleleTable kPedAllele;
// Test whether a character is an accepted .ped allele symbol.
// Params: c — allele character.
// Returns: true if c is a valid allele code, else false.
inline bool IsValidPedAllele(char c) {
  return kPedAllele.ok[static_cast<unsigned char>(c)];
}
}

// Read a PLINK .ped file (six pedigree columns then two chars per locus)
// into popInfo->indi[][], setting the locus and individual counts.
// Params: fichPed — path to the .ped file; popInfo — population struct
//   filled in place.
// Returns: true on success, false if the file can't be read or is malformed.
bool ReadPed(std::string fichPed, PopulationInfo* popInfo) {
  uint8_t base[MAXLOCI] = {0};
  uint8_t base1, base2;
  int contaLociBase = 0, nline = 0;
  int conta = 0, posi = 0, posi2 = 0, longi = 0;
  std::string line;

  std::ifstream entrada;
  entrada.open(fichPed, std::ios::in);
  if (!entrada.good()) {
    std::cerr << "Could not open \"" << fichPed << "\". Does the file exist?"
              << std::endl;
    return false;
  }
  while (std::getline(entrada, line)) {
    longi = static_cast<int>(line.length());
    if (longi < 12) {
      std::cerr << "Line too short in ped file" << std::endl;
      return false;
    }
    conta = 0;
    posi = 0;
    while ((posi < longi) && (conta < 6)) {
      posi2 = posi;
      posi=static_cast<int>(line.find_first_of(" \t",posi2));
      if (posi < 0) {
        break;
      }
      ++posi;
      ++conta;
    }
    if (conta == 6) {
      popInfo->numLoci = 0;
      if ((popInfo->haplotype == 2 ? popInfo->numIndi + 1
                                   : popInfo->numIndi) >= MAXIND) {
        std::cerr << "Reached limit of sample size (" << MAXIND << ")"
                  << std::endl;
        return false;
      }
      while (posi < longi) {
        base1 = line.at(posi);
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of(" \t", posi2));
        if (posi < 0) break;
        ++posi;
        base2 = line.at(posi);

        if (!IsValidPedAllele(base1) || !IsValidPedAllele(base2)) {
          std::cerr << "Wrong allele in line " << nline
                    << " of ped file (maybe also in other lines). Only bases "
                       "A, G, C, T and N (case-insensitive) and numbers from "
                       "0 to 9 are allowed."
                    << std::endl;
          exit(EXIT_FAILURE);
        }

        const bool missing = (base1 == '0' || base2 == '0' ||
                              base1 == 'N' || base2 == 'N' ||
                              base1 == 'n' || base2 == 'n');
        const int locusIdx = popInfo->numLoci;
        if (locusIdx >= MAXLOCI) {
          std::cerr << "Reached max number of loci (" << MAXLOCI << ")"
                    << std::endl;
          return false;
        }
        if (!missing) {
          if (base[locusIdx] == '\0') base[locusIdx] = base1;
          const uint8_t ref = base[locusIdx];
          if (base1 != ref) base1 = 'X';
          if (base2 != ref) base2 = 'X';

          if (popInfo->haplotype == 0 || popInfo->haplotype == 3) {
            if (base1 == base2) {
              popInfo->indi[popInfo->numIndi][locusIdx] =
                  (base1 == ref ? 0 : 2);
            } else {
              popInfo->indi[popInfo->numIndi][locusIdx] = 1;
            }
          } else if (popInfo->haplotype == 1) {
            popInfo->indi[popInfo->numIndi][locusIdx] =
                (base1 == ref ? 0 : 2);
          } else {
            popInfo->indi[popInfo->numIndi][locusIdx] =
                (base1 == ref ? 0 : 2);
            popInfo->indi[popInfo->numIndi + 1][locusIdx] =
                (base2 == ref ? 0 : 2);
          }
        } else {
          popInfo->indi[popInfo->numIndi][locusIdx] = 9;
          if (popInfo->haplotype == 2) {
            popInfo->indi[popInfo->numIndi + 1][locusIdx] = 9;
          }
        }
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of(" \t", posi2));
        if (posi < 0) posi = longi;
        ++posi;
        ++popInfo->numLoci;
        if (popInfo->numLoci > MAXLOCI) {
          std::cerr << "Reached max number of loci (" << MAXLOCI << ")"
                    << std::endl;
          std::cerr << "\nTo increase the number of loci, recompile the program "
                    << "increasing the maximum number of loci. Please note that "
                    << "doing so will also require more RAM to run. To see how "
                    << "to do it, you can run:"
                    << "\n    make info\n"
                    << std::endl;
          return false;
        }
      }

      if (popInfo->numIndi == 0) {
        contaLociBase = popInfo->numLoci;
      }

      if (popInfo->numLoci != contaLociBase) {
        std::cerr << "Some genomes in the sample are of different sizes"
                  << std::endl;
        return false;
      }

      popInfo->numIndi++;
      if (popInfo->haplotype == 2){
        popInfo->numIndi++;
      }
      if (popInfo->numIndi > MAXIND) {
        std::cerr << "Reached limit of sample size (" << MAXIND << ")"
                  << std::endl;
        std::cerr << "\nTo increase the sample size, recompile the program "
                  << "increasing the maximum number of individuals. Please note that "
                  << "doing so will also require more RAM to run. To see how "
                  << "to do it, you can run:"
                  << "\n    make info\n"
                  << std::endl;
        return false;
      }
    }
  }
  entrada.close();
  return true;
}

// Read a PLINK .map file (chromosome, SNP id, genetic position, bp) into
// popInfo's per-locus position and chromosome arrays and the total map
// length.
// Params: fichMap — path to the .map file; popInfo — population struct
//   filled in place.
// Returns: true on success, false if the file can't be read or is malformed.
bool ReadMap(std::string fichMap, PopulationInfo* popInfo) {
  int posi = 0, posi2 = 0, longi = 0;
  std::string line;

  std::string cromocod;
  std::string cromocodback = "laksjhbqne";

  std::ifstream entrada;
  entrada.open(fichMap, std::ios::in);
  if (!entrada.good()) {
    std::cerr << "Could not open \"" << fichMap << "\". Does the file exist?"
              << std::endl;
    return false;
  }
  int contalines = 0;
  int rangocromo[MAXCROMO] = {0};
  while (std::getline(entrada, line)) {
    longi = static_cast<int>(line.length());
    if (longi < 5) {
      std::cerr << "Line too short in map file" << std::endl;
      return false;
    }

    posi=static_cast<int>(line.find_first_of(" \t",0));
    if (posi <= 0) {
      std::cerr << "Empty line in map file" << std::endl;
      return false;
    }
    cromocod = line.substr(0, posi);
    if (cromocod != cromocodback) {
      cromocodback = cromocod;
      rangocromo[popInfo->numCromo] = contalines;
      ++popInfo->numCromo;
    }
    popInfo->cromo[contalines] = popInfo->numCromo;

    ++posi;
    posi2 = posi;
    posi=static_cast<int>(line.find_first_of(" \t",posi2));
    if (posi <= 0) {
      std::cerr << "Error in map file (1)" << std::endl;
      return false;
    }

    ++posi;
    posi2 = posi;
    posi=static_cast<int>(line.find_first_of(" \t",posi2));
    if (posi <= 0) {
      std::cerr << "Error in map file (2)" << std::endl;
      return false;
    }
    popInfo->posiCM[contalines] = std::stod(line.substr(posi2, posi - posi2));
    ++posi;

    if (posi > longi) {
      std::cerr << "Error in map file (3)" << std::endl;
      return false;
    }
    popInfo->posiBP[contalines] = std::stoi(line.substr(posi, longi - posi));
    ++contalines;
  }
  rangocromo[popInfo->numCromo] = contalines;
  entrada.close();
  popInfo->Mtot = 0;
  popInfo->Mbtot = 0;
  for (int conta = 0; conta < popInfo->numCromo; ++conta) {
    popInfo->Mtot += popInfo->posiCM[rangocromo[conta+1]-1]-popInfo->posiCM[rangocromo[conta]];
    popInfo->Mbtot += popInfo->posiBP[rangocromo[conta+1]-1]-popInfo->posiBP[rangocromo[conta]];
  }
  popInfo->Mtot /= 100.0;
  popInfo->Mbtot /= 1000000.0;
  if (popInfo->numLoci != contalines) {
    std::cerr << "Different number of loci in ped and map files" << std::endl;
    return false;
  }
  popInfo->numCromo = popInfo->numCromo;
  return true;
}

// Read a .ped genotype file together with its companion .map file into
// popInfo (calls ReadPed then ReadMap).
// Params: fichPed — path to the .ped file; fichMap — path to the .map file;
//   popInfo — population struct filled in place.
// Returns: true on success, false if either file fails to load.
bool ReadFile(std::string fichPed, std::string fichMap,
              PopulationInfo* popInfo) {
  if (!ReadPed(fichPed, popInfo)) {
    return false;
  }
  if (!ReadMap(fichMap, popInfo)) {
    return false;
  }
  return true;
}
