
#include "tped.hpp"
#include <ostream>

namespace {
struct TpedAlleleTable {
  bool ok[256];
  constexpr TpedAlleleTable() : ok{} {
    for (const char* p = "AGCTNagctn0123456789"; *p; ++p) {
      ok[static_cast<unsigned char>(*p)] = true;
    }
  }
};
const TpedAlleleTable kTpedAllele;
// Test whether a character is an accepted .tped allele symbol.
// Params: c — allele character.
// Returns: true if c is a valid allele code, else false.
inline bool IsValidTpedAllele(char c) {
  return kTpedAllele.ok[static_cast<unsigned char>(c)];
}
}

// Read a PLINK .tped file into popInfo: genotype matrix, SNP positions,
// chromosomes, and the genetic map carried in the .tped columns.
// Params: fichTped — path to the .tped file; popInfo — population struct
//   filled in place (loci, individuals, positions, genotypes).
// Returns: true on success, false if the file can't be read or is malformed.
bool ReadTped(std::string fichTped, PopulationInfo *popInfo) {
  char base1, base2;
  int contaIndBase = 0, nline = 0;
  int conta = 0, posi = 0, posi2 = 0, longi = 0;
  std::string line;
  std::string cromocod;
  std::string cromocodback = "laksjhbqne";
  std::vector<int>    rangocromo(MAXCROMO);
  std::vector<char>   base(MAXLOCI);
  std::vector<double> chrnum(MAXCROMO);
  std::vector<double> chrprop(MAXCROMO);

  std::ifstream entrada;
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

    ++contalines;

    ++posi;
    ++conta;

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
    posi2 = posi;
    posi = static_cast<int>(line.find_first_of(" \t", posi2));
    double position_centimorgans = std::stod(line.substr(posi2, posi-posi2));
    popInfo->posiCM[nLoci] = position_centimorgans;
    ++posi;
    ++conta;
    posi2 = posi;
    posi = static_cast<int>(line.find_first_of(" \t", posi2));
    int position_basepair = std::stoi(line.substr(posi2, posi-posi2));
    popInfo->posiBP[nLoci] = position_basepair;

    ++posi;
    ++conta;
    if (conta == 4) {
      if (nLoci >= MAXLOCI) {
        std::cerr << "Reached max number of loci (" << MAXLOCI << ")"
                  << std::endl;
        return false;
      }
      nIndi = 0;
      while (posi < longi) {
        base1 = line.at(posi);
        posi2 = posi;
        posi = static_cast<int>(line.find_first_of(" \t", posi2));
        if (posi < 0) {
          break;
        }
        ++posi;
        base2 = line.at(posi);
        if (!IsValidTpedAllele(base1) || !IsValidTpedAllele(base2)) {
          std::cerr << "Wrong allele in line " << nline
                    << " of tped file (maybe also in other lines). Only "
                       "bases A, G, C, T and N (case-insensitive) and "
                       "numbers from 0 to 9 are allowed."
                    << std::endl;
          exit(EXIT_FAILURE);
        }

        if ((popInfo->haplotype == 2 ? nIndi + 1 : nIndi) >= MAXIND) {
          std::cerr << "Reached limit of sample size (" << MAXIND << ")"
                    << std::endl;
          return false;
        }
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
            if (base1 == base2) {
              popInfo->indi[nIndi][nLoci] = (base1 == base[nLoci] ? 0 : 2);
            } else {
              popInfo->indi[nIndi][nLoci] = 1;
            }
          } 
          else if (popInfo->haplotype ==1){
            if (base1 == base[nLoci]) {
              popInfo->indi[nIndi][nLoci] = 0;
            } else {
              popInfo->indi[nIndi][nLoci] = 2;
            }
          }
          else {
            popInfo->indi[nIndi][nLoci] =
              base1 == base[nLoci] ? 0 : 2;
            popInfo->indi[nIndi + 1][nLoci] =
              base2 == base[nLoci] ? 0 : 2;
          }
        }
        else {
          popInfo->indi[nIndi][nLoci] = 9;
          if (popInfo->haplotype == 2) {
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
        if (popInfo->haplotype == 2){
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
  popInfo->Mtot /= 100.0;
  popInfo->Mbtot /= 1000000.0;
  return true;
}
