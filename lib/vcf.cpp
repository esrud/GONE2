#include "vcf.hpp"

namespace {
struct VcfAlleleTable {
  bool ok[256];
  constexpr VcfAlleleTable() : ok{} {
    for (const char* p = "AGCT,agct"; *p; ++p) {
      ok[static_cast<unsigned char>(*p)] = true;
    }
  }
};
const VcfAlleleTable kVcfAllele;
// Test whether a character is an accepted VCF allele symbol.
// Params: c — allele character.
// Returns: true if c is a valid allele code, else false.
inline bool IsValidVcfAllele(char c) {
  return kVcfAllele.ok[static_cast<unsigned char>(c)];
}
}

// Read genotype data from a VCF file into popInfo (genotype matrix, SNP
// positions, chromosomes).
// Params: fichVcf — path to the .vcf file; popInfo — population struct
//   filled in place.
// Returns: true on success, false if the file can't be read or is malformed.
bool ReadVcf(std::string fichVcf, PopulationInfo *popInfo) {
  char base1, base2;
  int contaIndBase = 0, nline = 0, i, j;
  int conta = 0, posi = 0, posi2 = 0, longi = 0;
  std::string line;
  std::string cromocod;
  std::string cromocodback = "laksjhbqne";
  std::vector<int>    rangocromo(MAXCROMO);
  std::vector<char>   base(MAXLOCI);
  std::vector<double> chrnum(MAXCROMO);
  std::vector<double> chrprop(MAXCROMO);
  bool hayalelos;

  std::ifstream entrada;
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

    posi2 = posi;
    posi  = int(line.find_first_of("\t", posi2));
    j = posi - posi2;
    bool rightletter = true;
    hayalelos = true;
    for (i = 0; i < j; ++i) {
      const char c = line[posi2 + i];
      if (c == '.') {
        hayalelos = false;
      } else if (!IsValidVcfAllele(c)) {
        rightletter = false;
        break;
      }
    }
    if (!rightletter) {
      std::cerr << "Wrong reference allele in line " << nline
                << " of vcf file. Only bases A, G, C and T "
                   "(case-insensitive) are allowed."
                << std::endl;
      exit(EXIT_FAILURE);
    }
    ++posi;
    ++conta;

    posi2 = posi;
    posi  = int(line.find_first_of("\t", posi2));
    j = posi - posi2;
    for (i = 0; i < j; ++i) {
      const char c = line[posi2 + i];
      if (c == '.') {
        hayalelos = false;
      } else if (!IsValidVcfAllele(c)) {
        rightletter = false;
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
    bool haygenotipos = line.substr(posi, 2) == "GT";

    if (haygenotipos && hayalelos) {
      if (cromocod != cromocodback) {
        cromocodback = cromocod;
        rangocromo[ncromos] = contalines;
        ++ncromos;
      }
      popInfo->cromo[contalines] = ncromos;
      ++contalines;
      posi2 = posi;
      posi = static_cast<int>(line.find_first_of("\t", posi2));
      if (posi < 0) {
        std::cerr << "Line too short in vcf file" << std::endl;
        return false;
      }
      ++posi;
      ++conta;
      if (nLoci >= MAXLOCI) {
        std::cerr << "Reached max number of loci (" << MAXLOCI << ")"
                  << std::endl;
        return false;
      }
      nIndi = 0;
      while (posi < longi) {
        base1 = line.at(posi);
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
        if ((popInfo->haplotype == 2 ? nIndi + 1 : nIndi) >= MAXIND) {
          std::cerr << "Reached limit of sample size (" << MAXIND << ")"
                    << std::endl;
          return false;
        }
        if ((base1 != '.') && (base2 != '.')) {
          if (base[nLoci] == '\0') {
            base[nLoci] = base1;
          }
          base1 = base1 != base[nLoci] ? 'X' : base[nLoci];
          base2 = base2 != base[nLoci] ? 'X' : base[nLoci];

          if ((popInfo->haplotype == 0) || (popInfo->haplotype == 3)) {
            if (base1 == base2) {
              popInfo->indi[nIndi][nLoci] =
                base1 == base[nLoci] ? 0 : 2;
            } else {
              popInfo->indi[nIndi][nLoci] = 1;
            }
          } 
          else if ((popInfo->haplotype == 1)) {
            popInfo->indi[nIndi][nLoci] = (base1 == base[nLoci] ? 0 : 2);
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
        posi = static_cast<int>(line.find_first_of("\t", posi2));
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
  popInfo->Mbtot = 0;
  for (int i = 0; i < popInfo->numCromo; ++i) {
    int idx = rangocromo[i];
    int next_idx = rangocromo[i+1] - 1;
    popInfo->Mbtot += popInfo->posiBP[next_idx] - popInfo->posiBP[idx];
  }
  popInfo->Mbtot /= 1000000.0;
  return true;
}
