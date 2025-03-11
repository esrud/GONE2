
#include "ped_map.hpp"
#include "./pool.hpp"

bool ReadPed(std::string fichPed, PopulationInfo* popInfo) {
  /*
   * Takes as input the file name and a pointer to the population
   */
  // Contiene las bases de referencia de cada locus
  uint8_t base[MAXLOCI] = {0};
  uint8_t base1[1], base2[1];
  int contaLociBase = 0, nline=0;
  int conta = 0, posi = 0, posi2 = 0, longi = 0;
  std::string line;
  std::string str_plink = "AGCTNagctn0123456789";
  std::string str_base = "";
  bool rightletter;

  // READING .ped DATA:
  std::ifstream entrada;
  entrada.open(fichPed, std::ios::in);  // Bucle de lectura del fichero ped
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
      while (posi < longi) {  // asigna genot.
        base1[0] = line.at(posi);
        posi2 = posi;
        posi=static_cast<int>(line.find_first_of(" \t",posi2));
        if (posi < 0) {
          break;
        }
        ++posi;
        base2[0] = line.at(posi);
        rightletter=true;
        str_base=base1[0];
        if (str_plink.find(str_base) == std::string::npos){rightletter=false;}
        str_base=base2[0];
        if (str_plink.find(str_base) == std::string::npos){rightletter=false;}
        if (!rightletter) {
            std::cerr << "Wrong allele in line "<<nline<<" of ped file (maybe also in other lines). Only bases A, G, C, T and N (case-insensitive) and numbers from 0 to 9 are allowed." << std::endl;
            exit(EXIT_FAILURE);
        }
        if ((base1[0]!='0') && (base2[0]!='0') && (base1[0]!='N') && (base2[0]!='N') && (base1[0]!='n') && (base2[0]!='n')) {
          if (base[popInfo->numLoci] == '\0') {
            base[popInfo->numLoci] = base1[0];
          }
          if (base1[0] != base[popInfo->numLoci])
          {
              base1[0] = 'X';
          }
          if (base2[0] != base[popInfo->numLoci])
          {
              base2[0] = 'X';
          }

          if ((popInfo->haplotype == 0) || (popInfo->haplotype == 3)){ // FASE DESCONOCIDA y Low coverage
            if (base1[0] == base2[0]) {
                if (base1[0] == base[popInfo->numLoci]) {
                  popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 0;
                } else {
                  popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 2;
                }
            } 
            else {
                    popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 1;
            }
          }
          else if (popInfo->haplotype == 1){ // HAPLOIDES
            if (base1[0] == base[popInfo->numLoci]) {
              popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 0;
            } else {
              popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 2;
            }
          }
          else{ // FASE CONOCIDA (2) genera dos individuos haploides
              if (base1[0] == base[popInfo->numLoci]) {
                popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 0;
              } 
              else {
                popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 2;
              }
              if (base2[0] == base[popInfo->numLoci]) {
                popInfo->indi[popInfo->numIndi+1][popInfo->numLoci] = 0;
              } 
              else {
                popInfo->indi[popInfo->numIndi+1][popInfo->numLoci] = 2;
              }
          }
        } 
        else {
          // '9' = Genotipo sin asignar
          if (popInfo->haplotype !=2){// Fase desconocida, haploides y low coverage
              popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 9;
          }
          else{// FASE CONOCIDA (2) genera dos individuos haploides
              popInfo->indi[popInfo->numIndi][popInfo->numLoci] = 9;
              popInfo->indi[popInfo->numIndi+1][popInfo->numLoci] = 9;
          }
        }
        posi2 = posi;
        posi=static_cast<int>(line.find_first_of(" \t",posi2));
        if (posi < 0) {
          posi = longi;
        }
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
          //break;
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
      if (popInfo->haplotype == 2){// Fase Conocida son el doble de individuos
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

bool ReadMap(std::string fichMap, PopulationInfo* popInfo) {
  int posi = 0, posi2 = 0, longi = 0;
  std::string line;

  std::string cromocod;
  std::string cromocodback = "laksjhbqne";

  std::ifstream entrada;
  entrada.open(fichMap, std::ios::in);  // Bucle de lectura del fichero ped
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

    ++posi;  // Nombre del SNP que no se lee
    posi2 = posi;
    posi=static_cast<int>(line.find_first_of(" \t",posi2));
    if (posi <= 0) {
      std::cerr << "Error in map file (1)" << std::endl;
      return false;
    }

    ++posi;  // Localizacion cM
    posi2 = posi;
    posi=static_cast<int>(line.find_first_of(" \t",posi2));
    if (posi <= 0) {
      std::cerr << "Error in map file (2)" << std::endl;
      return false;
    }
    popInfo->posiCM[contalines] = std::stod(line.substr(posi2, posi - posi2));
    ++posi;

    if (posi > longi) {  // Localizacion bp
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
  popInfo->Mtot /= 100.0;       // en Morgans
  popInfo->Mbtot /= 1000000.0;  // en megabases
  if (popInfo->numLoci != contalines) {
    std::cerr << "Different number of loci in ped and map files" << std::endl;
    return false;
  }
  popInfo->numCromo = popInfo->numCromo;
  return true;
}

bool ReadFile(std::string fichPed, std::string fichMap,
              PopulationInfo* popInfo) {
  /*
   * Takes as input the file name and a pointer to the population
   * matrix and returns the number of individuals in the file
   */

  //std::cout << "Reading PED file" << std::endl;
  if (!ReadPed(fichPed, popInfo)) {
    return false;
  }
  //std::cout << "Reading MAP file" << std::endl;
  if (!ReadMap(fichMap, popInfo)) {
    return false;
  }
  return true;
}
