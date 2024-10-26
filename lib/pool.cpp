
#include "./pool.hpp"

void SetMutationParameters(Pool* pool, double frecMut, double efectoMut,
                           double efectoMutLateral, double frecRec,
                           double frecNoRec, double frecMutLateral,
                           int maxSegmentos, bool solape) {
  /*
   * Simple wrapper to set the mutation parameters 
   */
  pool->mutParams.frecMut = frecMut;
  pool->mutParams.efectoMut = efectoMut;
  pool->mutParams.efectoMutLateral = efectoMutLateral;
  pool->mutParams.frecRec = frecRec;
  pool->mutParams.frecNoRec = frecNoRec;
  pool->mutParams.frecMutLateral = frecMutLateral;
  pool->mutParams.maxSegmentos = maxSegmentos;
  pool->mutParams.solape = solape;
}


void SetInitialPoolParameters(Pool* pool, double mincval) {
  /*
   * Sets the parameters that will be needed by the rest of the functions
   */

  // Set gmax
  //std::cout << "mincval= " << mincval << std::endl;
  int gmax = std::min(static_cast<int>(1.0 / mincval), kNumGenMax-10);
  //if (gmax > (kNumGenMax - 10)) {
  //  gmax = kNumGenMax - 10;
  //}
  int gmax2 = gmax * 26 / 40;
  int topeposigen = gmax2 - 2 * kResolucion;
  pool->poolParams.gmax[0] = gmax;
  pool->poolParams.gmax[1] = gmax2;
  pool->poolParams.gmax[2] = gmax2 + (gmax - gmax2) / 2;

  // Set topePosInicio
  pool->poolParams.topePosInicio[0] = std::min(topeposigen, kTopePosInicio1);
  pool->poolParams.topePosInicio[1] = std::min(topeposigen, kTopePosInicio2);
  pool->poolParams.topePosInicio[2] = std::min(topeposigen, kTopePosInicio3);
  pool->poolParams.topePosInicio[3] = std::min(topeposigen, kTopePosInicio4);

  // Set number of generations to run for
  pool->numGenerations = kNumGenerations;

  // Set the number of children per generation
  // TODO(me): load this from a file or something
  pool->numChildren = kNumDes;

  // Set initial mutation parameters
  // TODO(me): load this from a file or something
  pool->mutParams.frecMut = 0.3;
  pool->mutParams.efectoMut = 0.3;
  pool->mutParams.efectoMutLateral = 5.0;
  pool->mutParams.frecRec = 0.6;
  pool->mutParams.frecNoRec = 0.3;
  pool->mutParams.frecMutLateral = 0.2;
  pool->mutParams.maxSegmentos = 4;
  pool->mutParams.solape = true;
}

void PrePopulatePool(Pool* pool, GsampleInfo* sampleInfo) {
  /*
   * Populates the pool with the kNumDesInicio bichos, using an
   * educated guess based on the sample
   */

  // PRECARGA ALEATORIA CON SELECCION
  Bicho bicho = {0};

  bicho.efval = sampleInfo->fVal;      //++
  bicho.nSeg = 4;          //++
  bicho.segBl[0] = 0;      //++
  bicho.segBl[3] = pool->poolParams.gmax[2];  //++
  bicho.segBl[4] = pool->poolParams.gmax[0];   //++

  double r;
  double tope, Ne, NeCoeff;
  int posi1, posi2;

  int i, j;
  double maxSC, minSC;
  int idxMaxSC, idxMinSC;

  for (i = 0; i < kNumDesInicio; ++i) {
    bool isValidChild = false;
    while (!isValidChild) {
      r = xprng.uniforme01();
      if (r < 0.4) {
        posi1 = static_cast<int>(
                    xprng.uniforme01() *
                    (pool->poolParams.topePosInicio[1] - kResolucion)) +
                kResolucion;
        posi2 = static_cast<int>(xprng.uniforme01() *
                                 (pool->poolParams.topePosInicio[1] - kResolucion)) +
                kResolucion;
      } else if (r < 0.8) {
        posi1 = static_cast<int>(xprng.uniforme01() *
                                 (pool->poolParams.topePosInicio[2] - kResolucion)) +
                kResolucion;
        posi2 = static_cast<int>(xprng.uniforme01() *
                                 (pool->poolParams.topePosInicio[2] - kResolucion)) +
                kResolucion;
      } else {
        posi1 = static_cast<int>(xprng.uniforme01() *
                                 (pool->poolParams.topePosInicio[3] - kResolucion)) +
                kResolucion;
        posi2 = static_cast<int>(xprng.uniforme01() *
                                 (pool->poolParams.topePosInicio[3] - kResolucion)) +
                kResolucion;
      }

      if (posi1 > posi2) {
        std::swap(posi1, posi2);
        //int posi3 = posi1;
        //posi1 = posi2;
        //posi2 = posi3;
      }
      if ((posi2 - posi1) > MIN_NE_SIZE) {
        isValidChild = true;
      }
    }
    bicho.segBl[1] = posi1;  //++
    bicho.segBl[2] = posi2;  //++
    bicho.NeBl[0] = sampleInfo->NeMed;
    tope = 2 * xprng.uniforme01() * kTopeSalto;
    // tope=topesalto;  // ALTERNATIVA A LA LINEA ANTERIOR
    for (j = 1; j < 3; ++j) {
      NeCoeff = 1.0 + xprng.uniforme01() * tope;  // ERA2
      if (xprng.uniforme01() < 0.66) {
        NeCoeff = 1.0 / NeCoeff;
      }  // ERA0.5
      
      //Ne = sampleInfo->NeMed * NeCoeff;
      Ne = std::max(std::min(sampleInfo->NeMed * NeCoeff, MAX_NE_SIZE), MIN_NE_SIZE);
      
      bicho.NeBl[j] = Ne;  //++
    }
    NeCoeff = 1.0 + xprng.uniforme01() * tope / 2.0;  // era 0.5
    //Ne = sampleInfo->NeMed * NeCoeff;
    Ne = std::max(std::min(sampleInfo->NeMed * NeCoeff, MAX_NE_SIZE), MIN_NE_SIZE);
    
    bicho.NeBl[3] = Ne;  //++
    bicho.SCval = CalculaSC(&bicho, sampleInfo);
    // Selecciona los padres
    if (i < kTercio12) {
      pool->parents[i] = bicho;
    } else {
      // Busca el peor SC entre los futuros padres
      maxSC = 0;
      idxMaxSC = 0;
      for (j = 0; j < kTercio12; ++j) {
        if (pool->parents[j].SCval > maxSC) {
          idxMaxSC = j;
          maxSC = pool->parents[j].SCval;
        }
      }
      if (bicho.SCval < pool->parents[idxMaxSC].SCval) {
        pool->parents[idxMaxSC] = bicho;
      }
    }
  }

  // Ordena los padres resultantes
  // TODO(me): Maybe quicksort the sheit out of this
  for (i = 0; i < kTercio12 - 1; ++i) {
    idxMinSC = i;
    minSC = MAX_DOUBLE;
    for (j = i + 1; j < kTercio12; ++j) {
      if (pool->parents[j].SCval < minSC) {
        idxMinSC = j;
        minSC = pool->parents[j].SCval;
      }
    }
    if (idxMinSC != i) {  // swap
      bicho = pool->parents[i];
      pool->parents[i] = pool->parents[idxMinSC];
      pool->parents[idxMinSC] = bicho;
    }
  }
}

bool CheckParentMutations(Bicho* parent, const int posigen) {
  /*
   * Assert that the block resolution is valid after a mutation
   */
  int abs_val;
  int temp;
  for (int j = 0; j <= parent->nSeg; ++j) {
    abs_val  = (parent->segBl[j] - posigen);
    // A faster way to get the absolute value of an int of 32 bits
    // Does not work for double or float variables
    // Make a mask of the sign bit
    temp = abs_val >> 31;
    // Toggle the bits if the value is negative
    abs_val ^= temp;
    // Add one if it was negative
    abs_val += temp & 1;
    //if (abs(parent->segBl[j] - posigen) <= kResolucion) {
    if(abs_val <= kResolucion) {
      return false;
    }
  }
  return true;
}

bool MutateFusion(Bicho* bicho, Pool* pool) {
  /*
   * Fuse three? adjancent blocks together
   */
  int posiblock;
  double ancho,ancho1, ancho2;
  double anchoDivNe;
  if (bicho->nSeg > 3) {
    if ((xprng.uniforme01() < pool->mutParams.frecNoRec) ||
        (bicho->nSeg > pool->mutParams.maxSegmentos)) {
      // FUSIONA AL AZAR
      posiblock = static_cast<int>(xprng.uniforme01() * (bicho->nSeg - 2));

      ancho1 = bicho->segBl[posiblock + 1] - bicho->segBl[posiblock];
      ancho2 = bicho->segBl[posiblock + 2] - bicho->segBl[posiblock + 1];

      ancho = (ancho1 + ancho2) / (ancho1 / bicho->NeBl[posiblock] +
                                ancho2 / bicho->NeBl[posiblock + 1]);
      if (posiblock > 0) {
        anchoDivNe = ancho / bicho->NeBl[posiblock - 1];
        if (anchoDivNe < kInvTopeSalto2) {
          return false;
        }
      }

      if (posiblock < bicho->nSeg - 2) {
        anchoDivNe = ancho / bicho->NeBl[posiblock + 2];
        if (anchoDivNe > kTopeSalto2) {
          return false;
        }
      }

      bicho->NeBl[posiblock] = ancho;
      for (int j = posiblock + 1; j <= bicho->nSeg; ++j) {
        bicho->NeBl[j] = bicho->NeBl[j + 1];
        bicho->segBl[j] = bicho->segBl[j + 1];
      }
      --bicho->nSeg;
    }
  }
  // Control del tope de número de segmentos
  if ((bicho->nSeg > pool->mutParams.maxSegmentos) || (bicho->nSeg < 3)) {
    return false;
  }
  return true;
}

void Reproduce(Pool* pool, GsampleInfo* sampleInfo, int childIdx) {
  /*
   * Main reproduction block to generate a new individual
  */
  int ind1, ind2;
  //int ind3;
  int topePosInicio;
  double r;
  int posigen;
  int nSeg;
  bool isValidChild=false;
  Bicho child = {0};
  Bicho* parent1;
  Bicho* parent2;
  // Counters
  int j;

#if FLAG_POR_BLOQUES
  int posiblock;
#endif

  while (!isValidChild) {
    //isValidChild = true;  // var bool (b) declarada global (G). Tendría que
                          // hacerse copia para cada thread (T)
    // ind1: iGT (int Global Thread), tercio1:
    // iGnt (G no-t)
    ind1 = static_cast<int>(xprng.uniforme01() * kTercio1);
    // ind2: iGT , tercio12: iGnT
    ind2 = static_cast<int>(xprng.uniforme01() * kTercio12);
    if (xprng.uniforme01() < 0.5) {
      std::swap(ind1, ind2);
    }  // genera: dGT (double GT), ind3: iGT

    parent1 = &pool->parents[ind1];
    parent2 = &pool->parents[ind2];

    // Recombinacion ENTRE y DENTRO CON IGUAL PROBABILIDAD: genera una
    // serie a partir de dos o un padre
    if (xprng.uniforme01() < pool->mutParams.frecRec) {  // frecrec: dGnt
      r = xprng.uniforme01();                            // aa: dGT
      if (r < 0.1) {
        // topeposiinicio: iGT, topeposiinicio1: iGnT
        topePosInicio = pool->poolParams.topePosInicio[0];
      } else if (r < 0.5) {
        // topeposiinicio2: iGnT
        topePosInicio = pool->poolParams.topePosInicio[1];
      } else if (r < 0.9) {
        // topeposiinicio3: iGnT
        topePosInicio = pool->poolParams.topePosInicio[2];
      } else {
        // topeposiinicio4: iGnT
        topePosInicio = pool->poolParams.topePosInicio[3];
      }

#if FLAG_POR_BLOQUES
      do {
        // posiblock: iGT, bichoP: structGnT
        posiblock = static_cast<int>(xprng.uniforme01() * (parent1->nSeg - 1));
        // posigen: iGT, topeposiinicio: iGT
        posigen =
            parent1->segBl[posiblock] +
            static_cast<int>(xprng.uniforme01() * (parent1->segBl[posiblock + 1] -
                                                parent1->segBl[posiblock]));
      } while (posigen > topePosInicio);
#else
      // resolucion: iGnT
      posigen =
          kResolucion + static_cast<int>(xprng.uniforme01() * topePosInicio);
#endif
      // DENTRO: iguala individuos y centra posigen
      if (xprng.uniforme01() < 0.5) {
        parent2 = parent1;  // ind1: iGT, ind2: iGT
      }
      // Check that the mutations resulted in valid blocks
      // for parent1
      isValidChild = CheckParentMutations(parent1, posigen);
      if (!isValidChild) {
        continue;
      }

      // Check that the mutations resulted in valid blocks
      // for parent2
      isValidChild = CheckParentMutations(parent2, posigen);
      if (!isValidChild) {  // flag bGT
        continue;
      }
      nSeg = 0;

      // Inherit from parent 1
      for (j = 0; j < parent1->nSeg; ++j) {
        if (parent1->segBl[j] < posigen) {
          child.segBl[nSeg] = parent1->segBl[j];
          child.NeBl[nSeg] = parent1->NeBl[j];
          ++nSeg;
        } else {
          break;
        }
      }

      // Inherit from parent 2
      for (j = 1; j < parent2->nSeg + 1; ++j) {
        if (parent2->segBl[j] > posigen) {
          child.segBl[nSeg] = posigen;
          child.NeBl[nSeg] = parent2->NeBl[j - 1];
          ++nSeg;
          break;
        }
      }

      for (j = 1; j < parent2->nSeg; ++j) {
        if (parent2->segBl[j] > posigen) {
          child.segBl[nSeg] = parent2->segBl[j];
          child.NeBl[nSeg] = parent2->NeBl[j];
          ++nSeg;
        }
      }
      child.nSeg = nSeg;
      child.segBl[nSeg] = parent2->segBl[parent2->nSeg];
    } else {
      child = *parent1;
    }

    isValidChild = MutateFusion(&child, pool);
  }

  MutateNeRnd(&child, pool);
  MutateLateral(&child, pool);

  // Get the offspring score
  child.efval = parent1->efval;
  child.SCval = CalculaSC(&child, sampleInfo);

  // If the pool isn't full, add a new child
  // else, replace the worst child if the offspring is better
  if (childIdx < kNHijos) {
    pool->children[childIdx] = child;
  } else {
    double maxSC = -1;
    int idxMaxSC = 0;
    for (int i = 0; i < childIdx; ++i) {
      if (pool->children[i].SCval > maxSC) {
        maxSC = pool->children[i].SCval;
        idxMaxSC = i;
      }
    }
    // TODO(me): Shouldn't you also add the children that
    //           perform worse?
    if (child.SCval < maxSC) {
      pool->children[idxMaxSC] = child;
    }
  }
}

void MutateLateral(Bicho* bicho, Pool* pool) {
  /*
   * Mutación lateral: Cambia el límite de segmentos al azar
   */
  int posiblock;
  int efecto;

  for (posiblock = 1; posiblock < bicho->nSeg - 1; ++posiblock) {
    if (xprng.uniforme01() < pool->mutParams.frecMutLateral) {
      efecto = static_cast<int>(xprng.uniforme01() * pool->mutParams.efectoMutLateral) + 1;
    }  // efecto: dGT
    else {
      efecto = 1;
      if (xprng.uniforme01() < 0.5) {
        efecto = 0;
      }
    }
    // ERA0.6
    if (xprng.uniforme01() < 0.5) {
      if ((bicho->segBl[posiblock] - bicho->segBl[posiblock - 1]) >
          (kResolucion + efecto)) {
        bicho->segBl[posiblock] -= efecto;
      }
    } else {
      if ((bicho->segBl[posiblock + 1] - bicho->segBl[posiblock]) >
          (kResolucion + efecto)) {
        bicho->segBl[posiblock] += efecto;
      }
    }
  }
}

void MutateNeRnd(Bicho* bicho, Pool* pool) {
  /* 
   * Mutación: Cambia el valor de ne de segmentos al azar
   */
  int posiblock;
  double efecto;
  double NeMutado;
  double mutatedNeRate;
  int lastIdx, prevToLastIdx;

  for (posiblock = 0; posiblock < bicho->nSeg; ++posiblock) {
    if (xprng.uniforme01() < pool->mutParams.frecMut) {
      efecto = bicho->NeBl[posiblock] * xprng.uniforme01() * pool->mutParams.efectoMut;
    } else {
      efecto = bicho->NeBl[posiblock] * xprng.uniforme01() * kEfectoMutSuave;
    }

    if (xprng.uniforme01() < 0.5) {
      efecto = -efecto;
    }

    NeMutado = bicho->NeBl[posiblock] + efecto;  // aa: dMT

    // TODO(me): this if can be taken out of the loop at the expense of
    //           duplicating code
    if (posiblock > 0) {
      mutatedNeRate = NeMutado / bicho->NeBl[posiblock - 1];
      if (mutatedNeRate < kInvTopeSalto2) {
        continue;
      }
    }

    if (posiblock < (bicho->nSeg - 1)) {
      mutatedNeRate =
          NeMutado /
          bicho->NeBl[posiblock + 1];  // bicho: structGT, posiblock: iGT
      if (mutatedNeRate > kTopeSalto2) {
        continue;
      }
    }

    bicho->NeBl[posiblock] = std::max(std::min(NeMutado, MAX_NE_SIZE), MIN_NE_SIZE);
  }
  // Iguala dos dos ultimos:
  lastIdx = bicho->nSeg - 1;
  prevToLastIdx = bicho->nSeg - 2;
  if (bicho->NeBl[lastIdx] > (bicho->NeBl[prevToLastIdx] * 1.2)) {
    bicho->NeBl[lastIdx] = bicho->NeBl[prevToLastIdx] * 1.2;
  } else if (bicho->NeBl[lastIdx] < (bicho->NeBl[prevToLastIdx] / 1.2)) {
    bicho->NeBl[lastIdx] = bicho->NeBl[prevToLastIdx] / 1.2;
  }
}

void ReproducePool(Pool* pool, GsampleInfo* sampleInfo) {
  /*
   * Generate as many individuals as set on the pool from
   * an array of parents
   */
  bool foundBetterChild;
  double minSC;
  int idxMinSC;
  int i, j;

  for (int des = 0; des < pool->numChildren; ++des) {
    Reproduce(pool, sampleInfo, des);
  }

  // El primer tercio de los bichoP no se toca.
  if (pool->mutParams.solape) {
    // Mete los mejores hijos en el segundo tercio de los padres si es que
    // son mejores
    for (i = kTercio1; i < kTercio12; ++i) {
      // Busca el mínimo SC entre los hijos
      foundBetterChild = false;
      minSC = MAX_DOUBLE;
      idxMinSC = 0;
      for (j = 0; j < kNHijos; ++j) {
        if (pool->children[j].SCval < minSC) {
          minSC = pool->children[j].SCval;
          idxMinSC = j;
        }
      }
      // mira si ese mínimo es mejor que algún padre del
      // segundo tercio
      for (j = kTercio1; j < kTercio12; ++j) {
        if (minSC < pool->parents[j].SCval) {
          pool->parents[j] = pool->children[idxMinSC];
          pool->children[idxMinSC].SCval = MAX_DOUBLE;
          foundBetterChild = true;
          break;
        }
      }
      if (!foundBetterChild) {
        break;
      }
    }
  } else {
    // COPIA LOS HIJOS EN LOS PADRES:
    // TODO(me): this could be faster using pointers instead of copying
    //           so after each generation parents points to the array of
    //           children at generation n-1, and children points to the
    //           array of parents at generation n-1
    for (i = 0; i < kTercio12; ++i) {
      pool->parents[i] = pool->children[i];
    }
  }

  // Inversión: Intercambia dos ultimosNe adyacentes con probabilidad
  //            kFrecInversion
  Bicho* bicho;
  Bicho tempBicho;
  int posiblock;
  double NeRate, temp;
  if (xprng.uniforme01() < kFrecInversion) {
    for (i = 0; i < kTercio12; ++i) {
      bicho = &pool->parents[i];
      if (bicho->nSeg >= 4) {
        posiblock = bicho->nSeg - 2;
        NeRate = bicho->NeBl[posiblock + 1] / bicho->NeBl[posiblock - 1];
        if (NeRate > kInvTopeSalto && NeRate < kTopeSalto) {
          if (bicho->NeBl[posiblock + 1] != bicho->NeBl[posiblock]) {
            temp = bicho->NeBl[posiblock + 1];
            bicho->NeBl[posiblock + 1] = bicho->NeBl[posiblock];
            bicho->NeBl[posiblock] = temp;
            bicho->SCval = CalculaSC(bicho, sampleInfo);
          }
        }
      }
    }
  }

  // Ordena a los padres
  for (i = 0; i < kTercio12 - 1; ++i) {
    idxMinSC = i;
    for (j = i + 1; j < kTercio12; ++j) {
      if (pool->parents[j].SCval < pool->parents[idxMinSC].SCval) {
        idxMinSC = j;
      }
    }
    if (idxMinSC != i) {  // swap
      tempBicho = pool->parents[i];
      pool->parents[i] = pool->parents[idxMinSC];
      pool->parents[idxMinSC] = tempBicho;
    }
  }
}

void RunDbg(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal) {
  double SCmed, nSegMed;

  std::ofstream salida;
  std::string fichero_sal_NeH = fichSal + "_GONE2_Nebest";
  std::string fichero_sal_d2 = fichSal + "_GONE2_d2";
  std::string fichero_sal_log = fichSal + "_GONE2_log";
  std::string fichero_sal_rearrange = fichSal + "_GONE2_input";
  std::string fichero_sal_evol = fichSal + "_GONE2_evol";  // --debug
  std::string fichero_dbg = fichSal + "_GONE2_dbg";

  salida.open(fichero_dbg, std::ios::app);
  for (int j=0;j < pool->parents[0].nSeg;++j){
      salida << j << "\t" <<pool->parents[0].segBl[j] <<"\t" << pool->parents[0].NeBl[j]/2.0 <<"\n";}
  salida<<"\n";
  salida.close();

  for (int gen = 0; gen < pool->numGenerations; ++gen) {
    if (gen == static_cast<int>(static_cast<double>(gen) / 100.0)) {
      salida.open(fichero_dbg, std::ios::app);
      for (int j = 0; j < pool->parents[0].nSeg; ++j) {
        salida << "-> " << gen << "\t" << j << "\t" << pool->parents[0].segBl[j]
               << "\t" << pool->parents[0].NeBl[j] / 2.0 << "\n";
      }
      salida << "\n";
      salida.close();
    }
    switch (gen) {
      case 300:
        SetMutationParameters(pool, 0.2, 0.2, 2.0, 0.2, 0.5, 0.2,
                              pool->mutParams.maxSegmentos + 2, true);
        break;
      case 600:
        SetMutationParameters(pool, 0.5, 0.05, 1.0, 0.2, 0.5, 0.2,
                              pool->mutParams.maxSegmentos + 10, true);
        break;
      case 700:
        SetMutationParameters(pool, 0.5, 0.2, 1.0, 1.0, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 10, true);
        break;
      case 710:
        SetMutationParameters(pool, 0.5, 0.2, 1.0, 1.0, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 10, false);
        break;
      case 720:
        SetMutationParameters(pool, 0.5, 0.1, 1.0, 0.95, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 10, false);
        break;
      case 730:
        SetMutationParameters(pool, 0.5, 0.04, 1.0, 0.95, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 20, false);
        break;
      default:
        break;
    }
    ReproducePool(pool, sampleInfo);
    SCmed = nSegMed = 0;
    for (int i = 0; i < kTercio12; ++i) {
      SCmed += pool->parents[i].SCval;
      nSegMed += pool->parents[i].nSeg;
    }
    SCmed /= kTercio12;
    nSegMed /= kTercio12;

    salida.open(fichero_sal_evol, std::ios::app);
    salida << gen << "\t" << pool->parents[0].SCval << "\t" << SCmed << "\t"
           << pool->parents[0].nSeg << "\t" << nSegMed << "\n";
    salida.close();
  }

  salida.open(fichero_dbg, std::ios::app);
  for (int j = 0; j < pool->parents[0].nSeg; ++j) {
    salida << j << "\t" << pool->parents[0].segBl[j] << "\t"
           << pool->parents[0].NeBl[j] / 2.0 << "\n";
  }
  salida.close();
}

void Run(Pool* pool, GsampleInfo* sampleInfo, std::string fichSal) {
  for (int gen = 0; gen < pool->numGenerations; ++gen) {
    switch (gen) {
      case 300:
        SetMutationParameters(pool, 0.2, 0.2, 2.0, 0.2, 0.5, 0.2,
                              pool->mutParams.maxSegmentos + 2, true);
        break;
      case 600:
        SetMutationParameters(pool, 0.5, 0.05, 1.0, 0.2, 0.5, 0.2,
                              pool->mutParams.maxSegmentos + 10, true);
        break;
      case 700:
        SetMutationParameters(pool, 0.5, 0.2, 1.0, 1.0, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 10, true);
        break;
      case 710:
        SetMutationParameters(pool, 0.5, 0.2, 1.0, 1.0, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 10, false);
        break;
      case 720:
        SetMutationParameters(pool, 0.5, 0.1, 1.0, 0.95, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 10, false);
        break;
      case 730:
        SetMutationParameters(pool, 0.5, 0.04, 1.0, 0.95, 0.0, 0.2,
                              pool->mutParams.maxSegmentos + 20, false);
        break;

      default:
        break;
    }
    ReproducePool(pool, sampleInfo);
  }
}
