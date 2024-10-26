#pragma once

#include <time.h>
#include <math.h>
#include <stdlib.h>

#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <cstring>
#include <limits>

#include "./sample.hpp"
#include "./pool.hpp"
#include "./console.hpp"
#include "./population.hpp"

void gone(AppParams* params, std::string fichero, int argc, char* argv[], PopulationInfo *popInfo, SampleInfo* sInfo);
