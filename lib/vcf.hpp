#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "constants.hpp"
#include "population.hpp"

bool ReadVcf(std::string fichVcf, PopulationInfo* popInfo);
