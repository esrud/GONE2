#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "constants.hpp"
#include "population.hpp"

bool ReadTped(std::string fichTped, PopulationInfo* popInfo);
