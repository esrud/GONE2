#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include "./constants.hpp"
#include "./population.hpp"
bool ReadFile(std::string fichPed, std::string fichMap,
              PopulationInfo* popInfo);
bool ReadMap(std::string fichMap, PopulationInfo* popInfo);
