#include "ne_ref.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// Strip leading whitespace and trailing CR/whitespace from a string.
// Params: s — string trimmed in place.
void Trim(std::string& s) {
  while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(0, 1);
  while (!s.empty() && (s.back() == ' ' || s.back() == '\t' ||
                        s.back() == '\r' || s.back() == '\n')) s.pop_back();
}

// Test whether a string is a non-empty run of decimal digits.
// Params: s — string to test.
// Returns: true if s is non-empty and every char is 0-9, else false.
bool IsAllDigits(const std::string& s) {
  if (s.empty()) return false;
  for (char c : s) if (c < '0' || c > '9') return false;
  return true;
}

}

// Load a reference Ne(t) history from a file for the -f chart overlay.
// Accepts two formats: a GONE simu-params file (single-value parameter
// lines followed by `Ne nGen` demographic blocks of two integers each,
// oldest→most recent, expanded into a per-generation series in reverse), or
// a plain two-column `Gen Ne` table (sorted by generation, with gaps
// step-filled from the previous value). The simu-params format is detected
// by the presence of single-value numeric parameter lines, not by any
// comment in the file.
// Params: path — path to the reference file.
// Returns: a per-generation Ne vector (index 0 = most recent generation),
//   or an empty vector if the file cannot be opened or parsed.
std::vector<double> LoadReferenceNe(const std::string& path) {
  std::ifstream in(path);
  if (!in.good()) {
    std::cerr << "Could not open Ne reference file: " << path << "\n";
    return {};
  }

  std::vector<std::pair<long long, long long>> blocks;
  std::vector<std::pair<long long, double>> perGen;
  // A simu-params file carries single-value numeric parameter lines
  // (mutation rate, chromosome count, …) ahead of its demographic blocks;
  // a plain Gen/Ne table has none. Use that to tell the formats apart so
  // the two-integer demographic blocks are only read as `Ne nGen` pairs in
  // a simu-params file.
  bool sawParamLine = false;
  std::string line;
  while (std::getline(in, line)) {
    Trim(line);
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    std::istringstream iss(line);
    std::string t1, t2;
    if (!(iss >> t1)) continue;
    if (!std::isdigit(static_cast<unsigned char>(t1[0])) &&
        t1[0] != '-' && t1[0] != '+') {
      continue;
    }
    if (!(iss >> t2)) {
      sawParamLine = true;
      continue;
    }
    if (sawParamLine && IsAllDigits(t1) && IsAllDigits(t2)) {
      blocks.emplace_back(std::stoll(t1), std::stoll(t2));
    } else {
      try {
        const long long gen = std::stoll(t1);
        const double ne     = std::stod(t2);
        if (gen >= 0 && ne > 0) perGen.emplace_back(gen, ne);
      } catch (...) {
      }
    }
  }

  std::vector<double> out;
  if (!blocks.empty()) {
    for (auto it = blocks.rbegin(); it != blocks.rend(); ++it) {
      const double ne = static_cast<double>(it->first);
      const long long n_gen = it->second;
      for (long long i = 0; i < n_gen; ++i) out.push_back(ne);
    }
    return out;
  }
  if (!perGen.empty()) {
    std::sort(perGen.begin(), perGen.end());
    long long lastGen = -1;
    double lastNe = perGen[0].second;
    for (const auto& [gen, ne] : perGen) {
      while (lastGen + 1 < gen) {
        ++lastGen;
        out.push_back(lastNe);
      }
      lastGen = gen;
      lastNe = ne;
      out.push_back(ne);
    }
    return out;
  }
  std::cerr << "Could not parse Ne reference from: " << path
            << " (expected simu params or `Gen Ne` table)\n";
  return {};
}
