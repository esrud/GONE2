#pragma once

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

struct NeBand {
  std::vector<double> min;
  std::vector<double> max;
};

// `samples[g]` holds the per-round Ne values at generation g. Returns
// the per-generation min and max over those rounds; generations with
// no samples yield 0.0 (the renderer skips non-positive entries).
// Params: samples — per-generation vectors of each round's Ne value.
// Returns: a NeBand with the per-generation min and max curves.
inline NeBand ComputeNeBand(const std::vector<std::vector<double>>& samples) {
  NeBand band;
  band.min.assign(samples.size(), 0.0);
  band.max.assign(samples.size(), 0.0);
  for (std::size_t g = 0; g < samples.size(); ++g) {
    if (samples[g].empty()) continue;
    double mn = samples[g][0];
    double mx = samples[g][0];
    for (const double v : samples[g]) {
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }
    band.min[g] = mn;
    band.max[g] = mx;
  }
  return band;
}

// Long-format TSV of the per-round Ne behind the best-fit curve, for
// import into a spreadsheet: a "Generation\tNe" header followed by one
// row per sample (generation is 1-indexed, gen 1 = most recent), capped
// at generation 150. Generations with no samples are skipped.
// Params: samples — per-generation vectors of each round's Ne value.
// Returns: a "Generation\tNe" TSV string with one row per sample.
inline std::string FormatDensityTsv(
    const std::vector<std::vector<double>>& samples) {
  std::ostringstream out;
  out << "Generation\tNe\n";
  out << std::defaultfloat << std::setprecision(8);
  const std::size_t maxGen = std::min<std::size_t>(samples.size(), 150);
  for (std::size_t g = 0; g < maxGen; ++g) {
    for (const double v : samples[g]) {
      out << (g + 1) << "\t" << v << "\n";
    }
  }
  return out.str();
}
