#pragma once

enum class GASmoothKind {
  Quadratic = 0,
  L1        = 1,
  Truncated = 2,
};

struct GAConfig {
  double smoothLambda;
  GASmoothKind smoothKind;
  double smoothCutoff;
  double kickProb;
  double kickMag;
  double scoutFrac;
  const char* label;
};

// Build the default config from the compile-time -DGA_* flags that
// the rest of the codebase was built with. Used both as the
// configuration in non-combo builds and as one of the three combo
// presets.
// Returns: a GAConfig with the compile-time default smoothness/kick values.
inline GAConfig MakeDefaultGAConfig() {
  GAConfig c{};
#ifdef GA_SMOOTH_LAMBDA
  c.smoothLambda = static_cast<double>(GA_SMOOTH_LAMBDA);
#else
  c.smoothLambda = 0.0001;
#endif
#if defined(GA_SMOOTH_L1)
  c.smoothKind = GASmoothKind::L1;
#elif defined(GA_SMOOTH_CUTOFF)
  c.smoothKind = GASmoothKind::Truncated;
#else
  c.smoothKind = GASmoothKind::Quadratic;
#endif
#ifdef GA_SMOOTH_CUTOFF
  c.smoothCutoff = static_cast<double>(GA_SMOOTH_CUTOFF);
#else
  c.smoothCutoff = 0.5;
#endif
#ifdef GA_KICK_PROB
  c.kickProb = static_cast<double>(GA_KICK_PROB);
#else
  c.kickProb = 0.05;
#endif
#ifdef GA_KICK_MAG
  c.kickMag = static_cast<double>(GA_KICK_MAG);
#else
  c.kickMag = 10.0;
#endif
#ifdef GA_SCOUT_FRAC
  c.scoutFrac = static_cast<double>(GA_SCOUT_FRAC);
#else
  c.scoutFrac = 0.4;
#endif
  c.label = "";
  return c;
}

// The three presets used by the combo target. Order is the order
// of the sweep loop in libgone.cpp.
// Returns: the trunc05_kick preset (truncated smoothness + strong kicks).
inline GAConfig MakeComboTruncKickConfig() {
  GAConfig c = MakeDefaultGAConfig();
  c.smoothLambda = 0.0001;
  c.smoothKind   = GASmoothKind::Truncated;
  c.smoothCutoff = 0.5;
  c.kickProb     = 0.15;
  c.kickMag      = 30.0;
  c.scoutFrac    = 0.7;
  c.label        = "trunc05_kick";
  return c;
}

// Returns: the L2 preset (quadratic smoothness, gentle kicks).
inline GAConfig MakeComboL2Config() {
  GAConfig c = MakeDefaultGAConfig();
  c.smoothLambda = 0.0001;
  c.smoothKind   = GASmoothKind::Quadratic;
  c.smoothCutoff = 0.5;
  c.kickProb     = 0.05;
  c.kickMag      = 10.0;
  c.scoutFrac    = 0.4;
  c.label        = "L2";
  return c;
}

// Returns: the L1_kick preset (L1 smoothness + strong kicks).
inline GAConfig MakeComboL1KickConfig() {
  GAConfig c = MakeDefaultGAConfig();
  c.smoothLambda = 0.0001;
  c.smoothKind   = GASmoothKind::L1;
  c.smoothCutoff = 0.5;
  c.kickProb     = 0.15;
  c.kickMag      = 30.0;
  c.scoutFrac    = 0.7;
  c.label        = "L1_kick";
  return c;
}
