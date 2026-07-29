
#include "ncurses_tui.hpp"

#ifdef GONE_NCURSES_TUI

#include <atomic>
#include <algorithm>
#include <chrono>
#include <clocale>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <ncurses.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>

#include "progress.hpp"

namespace {

std::thread g_thread;
std::atomic<bool> g_running{false};
std::atomic<bool> g_stop{false};
std::atomic<bool> g_finished{false};
volatile std::sig_atomic_t g_sigintReceived = 0;
ProgressStatus* g_progress = nullptr;
std::mutex g_screenMu;
bool g_hasColors = false;
std::stringstream* g_cerrCapture = nullptr;
std::streambuf*    g_cerrOrigBuf = nullptr;

constexpr int kLeftMargin = 2;

enum ColorId {
  CP_TITLE      = 1,
  CP_BORDER     = 2,
  CP_DONE       = 3,
  CP_CURRENT    = 4,
  CP_PENDING    = 5,
  CP_BAR_FILL   = 6,
  CP_BAR_EMPTY  = 7,
  CP_AXIS       = 8,
  CP_CURVE      = 9,
  CP_REFERENCE  = 13,
  CP_STATUS     = 10,
  CP_HINT       = 11,
  CP_AXIS_LABEL = 12,
  CP_BAND_MIN   = 15,
  CP_BAND_MAX   = 16,
  CP_DENSITY    = 17,
};

// Initialise the ncurses colour pairs used by the TUI (no-op when the
// terminal has no colour support).
void InitColors() {
  if (!has_colors()) return;
  start_color();
  use_default_colors();
  init_pair(CP_TITLE,      COLOR_CYAN,    -1);
  init_pair(CP_BORDER,     COLOR_BLUE,    -1);
  init_pair(CP_DONE,       COLOR_GREEN,   -1);
  init_pair(CP_CURRENT,    COLOR_YELLOW,  -1);
  init_pair(CP_PENDING,    COLOR_WHITE,   -1);
  init_pair(CP_BAR_FILL,   COLOR_GREEN,   -1);
  init_pair(CP_BAR_EMPTY,  COLOR_WHITE,   -1);
  init_pair(CP_AXIS,       COLOR_BLUE,    -1);
  init_pair(CP_CURVE,      COLOR_CYAN,    -1);
  init_pair(CP_STATUS,     COLOR_WHITE,   -1);
  init_pair(CP_HINT,       COLOR_BLACK,   -1);
  init_pair(CP_AXIS_LABEL, COLOR_WHITE,   -1);
  init_pair(CP_REFERENCE,  COLOR_RED,     -1);
  init_pair(CP_BAND_MIN,   COLOR_GREEN,   -1);
  init_pair(CP_BAND_MAX,   COLOR_YELLOW,  -1);
  init_pair(CP_DENSITY,    COLOR_WHITE,   -1);
  g_hasColors = true;
}

struct Attr {
  int attr;
  Attr(int a) : attr(a) { attron(attr); }
  ~Attr() { attroff(attr); }
};

int Color(ColorId id) { return g_hasColors ? COLOR_PAIR(id) : 0; }

// Draw a rounded single-line box border.
// Params: top/left — top-left corner; height/width — box size in cells.
void DrawBox(int top, int left, int height, int width) {
  if (width < 2 || height < 2) return;
  Attr a(Color(CP_BORDER));
  mvaddstr(top, left, "╭");
  for (int x = 1; x < width - 1; ++x) mvaddstr(top, left + x, "─");
  mvaddstr(top, left + width - 1, "╮");
  for (int y = 1; y < height - 1; ++y) {
    mvaddstr(top + y, left, "│");
    mvaddstr(top + y, left + width - 1, "│");
  }
  mvaddstr(top + height - 1, left, "╰");
  for (int x = 1; x < width - 1; ++x) mvaddstr(top + height - 1, left + x, "─");
  mvaddstr(top + height - 1, left + width - 1, "╯");
}

// Draw a bold title centred (and clipped) within a width.
// Params: row — text row; width — available width; text — title string.
void DrawCenteredTitle(int row, int width, const char* text) {
  const int len = static_cast<int>(std::strlen(text));
  const int trimmed = std::min(len, width);
  const int col = std::max(0, (width - trimmed) / 2);
  Attr a(Color(CP_TITLE) | A_BOLD);
  mvaddnstr(row, col, text, trimmed);
}

static const char* kEighthBlocks[9] = {
  " ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█",
};

// Draw a progress bar in [col, col+width) using 1/8-block resolution.
// Params: row/col — bar start cell; width — bar width in cells; progress —
//   fill fraction in [0,1].
void DrawProgressBar(int row, int col, int width, double progress) {
  if (width < 1) return;
  const double clamped = std::max(0.0, std::min(1.0, progress));
  const int totalEighths = width * 8;
  const int filledEighths =
      std::max(0, std::min(totalEighths,
                           static_cast<int>(std::round(clamped * totalEighths))));
  const int fullCells = filledEighths / 8;
  const int rem = filledEighths % 8;

  {
    Attr a(Color(CP_BAR_FILL));
    for (int i = 0; i < fullCells; ++i) mvaddstr(row, col + i, "█");
    if (rem > 0 && fullCells < width) {
      mvaddstr(row, col + fullCells, kEighthBlocks[rem]);
    }
  }
  {
    Attr a(Color(CP_BAR_EMPTY) | A_DIM);
    const int emptyStart = fullCells + (rem > 0 ? 1 : 0);
    for (int i = emptyStart; i < width; ++i) mvaddstr(row, col + i, "·");
  }
}

// Draw a right-aligned percentage (e.g. " 12.3%").
// Params: row/col — start cell; pct — fraction in [0,1] shown as a percent.
void DrawPctRight(int row, int col, double pct) {
  char buf[16];
  std::snprintf(buf, sizeof(buf), "%5.1f%%", pct * 100.0);
  Attr a(Color(CP_PENDING));
  mvaddstr(row, col, buf);
}

// Draw the phase list — one row each with a status glyph, name, progress
// bar, and percentage.
// Params: top/left — top-left of the list; width — available width;
//   phases — the phases to render.
void DrawPhases(int top, int left, int width,
                const std::vector<ProgressPhaseView>& phases) {
  int row = top;
  const int barCol  = left + 32;
  const int pctCol  = std::max(barCol + 12, width - 9);
  const int barWidth = std::max(8, pctCol - barCol - 1);

  for (const ProgressPhaseView& phase : phases) {
    ColorId cid = CP_PENDING;
    const char* glyph = "○";
    int extraAttr = A_DIM;
    if (phase.state == ProgressPhaseState::Done) {
      cid = CP_DONE; glyph = "✔"; extraAttr = 0;
    } else if (phase.state == ProgressPhaseState::Current) {
      cid = CP_CURRENT; glyph = "▶"; extraAttr = A_BOLD;
    }
    {
      Attr a(Color(cid) | extraAttr);
      mvaddstr(row, left, glyph);
      mvaddch(row, left + 1, ' ');
      const int nameMax = barCol - left - 3;
      mvaddnstr(row, left + 2, phase.name.c_str(),
                std::min<int>(phase.name.size(), nameMax));
    }
    if (phase.state != ProgressPhaseState::Pending) {
      DrawProgressBar(row, barCol, barWidth, phase.progress);
      DrawPctRight(row, pctCol, phase.progress);
    }
    ++row;
  }
}

constexpr unsigned char kBrailleBit[4][2] = {
  {0x01, 0x08},
  {0x02, 0x10},
  {0x04, 0x20},
  {0x40, 0x80},
};

// Encode a braille bitmask into a 4-byte (UTF-8 + NUL) buffer.
// Params: mask — 8-bit braille dot bitmask; buf — receives the 3 UTF-8
//   bytes of the glyph plus a NUL terminator.
void EncodeBraille(unsigned char mask, char buf[4]) {
  const unsigned int cp = 0x2800u + mask;
  buf[0] = 0xE2;
  buf[1] = 0xA0 | ((cp >> 6) & 0x3F);
  buf[2] = 0x80 | (cp & 0x3F);
  buf[3] = 0;
}

// Draw the live Ne chart: bordered box, log10 Y axis with tick labels, the
// primary Ne curve in braille, the optional -f reference and -E band/
// density overlays, the gen-1 axis marker, and the pinned best-score box.
// Params: top/left — top-left corner; height/width — chart size in cells;
//   ne — the curve and overlays to draw; bestScore — value for the pinned
//   score box.
void DrawNeChart(int top, int left, int height, int width,
                 const NeEstimateView& ne, const BestScoreView& bestScore) {
  if (width < 32 || height < 8) return;

  DrawBox(top, left, height, width);
  {
    Attr a(Color(CP_TITLE) | A_BOLD);
    mvaddstr(top, left + 2, "─ Ne estimate (live) ");
  }
  if (!ne.unitLabel.empty()) {
    struct LegendSeg { std::string text; ColorId color; };
    std::vector<LegendSeg> segs;
    segs.push_back({ne.unitLabel, CP_CURVE});
    if (!ne.maxValues.empty() || !ne.minValues.empty()) {
      const std::string bandPrefix =
          (ne.unitLabel.rfind("N_T", 0) == 0) ? "N_T" : "Ne";
      if (!ne.maxValues.empty())
        segs.push_back({bandPrefix + "_max", CP_BAND_MAX});
      if (!ne.minValues.empty())
        segs.push_back({bandPrefix + "_min", CP_BAND_MIN});
    }
    int textLen = 0;
    for (const LegendSeg& s : segs) textLen += static_cast<int>(s.text.size());
    const int gaps = static_cast<int>(segs.size()) - 1;
    const int reserved = 2 + textLen + gaps + 2;
    const int col = std::max<int>(left + 24,
                                  left + width - 2 - reserved);
    int x = col;
    {
      Attr a(Color(CP_AXIS_LABEL) | A_DIM);
      mvaddstr(top, x, "─ ");
    }
    x += 2;
    for (std::size_t i = 0; i < segs.size(); ++i) {
      if (i > 0) {
        mvaddstr(top, x, " ");
        x += 1;
      }
      Attr a(Color(segs[i].color) | A_BOLD);
      mvaddstr(top, x, segs[i].text.c_str());
      x += static_cast<int>(segs[i].text.size());
    }
    {
      Attr a(Color(CP_AXIS_LABEL) | A_DIM);
      mvaddstr(top, x, " ─");
    }
  }

  const int yLabelWidth = 7;
  const int plotLeft = left + 2 + yLabelWidth + 1;
  const int plotRight = left + width - 3;
  const int plotCols = plotRight - plotLeft + 1;
  const int plotTop = top + 2;
  const int plotBottom = top + height - 4;
  const int plotRows = plotBottom - plotTop + 1;
  if (plotCols < 10 || plotRows < 4) return;

  constexpr int kBoxWidth  = 28;
  constexpr int kBoxHeight = 4;
  const bool showScoreBox = bestScore.hasValue &&
                            plotCols >= kBoxWidth + 24 &&
                            plotRows >= kBoxHeight + 4;
  const int boxRight  = left + width - 2;
  const int boxLeft   = boxRight - kBoxWidth + 1;
  const int boxTop    = top + 1;
  const int boxBottom = boxTop + kBoxHeight - 1;

  if (!ne.hasValues || ne.values.empty()) {
    Attr a(Color(CP_STATUS) | A_DIM);
    mvaddstr(top + height / 2, left + (width - 36) / 2,
             "waiting for the first GA round…");
    return;
  }

  constexpr int kMaxGen = 150;
  std::vector<double> v;
  v.reserve(std::min<std::size_t>(ne.values.size(), kMaxGen));
  for (std::size_t i = 0;
       i < ne.values.size() && i < static_cast<std::size_t>(kMaxGen); ++i) {
    v.push_back(std::max(1.0, ne.values[i]));
  }
  std::vector<double> vRef;
  vRef.reserve(std::min<std::size_t>(ne.referenceValues.size(), kMaxGen));
  for (std::size_t i = 0;
       i < ne.referenceValues.size() &&
       i < static_cast<std::size_t>(kMaxGen); ++i) {
    vRef.push_back(std::max(1.0, ne.referenceValues[i]));
  }
  std::vector<double> vMin, vMax;
  for (std::size_t i = 0;
       i < ne.minValues.size() && i < static_cast<std::size_t>(kMaxGen); ++i) {
    vMin.push_back(ne.minValues[i]);
  }
  for (std::size_t i = 0;
       i < ne.maxValues.size() && i < static_cast<std::size_t>(kMaxGen); ++i) {
    vMax.push_back(ne.maxValues[i]);
  }
  while (!vMin.empty() && vMin.back() <= 0.0) vMin.pop_back();
  while (!vMax.empty() && vMax.back() <= 0.0) vMax.pop_back();
  if (v.empty()) return;
  double mn = *std::min_element(v.begin(), v.end());
  double mx = *std::max_element(v.begin(), v.end());
  if (!vRef.empty()) {
    mn = std::min(mn, *std::min_element(vRef.begin(), vRef.end()));
    mx = std::max(mx, *std::max_element(vRef.begin(), vRef.end()));
  }
  for (const double x : vMin) if (x > 0.0) mn = std::min(mn, x);
  for (const double x : vMax) if (x > 0.0) mx = std::max(mx, x);
  int decadeMax = static_cast<int>(std::ceil(std::log10(mx) + 0.5));
  int decadeMin = decadeMax - 2;
  const int dataMinDecade = static_cast<int>(std::floor(std::log10(mn)));
  if (dataMinDecade < decadeMin) decadeMin = dataMinDecade;
  if (decadeMin < 0) decadeMin = 0;
  if (decadeMax > 7) decadeMax = 7;
  if (decadeMax - decadeMin < 2) {
    decadeMin = std::max(0, decadeMax - 2);
  }
  const double lmn = static_cast<double>(decadeMin);
  const double lmx = static_cast<double>(decadeMax);

  {
    Attr a(Color(CP_AXIS));
    for (int r = 0; r < plotRows; ++r) {
      mvaddstr(plotTop + r, plotLeft - 1, "│");
    }
  }
  {
    Attr a(Color(CP_AXIS_LABEL) | A_DIM);
    for (int d = decadeMin; d <= decadeMax; ++d) {
      const double frac = (lmx == lmn) ? 1.0 : (d - lmn) / (lmx - lmn);
      const int row =
          plotTop + static_cast<int>(std::round((1.0 - frac) * (plotRows - 1)));
      if (row < plotTop || row > plotBottom) continue;
      const double linv = std::pow(10.0, static_cast<double>(d));
      char buf[16];
      if (d >= 4)        std::snprintf(buf, sizeof(buf), "%5.0e", linv);
      else               std::snprintf(buf, sizeof(buf), "%5.0f", linv);
      mvaddstr(row, plotLeft - 1 - yLabelWidth + 2, buf);
      mvaddstr(row, plotLeft - 1, "┤");
    }
  }
  {
    const double gen1 = v[0];
    const double frac =
        (lmx == lmn) ? 1.0 : (std::log10(gen1) - lmn) / (lmx - lmn);
    const int row =
        plotTop + static_cast<int>(std::round((1.0 - frac) * (plotRows - 1)));
    if (row >= plotTop && row <= plotBottom) {
      char buf[16];
      if (gen1 >= 1e4) std::snprintf(buf, sizeof(buf), "%5.0e", gen1);
      else             std::snprintf(buf, sizeof(buf), "%5.0f", gen1);
      {
        Attr a(Color(CP_CURVE) | A_BOLD);
        mvaddstr(row, plotLeft - 1 - yLabelWidth + 2, buf);
        mvaddstr(row, plotLeft - 1, "►");
      }
    }
  }

  std::vector<unsigned char> cells(plotCols * plotRows, 0);
  const int pixW = plotCols * 2;
  const int pixH = plotRows * 4;
  const int N = static_cast<int>(v.size());

  auto sampleYAt = [&](double xPix) -> int {
    if (N <= 1) {
      const double frac = (std::log10(v[0]) - lmn) /
                          std::max(1e-9, (lmx - lmn));
      return std::max(0, std::min(pixH - 1,
                                  (int)std::round((1.0 - frac) * (pixH - 1))));
    }
    double t = xPix * (N - 1.0) / std::max(1, pixW - 1);
    if (t < 0) t = 0;
    if (t > N - 1) t = N - 1;
    const int i = std::min(N - 2, std::max(0, (int)t));
    const double frac = t - i;
    const double lv = std::log10(v[i]) +
                      frac * (std::log10(v[i + 1]) - std::log10(v[i]));
    const double yFrac = (lv - lmn) / std::max(1e-9, (lmx - lmn));
    int y = (int)std::round((1.0 - yFrac) * (pixH - 1));
    if (y < 0) y = 0;
    if (y >= pixH) y = pixH - 1;
    return y;
  };

  const int boxCol0 = showScoreBox ? boxLeft  - plotLeft : 0;
  const int boxCol1 = showScoreBox ? boxRight - plotLeft : -1;
  const int boxRow0 = showScoreBox ? boxTop    - plotTop : 0;
  const int boxRow1 = showScoreBox ? boxBottom - plotTop : -1;
  auto inScoreBox = [&](int textRow, int textCol) {
    return showScoreBox &&
           textRow >= boxRow0 && textRow <= boxRow1 &&
           textCol >= boxCol0 && textCol <= boxCol1;
  };

  auto setPix = [&](int x, int y) {
    if (x < 0 || x >= pixW || y < 0 || y >= pixH) return;
    const int textRow = y / 4;
    const int textCol = x / 2;
    if (inScoreBox(textRow, textCol)) return;
    cells[textRow * plotCols + textCol] |= kBrailleBit[y % 4][x % 2];
  };

  for (int xp = 0; xp < pixW; ++xp) {
    const double xLeft  = std::max(0.0, xp - 0.5);
    const double xRight = std::min((double)(pixW - 1), xp + 0.5);
    const int yL = sampleYAt(xLeft);
    const int yC = sampleYAt(xp);
    const int yR = sampleYAt(xRight);
    const int yMin = std::min({yL, yC, yR});
    const int yMax = std::max({yL, yC, yR});
    for (int y = yMin; y <= yMax; ++y) setPix(xp, y);
  }

  if (!vRef.empty()) {
    std::vector<unsigned char> refCells(plotCols * plotRows, 0);
    const int N_ref = static_cast<int>(vRef.size());
    auto sampleRefYAt = [&](double xPix) -> int {
      if (N_ref <= 1) {
        const double frac = (std::log10(vRef[0]) - lmn) /
                            std::max(1e-9, (lmx - lmn));
        return std::max(0, std::min(pixH - 1,
                                    (int)std::round((1.0 - frac) *
                                                    (pixH - 1))));
      }
      double t = xPix * (N_ref - 1.0) / std::max(1, pixW - 1);
      if (t < 0) t = 0;
      if (t > N_ref - 1) t = N_ref - 1;
      const int i = std::min(N_ref - 2, std::max(0, (int)t));
      const double frac = t - i;
      const double lv = std::log10(vRef[i]) +
                        frac * (std::log10(vRef[i + 1]) - std::log10(vRef[i]));
      const double yFrac = (lv - lmn) / std::max(1e-9, (lmx - lmn));
      int y = (int)std::round((1.0 - yFrac) * (pixH - 1));
      if (y < 0) y = 0;
      if (y >= pixH) y = pixH - 1;
      return y;
    };
    auto setRefPix = [&](int x, int y) {
      if (x < 0 || x >= pixW || y < 0 || y >= pixH) return;
      const int textRow = y / 4;
      const int textCol = x / 2;
      if (inScoreBox(textRow, textCol)) return;
      refCells[textRow * plotCols + textCol] |= kBrailleBit[y % 4][x % 2];
    };
    for (int xp = 0; xp < pixW; ++xp) {
      const double xLeft  = std::max(0.0, xp - 0.5);
      const double xRight = std::min((double)(pixW - 1), xp + 0.5);
      const int yL = sampleRefYAt(xLeft);
      const int yC = sampleRefYAt(xp);
      const int yR = sampleRefYAt(xRight);
      const int yMin = std::min({yL, yC, yR});
      const int yMax = std::max({yL, yC, yR});
      for (int y = yMin; y <= yMax; ++y) setRefPix(xp, y);
    }
    Attr a(Color(CP_REFERENCE) | A_BOLD);
    char buf[4];
    for (int r = 0; r < plotRows; ++r) {
      for (int c = 0; c < plotCols; ++c) {
        const unsigned char mask = refCells[r * plotCols + c];
        if (mask == 0) continue;
        EncodeBraille(mask, buf);
        mvaddstr(plotTop + r, plotLeft + c, buf);
      }
    }
  }

  {
    auto pixYOf = [&](double val) -> int {
      const double lv = std::log10(std::max(1.0, val));
      const double yFrac = (lv - lmn) / std::max(1e-9, (lmx - lmn));
      int y = (int)std::round((1.0 - yFrac) * (pixH - 1));
      if (y < 0) y = 0;
      if (y >= pixH) y = pixH - 1;
      return y;
    };
    if (!ne.densitySamples.empty()) {
      std::vector<unsigned char> dCells(plotCols * plotRows, 0);
      const int M = std::min<int>(static_cast<int>(ne.densitySamples.size()),
                                  kMaxGen);
      const int denom = std::max(1, M - 1);
      bool any = false;
      for (int g = 0; g < M; ++g) {
        const int xp =
            static_cast<int>(std::round(g * (double)(pixW - 1) / denom));
        if (xp < 0 || xp >= pixW) continue;
        for (const double s : ne.densitySamples[g]) {
          if (s <= 0.0) continue;
          const int y = pixYOf(s);
          const int textRow = y / 4, textCol = xp / 2;
          if (inScoreBox(textRow, textCol)) continue;
          dCells[textRow * plotCols + textCol] |= kBrailleBit[y % 4][xp % 2];
          any = true;
        }
      }
      if (any) {
        Attr a(Color(CP_DENSITY) | A_DIM);
        char buf[4];
        for (int r = 0; r < plotRows; ++r)
          for (int c = 0; c < plotCols; ++c) {
            const unsigned char mask = dCells[r * plotCols + c];
            if (mask == 0) continue;
            EncodeBraille(mask, buf);
            mvaddstr(plotTop + r, plotLeft + c, buf);
          }
      }
    }
    auto drawBandCurve = [&](const std::vector<double>& vv, int attr) {
      if (vv.empty()) return;
      const int Nc = static_cast<int>(vv.size());
      std::vector<unsigned char> bc(plotCols * plotRows, 0);
      auto yAt = [&](double xPix) -> int {
        if (Nc <= 1) return pixYOf(vv[0]);
        double t = xPix * (Nc - 1.0) / std::max(1, pixW - 1);
        if (t < 0) t = 0;
        if (t > Nc - 1) t = Nc - 1;
        const int i = std::min(Nc - 2, std::max(0, (int)t));
        const double frac = t - i;
        const double lv = std::log10(std::max(1.0, vv[i])) +
                          frac * (std::log10(std::max(1.0, vv[i + 1])) -
                                  std::log10(std::max(1.0, vv[i])));
        const double yFrac = (lv - lmn) / std::max(1e-9, (lmx - lmn));
        int y = (int)std::round((1.0 - yFrac) * (pixH - 1));
        if (y < 0) y = 0;
        if (y >= pixH) y = pixH - 1;
        return y;
      };
      auto setp = [&](int x, int y) {
        if (x < 0 || x >= pixW || y < 0 || y >= pixH) return;
        const int textRow = y / 4, textCol = x / 2;
        if (inScoreBox(textRow, textCol)) return;
        bc[textRow * plotCols + textCol] |= kBrailleBit[y % 4][x % 2];
      };
      for (int xp = 0; xp < pixW; ++xp) {
        const double xLeft  = std::max(0.0, xp - 0.5);
        const double xRight = std::min((double)(pixW - 1), xp + 0.5);
        const int yMin = std::min({yAt(xLeft), yAt(xp), yAt(xRight)});
        const int yMax = std::max({yAt(xLeft), yAt(xp), yAt(xRight)});
        for (int y = yMin; y <= yMax; ++y) setp(xp, y);
      }
      Attr a(attr);
      char buf[4];
      for (int r = 0; r < plotRows; ++r)
        for (int c = 0; c < plotCols; ++c) {
          const unsigned char mask = bc[r * plotCols + c];
          if (mask == 0) continue;
          EncodeBraille(mask, buf);
          mvaddstr(plotTop + r, plotLeft + c, buf);
        }
    };
    drawBandCurve(vMax, Color(CP_BAND_MAX) | A_BOLD);
    drawBandCurve(vMin, Color(CP_BAND_MIN) | A_BOLD);
  }

  {
    Attr a(Color(CP_CURVE) | A_BOLD);
    char buf[4];
    for (int r = 0; r < plotRows; ++r) {
      for (int c = 0; c < plotCols; ++c) {
        const unsigned char mask = cells[r * plotCols + c];
        if (mask == 0) continue;
        EncodeBraille(mask, buf);
        mvaddstr(plotTop + r, plotLeft + c, buf);
      }
    }
  }

  if (showScoreBox) {
    {
      Attr a(Color(CP_BORDER));
      for (int r = boxTop; r <= boxBottom; ++r) {
        for (int c = boxLeft; c <= boxRight; ++c) {
          mvaddstr(r, c, " ");
        }
      }
    }
    DrawBox(boxTop, boxLeft, kBoxHeight, kBoxWidth);
    const std::string boxTitle =
        bestScore.title.empty() ? "Best score" : bestScore.title;
    {
      Attr a(Color(CP_TITLE) | A_BOLD);
      const std::string titleBar = "─ " + boxTitle + " ";
      mvaddstr(boxTop, boxLeft + 2, titleBar.c_str());
    }
    char vbuf[40];
    const double absScore = std::fabs(bestScore.score);
    if (absScore != 0 && (absScore < 0.01 || absScore >= 1e5)) {
      std::snprintf(vbuf, sizeof(vbuf), "%.4g", bestScore.score);
    } else {
      std::snprintf(vbuf, sizeof(vbuf), "%.4f", bestScore.score);
    }
    {
      Attr a(Color(CP_CURVE) | A_BOLD);
      mvaddstr(boxTop + 1, boxLeft + 2, vbuf);
    }
    char rbuf[40];
    if (bestScore.subtitle.empty()) {
      std::snprintf(rbuf, sizeof(rbuf), "round %d / %d",
                    bestScore.completedRounds, bestScore.totalRounds);
    } else {
      std::snprintf(rbuf, sizeof(rbuf), "%s", bestScore.subtitle.c_str());
    }
    {
      Attr a(Color(CP_STATUS) | A_DIM);
      mvaddstr(boxTop + 2, boxLeft + 2, rbuf);
    }
  }

  const int xAxisRow = plotBottom + 1;
  {
    Attr a(Color(CP_AXIS));
    mvaddstr(xAxisRow, plotLeft - 1, "└");
    for (int c = 0; c < plotCols; ++c) {
      mvaddstr(xAxisRow, plotLeft + c, "─");
    }
  }
  {
    Attr a(Color(CP_AXIS_LABEL) | A_DIM);
    char lhs[16], rhs[16], mid[16];
    std::snprintf(lhs, sizeof(lhs), "gen 1");
    std::snprintf(rhs, sizeof(rhs), "gen %d", N);
    std::snprintf(mid, sizeof(mid), "gen %d", N / 2);
    mvaddstr(xAxisRow + 1, plotLeft, lhs);
    const int midCol = plotLeft + plotCols / 2 - (int)std::strlen(mid) / 2;
    mvaddstr(xAxisRow + 1, midCol, mid);
    const int rcol = std::max<int>(plotLeft, plotLeft + plotCols
                                              - (int)std::strlen(rhs));
    mvaddstr(xAxisRow + 1, rcol, rhs);
  }
}

// Render one full TUI frame (title, phase list, status/warning lines, and
// the Ne chart) from a progress snapshot.
// Params: s — snapshot of the progress state to draw.
void Render(const ProgressSnapshot& s) {
  std::lock_guard<std::mutex> lk(g_screenMu);
  erase();
  int rows = 0, cols = 0;
  getmaxyx(stdscr, rows, cols);
  if (rows < 12 || cols < 60) {
    Attr a(Color(CP_STATUS));
    mvaddstr(0, 0, "Terminal too small — resize to at least 60x12.");
    refresh();
    return;
  }

  DrawCenteredTitle(0, cols, "GONE2 — Genetic Optimization for Ne Estimation");

  {
    Attr a(Color(CP_STATUS));
    mvaddstr(2, kLeftMargin, "Overall");
  }
  const int gbarCol = kLeftMargin + 9;
  const int gbarRight = cols - 10;
  const int gbarWidth = std::max(10, gbarRight - gbarCol);
  DrawProgressBar(2, gbarCol, gbarWidth, s.globalProgress);
  DrawPctRight(2, cols - 8, s.globalProgress);

  const int phaseTop = 4;
  DrawPhases(phaseTop, kLeftMargin, cols, s.phases);

  int warnRow = phaseTop + (int)s.phases.size() + 1;
  int statusRow = warnRow;
  if (!s.warning.empty() && warnRow < rows - 2) {
    const int textCol = kLeftMargin + 2;
    const int textWidth = std::max(1, cols - textCol - 2);
    std::vector<std::string> lines;
    {
      std::string current;
      std::size_t i = 0;
      const std::string& msg = s.warning;
      while (i < msg.size()) {
        while (i < msg.size() && msg[i] == ' ') ++i;
        if (i >= msg.size()) break;
        std::size_t j = i;
        while (j < msg.size() && msg[j] != ' ') ++j;
        std::string word = msg.substr(i, j - i);
        i = j;
        if (current.empty()) {
          current = word;
        } else if ((int)(current.size() + 1 + word.size()) <= textWidth) {
          current += ' ';
          current += word;
        } else {
          lines.push_back(current);
          current = word;
        }
        while ((int)current.size() > textWidth) {
          lines.push_back(current.substr(0, textWidth));
          current = current.substr(textWidth);
        }
      }
      if (!current.empty()) lines.push_back(current);
    }
    Attr a(Color(CP_CURRENT) | A_BOLD);
    const int maxWarnRows = std::max(0, (rows - 2) - warnRow);
    const int drawRows = std::min<int>(lines.size(), maxWarnRows);
    for (int k = 0; k < drawRows; ++k) {
      if (k == 0) mvaddstr(warnRow, kLeftMargin, "⚠ ");
      mvaddnstr(warnRow + k, textCol, lines[k].c_str(),
                std::min<int>(lines[k].size(), textWidth));
    }
    statusRow = warnRow + drawRows;
  }

  if (statusRow < rows - 2 && !s.statusDetail.empty()) {
    Attr a(Color(CP_STATUS));
    mvaddstr(statusRow, kLeftMargin, "› ");
    mvaddnstr(statusRow, kLeftMargin + 2, s.statusDetail.c_str(),
              std::min<int>(s.statusDetail.size(), cols - kLeftMargin - 4));
  }

  const int chartTop = statusRow + 2;
  const int chartHeight = rows - chartTop - 2;
  if (chartHeight >= 8) {
    DrawNeChart(chartTop, kLeftMargin, chartHeight, cols - 2 * kLeftMargin,
                s.ne, s.bestScore);
  }

  {
    Attr a(Color(CP_HINT) | A_DIM);
    const char* hint = g_finished.load(std::memory_order_relaxed)
                           ? "press Ctrl-C to exit"
                           : "press Ctrl-C to abort";
    mvaddstr(rows - 1, std::max<int>(0, cols - (int)std::strlen(hint) - 1),
             hint);
  }
  refresh();
}

// Background render thread: repaint the screen ~10×/sec from the latest
// progress snapshot until asked to stop, draining keyboard input so resize
// events don't pile up.
void RendererLoop() {
  sigset_t set;
  sigemptyset(&set);
  sigaddset(&set, SIGINT);
  pthread_sigmask(SIG_BLOCK, &set, nullptr);

  while (!g_stop.load(std::memory_order_relaxed)) {
    ProgressSnapshot snap = g_progress->Snapshot();
    Render(snap);
    int ch;
    while ((ch = getch()) != ERR) { (void)ch; }
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  if (g_progress) Render(g_progress->Snapshot());
}

// SIGINT handler: record that Ctrl-C was received so WaitForSigint unblocks.
// Params: (signal number, unused).
extern "C" void HandleSigint(int) {
  g_sigintReceived = 1;
}

}

// atexit hook that restores the terminal on any exit path bypassing normal
// cleanup; idempotent, delegates to StopNcursesTui.
extern "C" void AtExitCleanupTui() {
  StopNcursesTui();
}

// Start the ncurses TUI: enter screen mode, set up colours, capture cerr,
// register the atexit cleanup, and launch the background render thread. A
// second call while already running is a no-op.
// Params: progress — the progress state the render thread reads from.
void StartNcursesTui(ProgressStatus* progress) {
  if (g_running.exchange(true)) return;
  g_progress = progress;
  g_stop.store(false);

  if (g_cerrOrigBuf == nullptr) {
    g_cerrCapture = new std::stringstream;
    g_cerrOrigBuf = std::cerr.rdbuf();
    std::cerr.rdbuf(g_cerrCapture->rdbuf());
  }

  std::setlocale(LC_ALL, "");
  initscr();
  cbreak();
  noecho();
  curs_set(0);
  nodelay(stdscr, TRUE);
  keypad(stdscr, TRUE);
  InitColors();

  static bool atexitRegistered = false;
  if (!atexitRegistered) {
    std::atexit(AtExitCleanupTui);
    atexitRegistered = true;
  }

  g_thread = std::thread(RendererLoop);
}

// Block the main thread until the user presses Ctrl-C, used to keep the
// finished chart on screen. No-op if the TUI isn't running.
void WaitForSigint() {
  if (!g_running.load()) return;
  g_finished.store(true, std::memory_order_relaxed);
  g_sigintReceived = 0;

  struct sigaction sa{};
  sa.sa_handler = HandleSigint;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;
  struct sigaction old{};
  sigaction(SIGINT, &sa, &old);

  sigset_t set;
  sigemptyset(&set);
  sigaddset(&set, SIGINT);
  pthread_sigmask(SIG_UNBLOCK, &set, nullptr);

  while (g_sigintReceived == 0) {
    pause();
  }

  sigaction(SIGINT, &old, nullptr);
}

// Stop the render thread, leave ncurses mode, restore the terminal and the
// captured cerr, and flush any buffered stderr. Idempotent.
void StopNcursesTui() {
  if (!g_running.exchange(false)) return;
  g_stop.store(true);
  if (g_thread.joinable()) g_thread.join();

  {
    std::lock_guard<std::mutex> lk(g_screenMu);
    endwin();
  }
  g_progress = nullptr;

  if (g_cerrOrigBuf != nullptr) {
    std::cerr.rdbuf(g_cerrOrigBuf);
    g_cerrOrigBuf = nullptr;
    if (g_cerrCapture != nullptr) {
      const std::string captured = g_cerrCapture->str();
      if (!captured.empty()) {
        std::cerr << captured;
        std::cerr.flush();
      }
      delete g_cerrCapture;
      g_cerrCapture = nullptr;
    }
  }
}

#endif
