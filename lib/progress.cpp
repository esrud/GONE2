#include "./progress.hpp"

#include <algorithm>

namespace {

// Clamp a value to the [0, 1] range.
// Params: value — value to clamp.
// Returns: value bounded to [0, 1].
double Clamp01(double value) {
  return std::max(0.0, std::min(1.0, value));
}

// Coerce a non-positive value to 1, for use as a safe denominator.
// Params: value — candidate value.
// Returns: value if > 0, otherwise 1.
double PositiveOrOne(double value) {
  return value > 0.0 ? value : 1.0;
}

}

// Replace the phase list with one named row per entry and refresh states.
// Params: phaseNames — ordered phase labels for the progress display.
void ProgressStatus::ConfigurePhases(
    const std::vector<std::string>& phaseNames) {
  std::lock_guard<std::mutex> lk(mu);
  phases.clear();
  phases.reserve(phaseNames.size());
  for (const std::string& name : phaseNames) {
    ProgressPhaseView phase;
    phase.name = name;
    phases.push_back(phase);
  }
  RefreshPhaseStates();
}

// Set the total number of phases and the legacy progress-file path.
// Params: nTasks — number of phases (coerced to ≥1); fname — path of the
//   text progress file to write (nullptr/empty disables it).
void ProgressStatus::InitTotalTasks(float nTasks, const char* fname) {
  std::lock_guard<std::mutex> lk(mu);
  totalTasks = PositiveOrOne(static_cast<double>(nTasks));
  progressFname = fname == nullptr ? "" : fname;
  EnsurePhaseList();
  RefreshPhaseStates();
}

// Switch to phase cTask, resetting the per-phase sub-window and floor.
// Params: cTask — phase index to activate (ignored if > totalTasks);
//   tName — label for the new phase.
void ProgressStatus::SetCurrentTask(float cTask, const char* tName) {
  std::unique_lock<std::mutex> lk(mu);
  if (cTask > totalTasks) {
    lk.unlock();
    std::cerr << "ERROR: Setting progress larger than max value\n";
    return;
  }
  currentTask = std::max(0.0, static_cast<double>(cTask));
  taskName = tName == nullptr ? "" : tName;
  taskWindowBase = 0.0;
  taskWindowSpan = 1.0;
  taskProgressFloor = 0.0;
  taskProgress = 0.0;
  currentStep = 0.0;
  globalProgress = Clamp01(currentTask / totalTasks);
  EnsurePhaseList();
  RefreshPhaseStates();
}

// Begin the current phase's step counter, starting the bar at the active
// sub-window base and never moving it backwards.
// Params: nSteps — number of steps in the phase (coerced to ≥1).
void ProgressStatus::InitCurrentTask(float nSteps) {
  std::lock_guard<std::mutex> lk(mu);
  totalSteps = PositiveOrOne(static_cast<double>(nSteps));
  currentStep = 0.0;
  taskProgress = std::max(Clamp01(taskWindowBase), taskProgressFloor);
  taskProgressFloor = taskProgress;
  globalProgress = Clamp01((currentTask + taskProgress) / totalTasks);
  RefreshPhaseStates();
}

// Confine subsequent SetTaskProgress updates to the sub-range
// [base, base+span] of the current phase, used to give each -x -T
// resampling round its own slice of the inference phase.
// Params: base — slice start in [0,1]; span — slice width (≥0).
void ProgressStatus::SetTaskWindow(double base, double span) {
  std::lock_guard<std::mutex> lk(mu);
  taskWindowBase = Clamp01(base);
  taskWindowSpan = std::max(0.0, span);
  taskProgress = std::max(taskWindowBase, taskProgressFloor);
  taskProgressFloor = taskProgress;
  globalProgress = Clamp01((currentTask + taskProgress) / totalTasks);
  RefreshPhaseStates();
}

// Advance the current phase's progress to step cStep (mapped through the
// active sub-window, never decreasing), then write the progress file.
// Params: cStep — current step within the phase (ignored if > totalSteps).
void ProgressStatus::SetTaskProgress(float cStep) {
  std::unique_lock<std::mutex> lk(mu);
  if (cStep > totalSteps) {
    lk.unlock();
    std::cerr << "ERROR: Setting progress bigger than max value\n";
    return;
  }
  currentStep = std::max(0.0, static_cast<double>(cStep));
  taskProgress = std::max(
      Clamp01(taskWindowBase + taskWindowSpan * (currentStep / totalSteps)),
      taskProgressFloor);
  taskProgressFloor = taskProgress;
  globalProgress = Clamp01((currentTask + taskProgress) / totalTasks);
  RefreshPhaseStates();
  lk.unlock();
  SaveProgress();
}

// Set the status-detail line shown under the phase list.
// Params: detail — status text.
void ProgressStatus::SetStatusDetail(const std::string& detail) {
  std::lock_guard<std::mutex> lk(mu);
  statusDetail = detail;
}

// Set the persistent warning line shown on the chart (empty clears it).
// Params: message — warning text.
void ProgressStatus::SetWarning(const std::string& message) {
  std::lock_guard<std::mutex> lk(mu);
  warning = message;
}

// Replace the primary Ne curve drawn on the chart.
// Params: values — per-generation Ne; unitLabel — legend label (e.g. "N_T").
void ProgressStatus::SetNeSnapshot(const std::vector<double>& values,
                                   const std::string& unitLabel) {
  std::lock_guard<std::mutex> lk(mu);
  ne.hasValues = !values.empty();
  ne.values = values;
  ne.unitLabel = unitLabel;
}

// Set the -f ground-truth Ne reference overlay (drawn in red).
// Params: values — per-generation reference Ne.
void ProgressStatus::SetReferenceNe(const std::vector<double>& values) {
  std::lock_guard<std::mutex> lk(mu);
  ne.referenceValues = values;
}

// Set the -E min/max bounding curves across the GA rounds.
// Params: minValues — per-gen lower bound; maxValues — per-gen upper bound.
void ProgressStatus::SetNeBand(const std::vector<double>& minValues,
                               const std::vector<double>& maxValues) {
  std::lock_guard<std::mutex> lk(mu);
  ne.minValues = minValues;
  ne.maxValues = maxValues;
}

// Set the -E per-round Ne density cloud.
// Params: samples — per-generation vectors of each round's Ne value.
void ProgressStatus::SetNeDensity(
    const std::vector<std::vector<double>>& samples) {
  std::lock_guard<std::mutex> lk(mu);
  ne.densitySamples = samples;
}

// Update the pinned best-score box, keeping the lowest score seen so the
// displayed value never moves backwards.
// Params: score — this round's score; completed — rounds done; total —
//   total rounds.
void ProgressStatus::SetBestScore(double score, int completed, int total) {
  std::lock_guard<std::mutex> lk(mu);
  if (!bestScore.hasValue || score < bestScore.score) {
    bestScore.score = score;
  }
  bestScore.hasValue = true;
  bestScore.completedRounds = completed;
  bestScore.totalRounds = total;
  bestScore.title.clear();
  bestScore.subtitle.clear();
}

// Clear the pinned best-score box.
void ProgressStatus::ResetBestScore() {
  std::lock_guard<std::mutex> lk(mu);
  bestScore = BestScoreView{};
}

// Set the best-score box to a labelled value (latest wins, not the min);
// used for the refining mix-mode Fst estimate.
// Params: score — value to show; title — box title; subtitle — second line.
void ProgressStatus::SetBestScoreLabeled(double score,
                                         const std::string& title,
                                         const std::string& subtitle) {
  std::lock_guard<std::mutex> lk(mu);
  bestScore.score = score;
  bestScore.hasValue = true;
  bestScore.title = title;
  bestScore.subtitle = subtitle;
  bestScore.completedRounds = 0;
  bestScore.totalRounds = 0;
}

// Take a thread-safe copy of the current progress state for rendering.
// Returns: a ProgressSnapshot with all progress/phase/chart fields copied.
ProgressSnapshot ProgressStatus::Snapshot() const {
  std::lock_guard<std::mutex> lk(mu);
  ProgressSnapshot snapshot;
  snapshot.totalTasks = totalTasks;
  snapshot.currentTask = currentTask;
  snapshot.totalSteps = totalSteps;
  snapshot.currentStep = currentStep;
  snapshot.phaseProgress = taskProgress;
  snapshot.globalProgress = globalProgress;
  snapshot.currentPhaseIndex =
      static_cast<std::size_t>(std::max(0.0, currentTask));
  snapshot.currentPhaseName = taskName;
  snapshot.statusDetail = statusDetail;
  snapshot.warning = warning;
  snapshot.phases = phases;
  snapshot.ne = ne;
  snapshot.bestScore = bestScore;
  return snapshot;
}

// Persist the current progress to the legacy text file.
void ProgressStatus::SaveProgress() {
  WriteLegacyProgressFile();
}

// Populate phases with generic "Task N" rows if none were configured.
void ProgressStatus::EnsurePhaseList() {
  if (!phases.empty()) return;
  const int count = static_cast<int>(totalTasks);
  phases.reserve(count);
  for (int i = 0; i < count; ++i) {
    ProgressPhaseView phase;
    phase.name = "Task " + std::to_string(i + 1);
    phases.push_back(phase);
  }
}

// Recompute each phase's Done/Current/Pending state and bar fill from the
// current task index and progress.
void ProgressStatus::RefreshPhaseStates() {
  const std::size_t current =
      static_cast<std::size_t>(std::max(0.0, currentTask));
  for (std::size_t i = 0; i < phases.size(); ++i) {
    phases[i].progress = 0.0;
    if (i < current) {
      phases[i].state = ProgressPhaseState::Done;
      phases[i].progress = 1.0;
    } else if (i == current) {
      phases[i].state = ProgressPhaseState::Current;
      phases[i].progress = taskProgress;
    } else {
      phases[i].state = ProgressPhaseState::Pending;
    }
  }
}

// Write the plain-text progress file (global %, completed phases, and the
// current phase %) that users can read via `cat <file>_GONE_progress.tmp`.
void ProgressStatus::WriteLegacyProgressFile() const {
  if (progressFname.empty()) return;
  std::ofstream progressFile(progressFname, std::ios::out);
  if (!progressFile.good()) return;
  progressFile << "    GLOBAL PROGRESS: " << std::fixed
               << std::setprecision(2) << globalProgress * 100.0 << "%\n";
  for (const ProgressPhaseView& phase : phases) {
    if (phase.state == ProgressPhaseState::Done) {
      progressFile << "    " << phase.name << " ... DONE\n";
    }
  }
  progressFile << "    " << taskName << " ... "
               << taskProgress * 100.0 << "%\n";
}

