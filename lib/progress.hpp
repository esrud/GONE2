#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

enum class ProgressPhaseState {
  Pending,
  Current,
  Done
};

struct ProgressPhaseView {
  std::string name;
  ProgressPhaseState state = ProgressPhaseState::Pending;
  double progress = 0.0;
};

struct NeEstimateView {
  bool hasValues = false;
  std::vector<double> values;
  std::string unitLabel;
  std::vector<double> referenceValues;
  std::vector<double> minValues;
  std::vector<double> maxValues;
  std::vector<std::vector<double>> densitySamples;
};

struct BestScoreView {
  bool hasValue = false;
  double score = 0.0;
  int completedRounds = 0;
  int totalRounds = 0;
  std::string title;
  std::string subtitle;
};

struct ProgressSnapshot {
  double totalTasks = 1.0;
  double currentTask = 0.0;
  double totalSteps = 1.0;
  double currentStep = 0.0;
  double phaseProgress = 0.0;
  double globalProgress = 0.0;
  std::size_t currentPhaseIndex = 0;
  std::string currentPhaseName;
  std::string statusDetail;
  std::string warning;
  std::vector<ProgressPhaseView> phases;
  NeEstimateView ne;
  BestScoreView bestScore;
};

class ProgressStatus {
 private:
  double totalTasks = 1.0;
  double currentTask = 0.0;
  double totalSteps = 1.0;
  double currentStep = 0.0;
  double taskProgress = 0.0;
  double globalProgress = 0.0;
  double taskWindowBase = 0.0;
  double taskWindowSpan = 1.0;
  double taskProgressFloor = 0.0;
  std::string taskName;
  std::string progressFname;
  std::string statusDetail;
  std::string warning;
  std::vector<ProgressPhaseView> phases;
  NeEstimateView ne;
  BestScoreView bestScore;
  mutable std::mutex mu;

  void EnsurePhaseList();
  void RefreshPhaseStates();
  void WriteLegacyProgressFile() const;

 public:
  void ConfigurePhases(const std::vector<std::string>& phaseNames);
  void InitTotalTasks(float nTasks, const char* fname);
  void SetCurrentTask(float cTask, const char* tName);
  void InitCurrentTask(float nSteps);
  void SetTaskProgress(float cStep);
  void SetTaskWindow(double base, double span);
  void SetStatusDetail(const std::string& detail);
  void SetWarning(const std::string& message);
  void SetNeSnapshot(const std::vector<double>& values,
                     const std::string& unitLabel);
  void SetReferenceNe(const std::vector<double>& values);
  void SetNeBand(const std::vector<double>& minValues,
                 const std::vector<double>& maxValues);
  void SetNeDensity(const std::vector<std::vector<double>>& samples);
  void SetBestScore(double score, int completed, int total);
  void ResetBestScore();
  void SetBestScoreLabeled(double score, const std::string& title,
                           const std::string& subtitle);
  ProgressSnapshot Snapshot() const;
  void SaveProgress();
};
