#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

class ProgressStatus {
 private:
  float totalTasks;
  float currentTask = 0;

  float totalSteps;
  float currentStep = 0;

  float taskProgress;
  float globalProgress;

  const char *taskName;
  const char *progressFname;

  std::vector<const char*> doneTasks;

 public:
  void InitTotalTasks(float nTasks, const char* fname);
  void SetCurrentTask(float cTask, const char* tName);
  void InitCurrentTask(float nSteps);
  void SetTaskProgress(float cStep);
  void SaveProgress();
  void PrintProgress();
};
