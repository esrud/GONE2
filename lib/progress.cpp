#include "./progress.hpp"

void ProgressStatus::InitTotalTasks(float nTasks, const char* fname) {
  totalTasks = nTasks;
  progressFname = fname;
}

void ProgressStatus::SetCurrentTask(float cTask, const char *tName) {
  if (cTask > totalTasks) {
    std::cerr << "ERROR: Setting progress larger than max value\n";
    return;
  }
  if (cTask != 0) {
    doneTasks.push_back(taskName);
  }
  currentTask = cTask;
  taskName = tName;
  taskProgress = 0.0;
  currentStep = 0.0;
  globalProgress = currentTask / totalTasks;
}

void ProgressStatus::InitCurrentTask(float nSteps) {
  totalSteps = nSteps;
}

void ProgressStatus::SetTaskProgress(float cStep) {
  if (cStep > totalSteps) {
    std::cerr << "ERROR: Setting progress bigger than max value\n";
    return;
  }
  currentStep = cStep;
  taskProgress = currentStep / totalSteps;
  globalProgress = (currentTask + taskProgress) / totalTasks;
  SaveProgress();
}

void ProgressStatus::SaveProgress() {
  // Save current progress to a file
  std::ofstream progressFile;
  progressFile.open(progressFname, std::ios::out);
  if (!progressFile.good()) {
    // TODO(me): maybe print some text? But then it will print every
    // time we save the progress
    return;
  }
  progressFile << "    GLOBAL PROGRESS: " << std::fixed << std::setprecision(2)
               << globalProgress * 100.0 << "%\n";
  for (const auto &task : doneTasks) {
    progressFile << "    " << task << " ... DONE\n";
  }
  progressFile << "    " << taskName << " ... " << taskProgress * 100.0 << "%\n";
  progressFile.close();
}

void ProgressStatus::PrintProgress() {
  std::ifstream progressFile;
  progressFile.open(progressFname, std::ifstream::in);

  if (!progressFile.good()) {
    return;
  }
  progressFile.seekg(0, progressFile.end);
  int length = progressFile.tellg();
  progressFile.seekg(0, progressFile.beg);
  char *pStatus = new char[length+1]{0};
  progressFile.read(pStatus, length);
  std::cout << pStatus;
  delete[] pStatus;
  progressFile.close();
}

