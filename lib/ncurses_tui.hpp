#pragma once

class ProgressStatus;

// Start the ncurses-backed terminal UI. Spawns a background thread
// that polls `progress->Snapshot()` and redraws at ~10 Hz. Idempotent:
// calling twice without a Stop in between is a no-op the second time.
void StartNcursesTui(ProgressStatus* progress);

// Block the calling thread until the user presses Ctrl-C. Installs a
// SIGINT handler that just sets a flag; the renderer thread already
// blocks SIGINT so the signal is reliably delivered here. Used at
// the end of a run to keep the final chart on screen until the user
// dismisses it.
void WaitForSigint();

// Signal the background thread to exit, tear down ncurses, and
// restore the terminal. Blocks until the renderer thread joins.
void StopNcursesTui();
