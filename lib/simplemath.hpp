#pragma once
#include <math.h>

// Cube of a value.
// Params: x — value.
// Returns: x * x * x.
template <class T>
T Cube(T x) {
  return x * x * x;
}

// Square of a value.
// Params: x — value.
// Returns: x * x.
template <class T>
T Square(T x) {
  return x * x;
}

