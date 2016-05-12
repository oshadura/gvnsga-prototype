#ifndef GENERATOR_H
#define GENERATOR_H

#include <stdlib.h>
#include <time.h>

namespace geantvmoop {

class Generator {
public:
  static Generator &GetInstance() {
    static Generator Instance;
    return Instance;
  }

  Generator() {}
  ~Generator() {}
  // Generator(Generator const&) = default;
  // void operator=(Generator const&) = default;

public:
  double RNGBool();
  double RNGDouble();
  double RNGDouble(int fMin, int fMax) {
    return fMin + (double)std::rand() / RAND_MAX * (fMax - fMin);
  }
  std::vector<double> RNGDoubleVector(int fMin, int fMax);
  int RNGInt();
  int RNGInt(int fMin, int fMax);
  std::vector<int> RNGIntVector(int fMin, int fMax);
  std::vector<int> RNGVector();
  // Special case for 2^x rnd generator
  int RNGSIMD(int fMin, int fMax);
  std::vector<int> RNGSIMDVector(int fMin, int fMax);
};
}

#endif