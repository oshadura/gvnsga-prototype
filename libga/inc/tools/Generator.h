#ifndef GENERATOR_H
#define GENERATOR_H

#include <stdlib.h>
#include <time.h>
#include <random>

namespace geantvmoop {

class Generator {
public:
  static Generator &GetInstance() { return Instance; }
  Generator() {}
  Generator(Generator const &);
  void operator=(Generator const &);
  ~Generator() {}
  double RNGBool();
  double RNGDouble();
  double RNGDouble(int fMin, int fMax) {
    return fMin + (double)std::rand() / RAND_MAX * (fMax - fMin);
  }
  std::vector<double> RNGDoubleVector(int fMin, int fMax);
  int RNGInteger();
  int RNGInteger(int fMin, int fMax);
  std::vector<int> RNGIntVector(int fMin, int fMax);
  std::vector<int> RNGVector();
  // Special case for 2^x rnd generator
  int RNGSIMD(int fMin, int fMax);
  std::vector<int> RNGSIMDVector(int fMin, int fMax);

private:
  static Generator Instance;
};
}

#endif
