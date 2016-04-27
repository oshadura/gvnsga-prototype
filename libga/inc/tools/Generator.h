#ifndef GENERATOR_H
#define GENERATOR_H

#include <stdlib.h>
#include <time.h>

// Re-think templated...

class Generator {

public:
  void GeneratorDouble();
  void GeneratorDouble(int fMin, int fMax);
  void GeneratorDoubleVector(int fMin, int fMax);
  void GeneratorInt();
  void GeneratorInt(int fMin, int fMax);
  void GeneratorIntVector(int fMin, int fMax);
  void GeneratorVector();
  void GeneratorSIMD(int fMin, int fMax);
  void GeneratorSIMDVector(int fMin, int fMax);
};

#endif