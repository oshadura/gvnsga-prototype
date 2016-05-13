#ifndef MOO_SINGLEPOINTCROSSOVER_H
#define MOO_SINGLEPOINTCROSSOVER_H

#include "Crossover.h"
#include "tools/Generator.h"

namespace geantvmoop {

class SimpleCrossover : public Crossover<SimpleCrossover> {

public:
  template <typename T>
  static T CrossoverGA(const T &fIndividual1, const T &fIndividual2, int fPoint = -1) {
    T fOffspring = fIndividual1;
    if (fPoint < 0 || fPoint > fIndividual1.size())
      fPoint = Generator::GetInstance()->RNGInteger(0, fIndividual2.size());
    for (int i = fPoint; i < fIndividual2.size(); ++i) {
      fOffspring[i] = fIndividual2[i];
    }
    return fOffspring;
  }
};
}
