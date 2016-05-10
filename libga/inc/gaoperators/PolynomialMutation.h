#ifndef __POLYNOMIALMUTATION__
#define __POLYNOMIALMUTATION__

#include "Mutation.h"
#include "tools/Generator.h"

class PolynomialMutation : public Mutation<PolynomialMutation> {

public:
  template <typename T> static void MutateGA(T &individual, double prob = -1) {
    if (prob == -1)
      prob = 1 / (double)individual.size();
    auto random = Generator::GetInstance();
    for (unsigned int i = 0; i < individual.size(); ++i) {
      if (random.RNGDouble() < prob)
        individual[i] = individual[i].random();
    }
  }
};

#endif
