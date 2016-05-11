#ifndef __POLYNOMIALMUTATION__
#define __POLYNOMIALMUTATION__

#include "Mutation.h"
#include "tools/Generator.h"

class PolynomialMutation : public Mutation<PolynomialMutation> {

public:
  template <typename T> static T MutateGA(T &individual, double fProbability = -1) {
    if (fProbability == -1)
      fProbability = 1 / (double)individual.size();
    auto random = Generator::GetInstance();
    for (unsigned int i = 0; i < individual.size(); ++i) {
      if (random.RNGDouble() < fProbability)
        individual[i] = individual[i].RandomSetup();
    }
  }
};

#endif
