#ifndef __POLYNOMIALMUTATION__
#define __POLYNOMIALMUTATION__

#include "Mutation.h"
#include "tools/Generator.h"

namespace geantvmoop{

class PolynomialMutation : public Mutation<PolynomialMutation> {

public:
  template <typename Individual> static void MutateGA(Individual &individual, double fProbability = -1) {
    if (fProbability == -1)
      fProbability = 1 / (double)individual.size();
    auto random = Generator::GetInstance();
    for (unsigned int i = 0; i < individual.size(); ++i) {
      if (random.RNGDouble() < fProbability)
        individual[i] = individual[i].RandomSetup();
    }
  }
};

}

#endif
