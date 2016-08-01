#pragma once

#ifndef __PROBLEMBINH__
#define __PROBLEMBINH__

#include "algorithms/GANSGA2.h"
#include "generic/Functions.h"
#include "generic/GADouble.h"
#include "generic/GAVector.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "output/HistogramManager.h"
#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <utility>

namespace geantvmoop {

class Binh : public Functions<Binh> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  // We need to add possibility to get constrained data generation,
  // here is: https://en.wikipedia.org/wiki/Test_functions_for_optimization

  static Output Evaluate(const Input &individual) {
    std::vector<double> fFitness, fParameters;
    fFitness.reserve(individual.size());
    fParameters.reserve(individual.size());
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
    auto it = fFitness.begin();
    fFitness.insert(it, (4 * fParameters[0] * fParameters[0] +
                         4 * fParameters[1] * fParameters[1]));
    fFitness.insert(it + 1, ((fParameters[0] - 5) * (fParameters[0] - 5) +
                             (fParameters[1] - 5) * (fParameters[1] - 5)));
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 2; ++i)
      vector.push_back(GADouble(-100, 100));
    return vector;
  }

  // Crap Here
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        std::sqrt(1 - parameter[0] * x[0] * x[0] - parameter[1] * x[1] * x[1] -
                  parameter[2] * x[2] * x[2]);
    return value;
  }

  static Output GetOutput() { return std::vector<double>(2); }
};
}

#endif
