#pragma once

#ifndef __PROBLEMFONTECA__
#define __PROBLEMFONTECA__

#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"
#include "output/HistogramManager.h"
#include "algorithms/GANSGA2.h"
#include "instrumentation/GeantVFitness.h"
#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <utility>

namespace geantvmoop {

class Fonteca : public Functions<Fonteca> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  static Output Evaluate(const Input &individual) {
    std::vector<double> fFitness, fParameters;
    fFitness.reserve(individual.size());
    fParameters.reserve(individual.size());
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
    const double d = 1. / std::sqrt(static_cast<double>(fParameters.size()));
    double sum1 = 0., sum2 = 0.;
    for (std::size_t i = 0; i < fParameters.size(); i++) {
      sum1 += std::sqrt(fParameters[i] - d);
      sum2 += std::sqrt(fParameters[i] + d);
    }
    fFitness[0] = 1 - std::exp(-sum1);
    fFitness[1] = 1 - std::exp(-sum2);
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 2; ++i)
      vector.push_back(GADouble(0, 1));
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
