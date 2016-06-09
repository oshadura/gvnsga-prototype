#pragma once

#ifndef __PROBLEMDTLZ7__
#define __PROBLEMDTLZ7__

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

class DTLZ7 : public Functions<DTLZ7> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  static double pi() { return std::atan(1) * 4; }

  static Output Evaluate(const Input &individual) {
    std::vector<double> fFitness, fParameters;
    fFitness.reserve(individual.size());
    fParameters.reserve(individual.size());
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());

    std::cout << "Vector input for evaluation function: ";
    for (auto i: fParameters)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;

    int n = 30;
    int m = 3;
    Int_t k = n - m + 1;
    Double_t g = 0.0;
    for (Int_t i = 0; i < m - 1; ++i) {
      fFitness[i] = fParameters[i];
    }
    for (Int_t i = m - 1; i < n; ++i) {
      g += fParameters[i];
    }
    g = 1 + 9 * g / k;
    Double_t h = m;
    for (Int_t j = 0; j < m - 1; ++j) {
      h -= fParameters[j] / (1 + g) * (1 + sin(3 * pi() * fParameters[j]));
    }
    fFitness[m - 1] = (1 + g) * h;
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 30; ++i)
      vector.push_back(GADouble(0, 1));
    return vector;
  }

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
