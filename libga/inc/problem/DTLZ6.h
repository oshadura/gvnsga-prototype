#pragma once

#ifndef __PROBLEMDTLZ6__
#define __PROBLEMDTLZ6__

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

class DTLZ6 : public Functions<DTLZ6> {

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
  
    int n = 7;
    int m = 3;
    int k = n - m + 1; // 5
    Double_t g = 0.0;
    for (int i = m - 1; i < n; ++i) {
    g += pow(fParameters[i], 0.1);
  }
  std::vector<double> theta(n);
  theta[0] = fParameters[0] * pi()/2;
  for (int i = 1; i < theta.size(); ++i) {
    theta[i] = pi()/2 / (2 * (1 + g)) * (1 + 2 * g * fParameters[i]);
  }
  for (int i = 0; i < m; ++i) {
    Double_t f = (1 + g);
    int j = 0;
    for (; i + m <= m - 2; ++j) {
      f *= cos(theta[j]);
    }
    if (m > 0) {
      f *= sin(theta[j]);
    }
      fFitness.push_back(f);
    }
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 7; ++i)
      vector.push_back(GADouble(0, 1));
    return vector;
  }

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
