#pragma once

#ifndef __PROBLEMDTLZ5__
#define __PROBLEMDTLZ5__

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

class DTLZ5 : public Functions<DTLZ5> {

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
    Int_t k = n - m + 1; // 5
    Double_t g = 0.0;
    for (Int_t i = m - 1; i < n; ++i) {
      g += pow(fParameters[i - 1] - 0.5, 2) -
           cos(20 * pi() * (fParameters[i - 1] - 0.5));
    }
    g = 100 * (k + g);
    for (Int_t i = 0; i < m; ++i) {
      Double_t f = 0.5 * (1 + g);
      size_t j = 0;
      for (; m >= 2 + i && j <= m - 2 - i; ++j) {
        f *= fParameters[j];
      }
      if (i > 0) {
        f *= (1 - fParameters[j]);
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
