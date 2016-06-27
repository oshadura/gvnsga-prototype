#pragma once

#ifndef __PROBLEMDTLZ2__
#define __PROBLEMDTLZ2__

#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"
#include "output/HistogramManager.h"
#include "algorithms/GANSGA2.h"
#include "instrumentation/GeantVFitness.h"
#include <boost/math/constants/constants.hpp>

#include "TMath.h"

#include <cmath>
#include <utility>

namespace geantvmoop {

class DTLZ2 : public Functions<DTLZ2> {

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
    /*
    std::cout << "Vector input for evaluation function: ";
    for (auto i : fParameters)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    */
    int n = 12;
    int m = 3;
    int k = n - m + 1; // 10
    double g = 0.0;
    for (int i = n - k; i < n; ++i) {
      g += std::pow(fParameters[i] - 0.5, 2);
    }
    for (int i = 0; i < m; ++i) {
      double f = 1 + g;
      for (std::size_t j = 0; j < m - i - 1; ++j) {
        f *= std::cos(fParameters[j] * pi() / 2);
      }
      if (i > 0) {
        f *= std::sin(fParameters[m - i - 1] * pi() / 2);
      }
      auto it = fFitness.begin();
      fFitness.insert(it + i, f);
    }
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 12; ++i)
      vector.push_back(GADouble(0, 1));
    return vector;
  }

  // ROOT Fitting to true Pareto front
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        std::sqrt(parameter[0] - parameter[1] * x[0] * x[0] - parameter[2] * x[1] * x[1] /*-
                  parameter[2] * x[2] * x[2]*/);
    return value;
  }

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
