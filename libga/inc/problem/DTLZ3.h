#pragma once

#ifndef __PROBLEMDTLZ3__
#define __PROBLEMDTLZ3__

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

class DTLZ3 : public Functions<DTLZ3> {

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

    int n = 7;
    int m = 3;
    int k = n - m + 1; // 5
    double g = 0.0;
    for (int i = n - k + 1; i <= n; i++) {
      g += pow(fParameters[i - 1] - 0.5, 2) -
           cos(20 * pi() * (fParameters[i - 1] - 0.5));
    }
    g = 100 * (k + g);
    for (int i = 1; i <= m; i++) {
      double f = (1 + g);
      for (int j = m - i; j >= 1; j--) {
        f *= cos(fParameters[j - 1] * pi() / 2);
      }
      if (i > 1) {
        f *= sin(fParameters[m - i + 1] - 1) * pi() / 2;
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

  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        std::sqrt(1 - parameter[0] * x[0] * x[0] - parameter[1] * x[1] * x[1] -
                  parameter[2] * x[2] * x[2]);
    return value;
  }

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
