#pragma once

#ifndef __PROBLEMDTLZ4__
#define __PROBLEMDTLZ4__

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

class DTLZ4 : public Functions<DTLZ4> {

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

    int alpha = 100;
    int n = 12;
    int m = 3;
    int k = n - m + 1; // 5
    Double_t g = 0.0;
    for (int i = m - 1; i < n; ++i) {
      g += pow(fParameters[i] - 0.5, 2);
    }
    for (int i = 0; i < m; ++i) {
      double f = (1 + g);
      int j = 0;
      for (; j + m <= m - 2; ++j) {
        f *= cos(pow(fParameters[j], alpha) * pi() / 2);
      }
      if (i > 0) {
        f *= sin(pow(fParameters[j], alpha) * pi() / 2);
      }
      fFitness.push_back(f);
    }
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 12; ++i)
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
