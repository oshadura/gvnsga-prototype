#pragma once

#ifndef __PROBLEMSUFC9__
#define __PROBLEMSUFC9__

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

class UFC9 : public Functions<UFC9> {

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
    unsigned int j, count1, count2, count3;
    double sum1, sum2, sum3, yj, E;
    E = 0.1;
    sum1 = sum2 = sum3 = 0.0;
    count1 = count2 = count3 = 0;
    for (j = 3; j <= fParameters.size(); j++) {
      yj = fParameters[j - 1] -
           2.0 * fParameters[1] *
               sin(2.0 * pi() * fParameters[0] + j * pi() / fParameters.size());
      if (j % 3 == 1) {
        sum1 += yj * yj;
        count1++;
      } else if (j % 3 == 2) {
        sum2 += yj * yj;
        count2++;
      } else {
        sum3 += yj * yj;
        count3++;
      }
    }
    yj = (0.5 + E) * (1.0 - 4.0 * (2.0 * fParameters[0] - 1.0) *
                                (2.0 * fParameters[0] - 1.0));
    if (yj < 0.0)
      yj = 0.0;
    auto it = fFitness.begin();
    fFitness.insert(it, (0.5 * (yj + 2 * fParameters[0]) * fParameters[1] +
                         2.0 * sum1 / (double)count1));
    fFitness.insert(it + 1,
                    (0.5 * (yj - 2 * fParameters[0] + 2.0) * fParameters[1] +
                     2.0 * sum2 / (double)count2));
    fFitness.insert(it +2, (1.0 - fParameters[1] + 2.0 * sum3 / (double)count3));
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 30; ++i)
      vector.push_back(GADouble(0, 1));
    return vector;
  }

  // Crap here!
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
