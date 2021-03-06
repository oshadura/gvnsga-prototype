#pragma once

#ifndef __PROBLEMDTLZ7__
#define __PROBLEMDTLZ7__

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

class DTLZ7 : public Functions<DTLZ7> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  /*
  private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
      ar & boost::serialization::base_object<Functions<DTLZ7>>(*this);
    }

  public:
    */

  static double pi() { return std::atan(1) * 4; }

  static Output Evaluate(const Input &individual) {
    std::vector<double> fFitness, fParameters;
    fFitness.reserve(individual.size());
    fParameters.reserve(individual.size());
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
    int n = 11;
    int m = 3;
    Int_t k = n - m + 1;
    double g = 0.0;
    auto it = fFitness.begin();
    for (int i = n - k + 1; i <= m; ++i) {
      fFitness.insert(it + i, fParameters[i]);
    }
    for (int i = n - k + 1; i <= n; ++i) {
      g += fParameters[i - 1];
    }
    g = 1 + 9 * g / k;
    double h = 0;

    for (int j = 1; j <= m - 1; ++j) {
      h += fParameters[j - 1] / (1 + g) *
           (1 + std::sin(3 * pi() * fParameters[j - 1]));
    }
    fFitness.insert(it + m - 1, (1 + g) * h);
    std::cout << "Vector output for evaluation function: ";
    for (auto i : fFitness)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 11; ++i)
      vector.push_back(GADouble(0, 1));
    return vector;
  }

  // Crap
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
