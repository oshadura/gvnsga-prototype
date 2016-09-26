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

/*
private:
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & boost::serialization::base_object<Functions<DTLZ3>>(*this);
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
    int n = 7;
    int m = 3;
    int k = n - m + 1; // 5
    double g = k;
    for (std::size_t i = n - k; i < sizeof(n); i++)
      g += pow(fParameters[i] - 0.5, 2) -
           std::cos(20.0 * pi() * (fParameters[i] - 0.5));
    g *= 100;
    for (std::size_t i = 0; i < sizeof(m); i++) {
      double f = 1.0 + g;
      for (std::size_t j = 0; j < sizeof(m - i - 1); ++j)
        f *= std::cos(fParameters[j]) * pi() / 2;
      if (i > 0)
        f *= std::sin(fParameters[m - i - 1]) * pi() / 2;
      auto it = fFitness.begin();
      fFitness.insert(it + i, f);
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
