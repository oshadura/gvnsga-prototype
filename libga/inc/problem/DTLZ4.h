#pragma once

#ifndef __PROBLEMDTLZ4__
#define __PROBLEMDTLZ4__

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

class DTLZ4 : public Functions<DTLZ4> {

public:
  typedef GAVector<GADouble> Input;
  typedef std::vector<double> Output;

  /*
  private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
      ar & boost::serialization::base_object<Functions<DTLZ4>>(*this);
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
    int alpha = 100;
    int n = 12;
    int m = 3;
    int k = n - m + 1; // 10
    double g = 0.0;
    for (std::size_t i = /* n - k*/ m; i <= n; ++i) {
      g += std::pow(fParameters[i] - 0.5, 2);
    }
    for (int i = 0; i < m; ++i) {
      double f = 1 + g;
      for (std::size_t j = 0; j < m - i - 1; ++j) {
        f *= std::cos(std::pow(fParameters[j], alpha) * pi() / 2);
      }
      if (i > 0) {
        f *= std::sin(std::pow(fParameters[m - i - 1], alpha) * pi() / 2);
      }
      auto it = fFitness.begin();
      fFitness.insert(it + i, f);
    }
  //  std::cout << "Vector output for evaluation function: ";
  //  for (auto i : fFitness)
  //    std::cout << i << ' ';
  //  std::cout << ' ' << std::endl;
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
