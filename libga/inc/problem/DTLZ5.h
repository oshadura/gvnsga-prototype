#pragma once

#ifndef __PROBLEMDTLZ5__
#define __PROBLEMDTLZ5__

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

class DTLZ5 : public Functions<DTLZ5> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  /*
  private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
      ar & boost::serialization::base_object<Functions<DTLZ5>>(*this);
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
    Int_t k = n - m + 1; // 5
    Double_t g = 0.0;
    std::vector<double> phi(m);
    for (std::size_t i = n - k + 1; i <= n; ++i) {
      g += std::pow(fParameters[i - 1] - 0.5, 2);
    }
    double t = pi() / (4 * (1 + g));
    phi[0] = fParameters[0] * pi() / 2;
    for (std::size_t i = 2; i <= m - 1; ++i) {
      phi[i - 1] = t * (1 + 2 * g * fParameters[i - 1]);
    }
    for (std::size_t i = 1; i <= m; ++i) {
      double f = 1 + g;
      for (std::size_t j = m - i; j >= 1; --j) {
        f *= std::cos(phi[j - 1]);
      }
      if (i > 1) {
        f *= std::sin(phi[m - i + 1] - 1);
      }
      // auto it = fFitness.begin();
      // fFitness.insert(it + i - 1, f);
      fFitness.push_back(f);
    }
//    std::cout << "Vector output for evaluation function: ";
//    for (auto i : fFitness)
//      std::cout << i << ' ';
//    std::cout << ' ' << std::endl;
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 7; ++i)
      vector.push_back(GADouble(0, 1));
    return vector;
  }

  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value = std::sqrt(parameter[0] - parameter[1] * x[0] * x[0] -
                               parameter[2] * x[1] * x[1]);
    return value;
  }

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
