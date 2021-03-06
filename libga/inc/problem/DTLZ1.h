#pragma once

#ifndef __PROBLEMDTLZ1__
#define __PROBLEMDTLZ1__

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

class DTLZ1 : public Functions<DTLZ1> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  DTLZ1() : fNGenes(7), fNObjectives(3) {}

  virtual ~DTLZ1() {}
  /*
  private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
      ar & boost::serialization::base_object<Functions<DTLZ1>>(*this);
    }

  public:
    */

  void SetNGenes(const int i) { fNGenes = i; };

  void SetNObjectives(const int i) { fNObjectives = i; };

  int GetNGenes() const { return fNGenes; };

  //  int GetNObjectives() const { return fNObjectives; };

private:
  int fNGenes, fNObjectives;

public:
  static double pi() { return std::atan(1) * 4; }

  static Output Evaluate(const Input &individual) {
    std::vector<double> fFitness, fParameters;
    fFitness.reserve(3);
    fParameters.reserve(individual.size());
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
    std::cout << "Vector input for evaluation function: ";
    for (auto i : fParameters)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    int n = 7;
    int m = 3;
    int k = n - m + 1; // 5
    double g = k;
    for (std::size_t w = /*n - k*/ m; w <= n; w++)
      g += (k + pow(fParameters[w] - 0.5, 2) -
           std::cos(20.0 * pi() * (fParameters[w] - 0.5)));
    g *= 100;
    for (std::size_t i = 0; i < m; i++) {
      double f = 0.5 * (1.0 + g);
      for (std::size_t j = 0; j < m - i - 1; ++j)
        f *= fParameters[j];
      if (i > 0)
        f *= 1 - fParameters[m - i - 1];
      auto it = fFitness.begin();
      fFitness.insert(it + i, f);
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

  // ROOT Fitting to true Pareto front
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value = parameter[0] * x[0] + parameter[1] * x[1] - parameter[2];
    return value;
  }

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
