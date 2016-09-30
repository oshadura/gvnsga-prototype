#pragma once

#ifndef __PROBLEMKURSAWE__
#define __PROBLEMKURSAWE__

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

class Kursawe : public Functions<Kursawe> {

public:
  typedef GAVector<GADouble, 2> Input;

  typedef boost::container::static_vector<double, 3> Output;

/*
private:
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & boost::serialization::base_object<Functions<Kursawe>>(*this);
  }
  
public:
*/

  static Output Evaluate(const Input &individual) {
    boost::container::static_vector<double, 3> fFitness;
    boost::container::static_vector<double, 2> fParameters;
    fFitness.reserve(3);
    fParameters.reserve(2);
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
    double aux, xi, xj;
    fFitness[0] = 0.0;
    for (std::size_t var = 0; var < individual.size() - 1; var++) {
      xi = fParameters[var] * fParameters[var];
      xj = fParameters[var + 1] * fParameters[var + 1];
      aux = (-0.2) * sqrt(xi + xj);
      fFitness[0] += (-10.0) * exp(aux);
    }
    fFitness[1] = 0.0;
    for (std::size_t var = 0; var < individual.size(); var++) {
      fFitness[1] += pow(fabs(fParameters[var]), 0.8) +
                 5.0 * sin(pow(fParameters[var], 3.0));
    }
    auto it = fFitness.begin();
    fFitness.insert(it, fFitness[0]);
    fFitness.insert(it + 1, fFitness[2]);
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 5; ++i)
      vector.push_back(GADouble(-5, 5));
    return vector;
  }

  // No true
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        std::sqrt(-parameter[0] - parameter[0] * x[0] - parameter[1] * x[1]);
    return value;
  }

  static Output GetOutput() { return boost::container::static_vector<double, 3>(3); }
};
}

#endif
