#pragma once

#ifndef __PROBLEMSHWEFEL__
#define __PROBLEMSHWEFEL__

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

// Only 2 parameters.....

namespace geantvmoop {

class Shwefel : public Functions<Shwefel> {

public:
  typedef GAVector<GADouble, 3> Input;

  typedef boost::container::static_vector<double, 3> Output;

/*
private:
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & boost::serialization::base_object<Functions<Shwefel>>(*this);
  }
  
public:
*/

  static Output Evaluate(const Input &individual) {
    boost::container::static_vector<double, 3> fFitness;
    boost::container::static_vector<double, 3> fParameters;
    fFitness.reserve(3);
    fParameters.reserve(3);
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
    auto it = fFitness.begin();
    fFitness.insert(it, fParameters[0] * fParameters[0]);
    fFitness.insert(it + 1, ((fParameters[1] - 2) * (fParameters[1] - 2)));
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 2; ++i)
      vector.push_back(GADouble(-100, 100));
    return vector;
  }

  // Crap Here
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        std::sqrt(1 - parameter[0] * x[0] * x[0] - parameter[1] * x[1] * x[1] -
                  parameter[2] * x[2] * x[2]);
    return value;
  }

  static Output GetOutput() { return boost::container::static_vector<double, 2>(2); }
};
}

#endif
