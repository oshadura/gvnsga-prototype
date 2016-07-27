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

namespace geantvmoop {

/*! \brief Multi-objective optimization benchmark function CIGTAB 1.
*
*  The function is described in
*
*  Christian Igel, Nikolaus Hansen, and Stefan Roth.
*  Covariance Matrix Adaptation for Multi-objective Optimization.
*  Evolutionary Computation 15(1), pp. 1-28, 2007
*/

class CigTab01 : public Functions<CigTab01> {

public:
  /*
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  static Output Evaluate(const Input &individual) {
    rotationMatrix = blas::randomRotationMatrix(individual.size());
    ResultType y = prod(rotationMatrix, x);
    double result = sqr(y(0)) + sqr(m_a) * sqr(y(individual.size() - 1));
    for (unsigned i = 1; i < individual.size() - 1; i++) {
      result += m_a * sqr(y(i));
    }
    fFitness.insert(it, (result / ( sqr(m_a) * individual.size());

    result = sqr(y(0) - 2) + sqr(m_a) * sqr(y(individual.size() - 1) - 2);
    for (unsigned i = 1; i < individual.size() - 1; i++) {
      result += m_a * sqr(y(i) - 2);
    }
    fFitness.insert(it + 1, (result / ( sqr(m_a) * individual.size()));
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 4; ++i)
      vector.push_back(GADouble(-10, 10));
    return vector;
  }

  // Crap Here
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        std::sqrt(1 - parameter[0] * x[0] * x[0] - parameter[1] * x[1] * x[1] -
                  parameter[2] * x[2] * x[2]);
    return value;
  }

  static Output GetOutput() { return std::vector<double>(2); }
  */
};
}

#endif
