#ifndef __PROBLEMDTLZ2__
#define __PROBLEMDTLZ2__

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

class DTLZ2 : public Functions<DTLZ2> {

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
    
    std::cout << "Vector input for evaluation function: ";
    for (auto i: fParameters)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    
    int n = 12;
    int m = 3;
    int k = n - m + 1; // 10
    double g = 0.0;
    for (int i = m - 1; i < n; ++i) {
      g += pow(fParameters[i] - 0.5, 2);
    }
    for (int i = 0; i < m; ++i) {
      Double_t f = (1 + g);
      size_t j = 0;
      for (; i + m <= m - 2; ++j) {
        f *= cos(fParameters[j] * pi()/2);
      }
      if (m > 0) {
        f *= sin(fParameters[j] * pi()/2);
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

  static Output GetOutput() { return std::vector<double>(3); }
};
}

#endif
