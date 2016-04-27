#ifndef PROBLEMDTLZ_H
#define PROBLEMDTLZ_H

#include "generic/TGenes.h"
#include "generic/Problem.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "output/HistogramManager.h"
#include "algorithms/AlgorithmNSGA.h"
#include "instrumentation/GeantVFitness.h"

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <utility>

class ProblemDTLZ1 : public Problem<ProblemDTLZ1> {

public:
  typedef Genes<double> InputType;
  typedef std::vector<double> OutputType;
  constexpr static const double pi = boost::math::constants::pi<double>();

  static OutputType Evaluate(InputType &individual) {
    Int_t n = individual.GetSetup()->GetNParam();      // 7
    Int_t m = individual.GetSetup()->GetNObjectives(); // 3
    Int_t k = n - m + 1;                               // 5

    Double_t g = 0.0;

    for (Int_t i = m - 1; i < n; ++i) {
      g += pow(individual.GetGene(i - 1) - 0.5, 2) -
           cos(20 * pi * (individual.GetGene(i - 1) - 0.5));
    }
    g = 100 * (k + g);
    for (Int_t i = 0; i < m; ++i) {
      Double_t f = 0.5 * (1 + g);
      size_t j = 0;
      for (; m >= 2 + i && j <= m - 2 - i; ++j) {
        f *= individual.GetGene(j);
      }
      if (i > 0) {
        f *= (1 - individual.GetGene(j));
      }
      individual.SetFitness(i, f);
    }
    return individual.GetFitnessVector();
  }

  //static InputType GetInput() { return 0; }

  //static OutputType GetOutput() { return 0; }
};

#endif