//===--- CMAES.h - LibGA -------------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Nsga2.h
 * @brief Implementation of CMAES algorithms for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __CMAES__
#define __CMAES__

#include "generic/GAAlgorithm.h"
//#include "generic/TGenes.h"
//#include "generic/Population.h"
//#include "generic/PF.h"
#include "cmaes.h"

using namespace libcmaes;

namespace geantvmoop {

template <typename F> class CMAES : public GAAlgorithm<CMAES<F>, F> {

public:
  CMAES(F problem) : GAAlgorithm<CMAES<F>, F>(problem) {}

  void Initialize() {
    int dim = 10; // problem dimensions.
    std::vector<double> x0(dim, 1.0);
    double sigma = 0.1;
    // int lambda = 100; // offsprings at each generation.
    GenoPheno<> gp(genof, phenof);
    CMAParameters<> cmaparams(x0, sigma, -1, 0, gp);
  }

  void Evolution() { CMASolutions cmasols = cmaes<>(F, cmaparams); }

  void Print(std::ostream &os) { os << fGen << std::endl; }

  // PF<F> GetParetoFront() { return fFront; }

private:
  Population<F> pop;

  PF<F> fFront;

  int fGen;
};
}

#endif
