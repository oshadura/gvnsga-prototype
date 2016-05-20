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
#ifdef ENABLE_CMAES

#ifndef __CMAES__
#define __CMAES__

#include "generic/GAAlgorithm.h"
//#include "generic/TGenes.h"
//#include "generic/Population.h"
//#include "generic/PF.h"
#include "cmaes.h"

namespace geantvmoop {

template <typename F> class CMAES : public GAAlgorithm<CMAES<F>, F> {

public:
  CMAES(F problem) : GAAlgorithm<CMAES<F>, F>(problem) {}

  void InitializeImpl() {
    int dim = 20; // problem dimensions.
    std::vector<double> x0(dim, 1.0);
    double sigma = 0.1;
    // int lambda = 100; // offsprings at each generation.
    libcmaes::CMAParameters<> cmaparams(x0, sigma);
    cmaparams.set_fplot("cmaes.dat");
    // TBD: LibCMaes Fitfunc -> LibGa Function!
    // libcmaes::CMASolutions cmasols = libcmaes::cmaes<>(
    //    F, cmaparams, libcmaes::CMAStrategy<libcmaes::CovarianceUpdate>::_defaultPFunc,
    //    nullptr, libcmaes::CMASolutions(), PlotFunction);
  }

  /////////////////////////////
  // Plotting test
  libcmaes::PlotFunc<libcmaes::CMAParameters<>, libcmaes::CMASolutions>
  PlotFunction = [](const libcmaes::CMAParameters<> &cmaparams,
                    const libcmaes::CMASolutions &cmasols,
                    std::ofstream &fplotstream) {
    fplotstream << "Kappa = " << cmasols.max_eigenv() / cmasols.min_eigenv()
                << std::endl; // storing covariance matrix condition number to
                              // file.
    return 0;
  };
  /////////////////////////////

  void EvolutionImpl() {/*Nothing to do*/}

  void PrintImpl(std::ostream &os) { os << fGen << std::endl; }

  PF<F> GetParetoFrontImpl() {/* Nothing to do*/}

private:
  Population<F> pop;

  PF<F> fFront;
};
}

#endif

#endif
