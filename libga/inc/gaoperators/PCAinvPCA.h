#pragma once

#ifndef __PCAINVPCA__
#define __PCAINVPCA__

#include "gaoperators/GANoiseReduction.h"
#include "mva/PCA.h"
#include "mva/LPCA.h"
#include "mva/KPCA.h"
#include "mva/RobustPCA.h"
#include "mva/UncenteredLPCA.h"

namespace geantvmoop {

// TBD templated!
class PCAinvPCA : public NoiseReduction<PCAinvPCA> {

public:
  /*
  template <typename F>
  Population<F> NoiseReductionImpl(Population<F> &population) {
    Population<F> result;
    LPCA lpca;
    result = lpca.MVA(population);
    return result;
  }
  */

  /*
    template <typename F>
    Population<F> NoiseReductionImpl(Population<F> &population) {
      Population<F> result;
      RobustPCA robustpca;
      result = robustpca.MVA(population);
      return result;
    }
  */

  template <typename F>
  Population<F> NoiseReductionImpl(Population<F> &population) {
    Population<F> result;
    UncenteredLPCA ulpca;
    result = ulpca.MVA(population);
    return result;
  }
};
}

#endif
