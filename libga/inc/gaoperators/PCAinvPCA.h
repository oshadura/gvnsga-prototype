#pragma once

#ifndef __PCAINVPCA__
#define __PCAINVPCA__

#include "gaoperators/GANoiseReduction.h"
#include "mva/PCA.h"
#include "mva/LPCA.h"
#include "mva/KPCA.h"
#include "mva/RobustPCA.h"
#include "mva/RobustTrickPCA.h"
#include "mva/UncenteredLPCA.h"
#include "mva/UncenteredTrickLPCA.h"
#include "mva/SVDRepresentation.h"

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
    template <typename F>
    Population<F> NoiseReductionImpl(Population<F> &population) {
      Population<F> result;
      RobustPCA robustpca;
      result = robustpca.MVA(population);
      return result;
    }

        template <typename F>
    Population<F> NoiseReductionImpl(Population<F> &population) {
      Population<F> result;
      RobustTrickPCA robustpca;
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
/*
  template <typename F>
  Population<F> NoiseReductionImpl(Population<F> &population) {
    Population<F> result;
    UncenteredTrickLPCA ulpca;
    result = ulpca.MVA(population);
    return result;
  }
*/
  // Debug purposes...
  /*
  template <typename F>
  Population<F> NoiseReductionImpl(Population<F> &population) {
    Population<F> result;
    SVDRepresentation svd;
    result = svd.MVA(population);
    return result;
  }
  */
};
}

#endif
