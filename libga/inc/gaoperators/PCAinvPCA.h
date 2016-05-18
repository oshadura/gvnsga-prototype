#ifndef __PCAINVPCA__
#define __PCAINVPCA__

#include "NoiseReduction.h"

namespace geantvmoop {

class PCAinvPCA : public NoiseReduction<PCAinvPCA> {

public:

  template <typename F>
  Population<F> NoiseReductionImpl(const Population<F> &population) {
    Population<F> result;


    return result;
  }

  template <typename F> void CovAnalisys(){}

  template <typename F> void PCA(){}

  template <typename F> void InversePCA(){}

  template <typename F> void ReEvaluation(){}

};
}

#endif
