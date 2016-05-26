#ifndef __PCAINVPCA__
#define __PCAINVPCA__

#include "NoiseReduction.h"

namespace geantvmoop {

class PCAinvPCA : public NoiseReduction<PCAinvPCA> {

public:
  template <typename F>
  Population<F> NoiseReductionImpl(const Population<F> &population) {
    Population<F> result;
    CovAnalisys();
    // if enough variance level go ahead..
    PCA();
    InversePCA();
    return result;
  }

  template <typename F> void CovAnalisys() {}

  template <typename F> void PCA() {
    // Distribution is normal then -> MVA::LPCA() else -> MVA::KPCA()
  
  }

  template <typename F> void InversePCA() {}
};
}

#endif
