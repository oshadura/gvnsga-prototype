#ifndef __CMAES__
#define __CMAES__

#include "generic/GAAlgorithm.h"
#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/PF.h"

namespace geantvmoop {

template <typename F> class CMAES : public GAAlgorithm<CMAES<F>, F> {

public:
  CMAES(F problem) : GAAlgorithm<CMAES<F>, F>(problem) {}

  void Initialize() {}

  void Evolution() {}

  void Print(std::ostream &os) { os << fGen << std::endl; }

  PF<F> GetParetoFront() { return fFront; }

private:
  Population<F> pop;

  PF<F> fFront;

  int fGen;
};
}

#endif
