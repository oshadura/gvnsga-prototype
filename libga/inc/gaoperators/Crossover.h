#ifndef __CROSSOVER__
#define __CROSSOVER__

#include "generic/TGenes.h"

template <typename Derived> class Crossover {

public:
  template <typename F>
  static Genes<F> 
  CrossoverGA(Genes<F>  &p1, Genes<F>  &p2) {
    typename F::Input input =
        Derived::CrossoverGA(p1->GetInput(), p2->GetInput());
    return std::make_shared<Genes<F> >(input);
  }
};

#endif
