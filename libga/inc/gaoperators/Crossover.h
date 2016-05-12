#ifndef __CROSSOVER__
#define __CROSSOVER__

#include "generic/TGenes.h"

namespace geantvmoop {

template <typename Derived> class Crossover {

public:
  template <typename F>
  static individual_t<F> CrossoverGA(individual_t<F> &p1, individual_t<F> &p2) {
    typename F::Input input =
        Derived::CrossoverGA(p1->GetInput(), p2->GetInput());
    return std::make_shared<geantvmoop::Genes<F> >(input);
  }
};

} // end of namespace geantvmooop

#endif
