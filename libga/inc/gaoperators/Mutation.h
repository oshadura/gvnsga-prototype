#ifndef __MUTATION__
#define __MUTATION__

#include "generic/TGenes.h"

namespace geantvmoop {

template <typename Derived> class Mutation {

public:
  template <typename F>
  static individual_t<F> MutationGA(individual_t<F> &ind1) {
    typename F::Input input = ind1->GetInput();
    Derived::MutationGA(input);
    return std::make_shared<geantvmoop::Genes<F> >(input);
  }
};

} // end of namespace geantvmooop

#endif
