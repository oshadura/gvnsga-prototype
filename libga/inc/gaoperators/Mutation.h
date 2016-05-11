#ifndef __MUTATION__
#define __MUTATION__

#include "generic/TGenes.h"

template <typename Derived> class Mutation {

public:
  template <typename F>
  static individual_t<F> MutationGA(const std::shared_ptr<Genes<F>> &ind1) {
    typename F::Input input = Derived::MutationGA(ind1->GetInput());
    return std::make_shared<Genes<F>>(input);
  }
};

#endif
