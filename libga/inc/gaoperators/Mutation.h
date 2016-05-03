#ifndef __MUTATION__
#define __MUTATION__

#include "generic/TGenes.h"

template <typename Derived> class Mutation {

public:
  template <typename T>
  void MutationGA(const Genes<T> &ind1, const Genes<T> &ind2) {
    Derived::MutationGA(ind1, ind2);
  }
};

#endif