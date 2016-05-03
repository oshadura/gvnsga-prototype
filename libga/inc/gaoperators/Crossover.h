#ifndef __CROSSOVER__
#define __CROSSOVER__

#include "generic/TGenes.h"

template <typename Derived> class Crossover {

public:
  template <typename T>
  static void CrossoverGA(Genes<T> &p1, Genes<T> &p2, Genes<T> &ch1,
                            Genes<T> &ch2) {
    Derived::CrossoverGA(p1, p2, ch1, ch2);
  }
};

#endif
