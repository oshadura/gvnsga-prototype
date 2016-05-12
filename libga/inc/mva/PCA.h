#ifndef __PCA__
#define __PCA__

#include "generic/Population.h"

namespace geantvmoop{

template <typename DerivedClass> class PCA {
public:

	template <typename F> Population<F> MVA(const Population<F> &population) {
            return static_cast<DerivedClass*>(this)->MVA(population);
        }

};

}

#endif
