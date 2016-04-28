#ifndef __PCA__
#define __PCA__

#include "generic/Population.h"

template <typename DerivedClass> class PCA {
public:

	template <typename T> Population<T> select(const Population<T> &population) {
            return static_cast<DerivedClass*>(this)->Select(population);
        }

};

#endif
