#ifndef __POPULATION__
#define __POPULATION__


#include "generic/TGenes.h"
#include <vector>
#include <memory>

template <typename F>
class Population : public std::vector<std::shared_ptr<Genes<F>> > {

};

#endif