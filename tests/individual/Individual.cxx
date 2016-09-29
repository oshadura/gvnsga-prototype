#include "GATest.h"
#include "problem/DTLZ2.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"

#include <boost/interprocess/anonymous_shared_memory.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

class Individual : public GATest {
public:
};

using namespace boost::interprocess;

TEST_F(Individual, GeneratingIndividual) {
  geantvmoop::TGenes<geantvmoop::DTLZ2> i;
}

TEST_F(Individual, GenerationPopulation) {
  mapped_region region(anonymous_shared_memory(10000000000));
  geantvmoop::Population<geantvmoop::DTLZ2, 10> *pop = new (
      region.get_address()) geantvmoop::Population<geantvmoop::DTLZ2, 10>();
  //pop->InitSimplePopulation(*pop);
  pop->InitSharedMemPopulation(*pop);
  std::cout << (*pop);
}
