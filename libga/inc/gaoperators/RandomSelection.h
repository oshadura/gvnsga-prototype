#ifndef __RANDOMSELECTION__
#define __RANDOMSELECTION__

#include "gaoperators/Selection.h"
#include "gacounters/GAComparator.h"

namespace geantvmoop {

class RandomSelection : public Selection<RandomSelection> {

public:
  template <typename F>
  std::vector<std::shared_ptr<F> >
  SelectionGAUnary(const Population<F> &population) {
    return *(SelectRamdomly(population.begin(), population.end()));
  }

  template <typename F>
  Population<F> SelectionGAMultiple(const Population<F> &population, int n) {
    Population<F> fResultingPopulation;
    for (int i = 0; i < n; ++i) {
      fResultingPopulation.push_back(SelectionGAUnary(population));
    }
  }

private:
  template <typename Iterator, typename RandomGenerator>
  Iterator SelectRandomly(Iterator fStart, Iterator fEnd,
                          RandomGenerator &generator) {
    std::uniform_int_distribution<> distribution(0,
                                                 std::distance(fStart, fEnd) - 1);
    std::advance(fStart, distribution(generator));
    return fStart;
  }

  template <typename Iterator>
  Iterator SelectRandomly(Iterator fStart, Iterator fEnd) {
    static std::random_device fRandomDevice;
    static std::mt19937 generator(fRandomDevice());
    return select_randomly(fStart, fEnd, generator);
  }
};
} // namespace

#endif
