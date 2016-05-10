#ifndef __TOURNAMENTSELECTION__
#define __TOURNAMENTSELECTION__

#include "gaoperators/Selection.h"
#include "gacounters/GAComparator.h"

template <typename GAComparator>
class TournamentSelection: public Selection<TournamentSelection<GAComparator>> {

private:
  GAComparator comparator;

public:
  TournamentSelection(const GAComparator & comparator)
      : comparator(comparator) {}

  template <typename F>
  std::vector<std::shared_ptr<F> > SelectionGAUnary(const Population<F> &
                                                  population) {
    throw std::runtime_error("Tournament Selection does not allow to "
                             "select only single individuals!");
  }

  template <typename F>
  Population<F> SelectionGAMultiple(const Population<F> & population, int n) {
    Population<F> fResult;
    Population<F> fPool;
    while (fResult.size() < n) {
      auto fIndex = population.GetIndex();
      std::random_shuffle(fIndex.begin(), fIndex.end());
      for (unsigned int i = 0; i < fIndex.size() - 1; i += 2) {
        fResult.push_back(
            std::min(population[fIndex[i]], population[fIndex[i + 1]]));
        if (fResult.size() >= n)
          break;
      }
    }
    return fResult;
  }
};

#endif
