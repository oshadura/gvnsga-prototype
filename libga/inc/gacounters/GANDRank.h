#ifndef __GANDRANK__
#define __GANDRANK__

#include <vector>

#include "generic/PF.h"
#include "GACounter.h"

namespace geantvmoop {

class GANDRank : public GACounter<GANDRank, int> {

  template <typename F> using Map = std::unordered_map<individual_t<F>, int>;

public:
  GANDRank() {}

  ~GANDRank() {}

  GANDRank(GANDRank const &) = delete;

  void operator=(GANDRank const &) = delete;

  template <typename F>
  static std::unordered_map<std::vector<individual_t<F> >, int>
  CalculateRank(Population<F> pop, int fBestIndividuals = -1) {
    if (fBestIndividuals == -1)
      fBestIndividuals = pop.size();
    Map<F> fMap;
    int fRank = 0;
    int fNumberOfIndividuals = 0;
    while (!pop.empty() && (fNumberOfIndividuals < fBestIndividuals)) {
      auto fFront = PF<F>::GetPF(pop);
      pop.Remove(fFront);
      for (individual_t<F> entry : fFront) {
        fMap[entry] = fRank;
        ++fNumberOfIndividuals;
      }
      ++fRank;
    }
    for (individual_t<F> entry : pop)
      fMap[entry] = std::numeric_limits<int>::max();
    return fMap;
  }
};

} // end of namespace geantvmooop

#endif
