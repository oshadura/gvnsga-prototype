#ifndef __GANDRANK__
#define __GANDRANK__

#include <vector>

#include "generic/PF.h"
#include "GACounter.h"

class GANDRank : public GACounter<GANDRank, double> {
public:
  GANDRank();

  ~GANDRank();

  GANDRank(GANDRank const &) = delete;

  void operator=(GANDRank const &) = delete;

  template <typename F>
  static std::unordered_map<std::vector<std::shared_ptr<Genes<F> > >, int>
  CalculateRank(Population<F> pop, int fBestIndividuals = -1) {
    if (fBestIndividuals == -1)
      fBestIndividuals = pop.size();
    std::unordered_map<std::vector<std::shared_ptr<Genes<F> > >, int> fMap;
    int fRank = 0;
    int fNumberOfIndividuals = 0;
    while (!pop.empty() && (fNumberOfIndividuals < fBestIndividuals)) {
      auto fFront = PF<F>::GetPF(pop);
      pop.Remove(fFront);
      for (std::vector<std::shared_ptr<Genes<F> > > entry : fFront) {
        fMap[entry] = fRank;
        ++fNumberOfIndividuals;
      }
      ++fRank;
    }
    for (std::vector<std::shared_ptr<Genes<F> > > entry : pop)
      fMap[entry] = std::numeric_limits<int>::max();
    return fMap;
  }
};

#endif