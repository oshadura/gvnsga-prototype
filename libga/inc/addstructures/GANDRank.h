//===--- GANDRank.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GANDRank.h
 * @brief Implementation of fast non-dominant ranking for algorithms for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __GANDRANK__
#define __GANDRANK__

#include "GAGenericIndicator.h"
#include <vector>

namespace geantvmoop {

class GANDRank : public GAGenericIndicator<GANDRank, int> {

public:
  /*
  template <typename F>
  static std::unordered_map<individual_t<F>, int>
  CalculateIndicatorImpl(Population<F> pop, int fBest = -1) {

    int fObserved = 0;

    if (fBest == -1)
      fBest = pop.size();

    std::unordered_map<individual_t<F>, int> fMap;

    for (auto entry : pop)
      fMap[entry] = std::numeric_limits<int>::max();

    std::vector<std::vector<int> > fMapIsDominating(pop.size(),
                                                    std::vector<int>());

    std::vector<int> fMapDominated(pop.size(), 0);

    std::vector<int> fFront;

    for (int i = 0; i < pop.size(); ++i) {
      for (int j = i + 1; j < pop.size(); ++j) {
        if (pop[i]->IsDominating(*pop[j])) {
          fMapIsDominating[i].push_back(j);
          fMapDominated[j] += 1;
        } else if (pop[i]->IsDominated(*pop[j])) {
          fMapIsDominating[j].push_back(i);
          fMapDominated[i] += 1;
        }
      }
      if (fMapDominated[i] == 0) {
        fFront.push_back(i);
        fMap[pop[i]] = 0;
      }
    }
    fObserved += fFront.size();
    int AssignedRank = 1;
    while (!fFront.empty() && (fObserved < fBest)) {
      std::vector<int> fNextFront;
      for (int i : fFront) {
        for (int j : fMapIsDominating[i]) {
          fMapDominated[j] -= 1;
          if (fMapDominated[j] == 0) {
            fNextFront.push_back(j);
            fMap[pop[j]] = AssignedRank;
          }
        }
      }
      fFront = fNextFront;
      fObserved += fFront.size();
      ++AssignedRank;
    }
    return fMap;
  }
  */

  template <typename F>
  static std::unordered_map<individual_t<F>, int>
  CalculateIndicatorImpl(Population<F> pop, int fBest = -1) {
    if (fBest == -1)
      fBest = pop.size();
    std::unordered_map<individual_t<F>, int> fMap;
    int AssignedRank = 0;
    int fObserved = 0;
    while (!pop.empty() && (fObserved < fBest)) {
      auto fFront = PF<F>::ParetoFrontND(pop);
      pop.Remove(fFront);
      for (individual_t<F> entry : fFront) {
        fMap[entry] = AssignedRank;
        ++fObserved;
      }
      ++AssignedRank;
    }
    for (individual_t<F> entry : pop)
      fMap[entry] = std::numeric_limits<int>::max();
    return fMap;
  }
};
}

#endif
