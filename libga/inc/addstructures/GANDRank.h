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

#ifndef MOO_FASTGANDRank_H
#define MOO_FASTGANDRank_H

#include "GAGenericIndicator.h"
#include <vector>

namespace geantvmoop {

class GANDRank : public GAGenericIndicator<GANDRank, int> {

public:
  template <typename T>
  static std::unordered_map<individual_t<T>, int>
  CalculateIndicator(Population<T> pop, int fBest = -1) {
    int fObserved = 0;
    if (fBest == -1)
      fBest = pop.size();
    std::unordered_map<individual_t<T>, int> fMap;
    for (auto entry : pop)
      fMap[entry] = std::numeric_limits<int>::max();
    std::vector<std::vector<int>> CurrentDominates(pop.size(),
                                                   std::vector<int>());
    std::vector<int> DominatedByOther(pop.size(), 0);
    std::vector<int> fFront;
    for (int i = 0; i < pop.size(); ++i) {
      for (int j = i + 1; j < pop.size(); ++j) {
        if (pop[i]->IsDominating(*pop[j])) {
          CurrentDominates[i].push_back(j);
          DominatedByOther[j] += 1;
        } else if (pop[i]->IsDominated(*pop[j])) {
          CurrentDominates[j].push_back(i);
          DominatedByOther[i] += 1;
        }
      }
      if (DominatedByOther[i] == 0) {
        fFront.push_back(i);
        fMap[pop[i]] = 0;
      }
    }
    fObserved += fFront.size();
    int AssignedRank = 1;
    while (!fFront.empty() && (fObserved < fBest)) {
      std::vector<int> fNextFront;
      for (int i : fFront) {
        for (int j : CurrentDominates[i]) {
          DominatedByOther[j] -= 1;
          if (DominatedByOther[j] == 0) {
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
};
}

#endif
