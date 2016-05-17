//===--- GAComparator.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GAComparator.h
 * @brief Implementation of comparison operator for LibGA prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __GACOMPAROR__
#define __GACOMPAROR__

#include "generic/TGenes.h"
#include <unordered_map>

template <typename F> class GAComparator {

public:
  std::unordered_map<geantvmoop::individual_t<F>, int> *fRank;
  std::unordered_map<geantvmoop::individual_t<F>, double> *fCrowDist;

  GAComparator(
      std::unordered_map<geantvmoop::individual_t<F>, int> *fIndRank,
      std::unordered_map<geantvmoop::individual_t<F>, double> *fIndCrowDist) {
    fRank = fIndRank;
    fCrowDist = fIndCrowDist;
  }

  bool operator()(geantvmoop::individual_t<F> a,
                  geantvmoop::individual_t<F> b) {
    if (fRank->find(a) == fRank->end() || fRank->find(b) == fRank->end())
      throw std::runtime_error("Please calculate the fRank indicator first!");
    if ((*fRank)[a] < (*fRank)[b])
      return true;
    else if ((*fRank)[a] > (*fRank)[b])
      return false;
    else {
      if (fCrowDist->find(a) == fCrowDist->end() ||
          fCrowDist->find(b) == fCrowDist->end())
        throw std::runtime_error(
            "Please calculate the crowding indicator first!");
      if ((*fCrowDist)[a] > (*fCrowDist)[b])
        return true;
      else if ((*fCrowDist)[a] < (*fCrowDist)[b])
        return false;
      else
        return std::tie((*fCrowDist)[a]) < std::tie((*fCrowDist)[b]);
    }
  }
};

#endif
