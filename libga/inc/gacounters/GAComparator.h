#ifndef __GACOMPARATOR__
#define __GACOMPARATOR__

template <typename F> class GAComparator {

public:
  std::unordered_map<std::shared_ptr<Genes<F>>, int> *fRank;
  std::unordered_map<std::shared_ptr<Genes<F>>, double> *fCrowDist;

  GAComparator(
      std::unordered_map<std::shared_ptr<Genes<F>>, int> *fIndividualRank,
      std::unordered_map<std::shared_ptr<Genes<F>>, double> *fIndividualCrowDist) {
    fRank = fIndividualRank;
    fCrowDist = fIndividualCrowDist;
  }

  bool operator()(std::shared_ptr<Genes<F>> individual1,
                  std::shared_ptr<Genes<F>> individual2) {
    if (fRank->find(individual1) == fRank->end() ||
        fRank->find(individual2) == fRank->end())
      throw std::runtime_error(
          "Error in rank calculation. Calculate it again!");
    if ((*fRank)[individual1] < (*fRank)[individual2])
      return true;
    else if ((*fRank)[individual1] > (*fRank)[individual2])
      return false;
    else {
      if (fCrowDist->find(individual1) == fCrowDist->end() ||
          fCrowDist->find(individual2) == fCrowDist->end())
        throw std::runtime_error(
            "Error in crowding distance. Please recalculate it!");
      if ((*fCrowDist)[individual1] > (*fCrowDist)[individual2])
        return true;
      else if ((*fCrowDist)[individual1] < (*fCrowDist)[individual2])
        return false;
      else
        return std::tie((*fCrowDist)[individual1]) <
               std::tie((*fCrowDist)[individual2]);
    }
  }
};
#endif
