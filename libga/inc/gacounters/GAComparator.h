#ifndef __GACOMPARATOR__
#define __GACOMPARATOR__

namespace geantvmoop{

template <typename F> class GAComparator {

public:
  std::unordered_map<individual_t<F>, int> *fRank;
  std::unordered_map<individual_t<F>, double> *fCrowDist;

  GAComparator(
      std::unordered_map<individual_t<F>, int> *fIndividualRank,
      std::unordered_map<individual_t<F>, double> *fIndividualCrowDist) {
    fRank = fIndividualRank;
    fCrowDist = fIndividualCrowDist;
  }

  ~GAComparator() {}

  GAComparator(GAComparator const &) = default;

  void operator=(GAComparator const &);

  bool operator()(individual_t<F> individual1, individual_t<F> individual2) {
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

} // end of namespace geantvmooop

#endif
