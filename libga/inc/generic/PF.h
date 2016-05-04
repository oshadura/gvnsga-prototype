#ifndef __PF__
#define __PF__

#include <list>

template <typename F> class PF {

  std::list<std::shared_ptr<Genes<F> > > fFront;

public:
  Population<F> GetPopulation() const {
    Population<F> fResultPop;
    for (auto ind : fFront)
      fResultPop.push_back(ind);
    return fResultPop;
  }

  bool AddIndToPF(const std::shared_ptr<Genes<F> > &ind) {
    // for every element of front
    for (auto it = fFront.begin(); it != fFront.end();) {
      // of one elements dominates ind -> does not belong to front
      if ((*it)->isDominating(*ind) || (*it)->isEqual(*ind))
        return false;
      // else remove all elements that are dominated by ind
      if ((*it)->isDominatedBy(*ind))
        fFront.erase(it++);
      else
        ++it;
    }
    fFront.push_back(ind);
    return true;
  }
};

#endif