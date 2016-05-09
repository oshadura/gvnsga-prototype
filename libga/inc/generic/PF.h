#ifndef __PF__
#define __PF__

#include <list>
#include <generic/Population.h>

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
    for (auto it = fFront.begin(); it != fFront.end();) {
      if ((*it)->IsDominating(*ind) || (*it)->IsEqual(*ind))
        return false;
      if ((*it)->IsDominated(*ind))
        fFront.erase(it++);
      else
        ++it;
    }
    fFront.push_back(ind);
    return true;
  }

  static Population<F> GetPF(const Population<F> &pop) {
    std::list<std::shared_ptr<Genes<F>>> fFront;
    if (pop.empty())
      return Population<F>();
    auto fFunction = [&fFront](const std::shared_ptr<Genes<F> > &ind) {
      for (auto it = fFront.begin(); it != fFront.end();) {
        if ((*it)->IsDominating(*ind))
          return false;
        if ((*it)->IsDominatedBy(*ind))
          fFront.erase(it++);
        else
          ++it;
      }
      fFront.push_back(ind);
      return true;
    };
    for (unsigned int i = 0; i < pop.size(); ++i)
      fFunction(pop[i]);
    Population<F> fResult;
    for (auto ind : fFront)
      fResult.push_back(ind);
    return fResult;
  }
};

#endif