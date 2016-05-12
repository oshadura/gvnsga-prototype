#ifndef __PF__
#define __PF__

#include <list>
#include <generic/Population.h>

namespace geantvmoop{

template <typename F> class PF {

  std::list<individual_t<F>> fFront;

public:
  Population<F> GetPopulation() const {
    Population<F> fResultPop;
    for (auto ind : fFront)
      fResultPop.push_back(ind);
    return fResultPop;
  }

  bool AddIndToPF(const individual_t<F> &ind) {
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
    std::list<individual_t<F>> fFront;
    if (pop.empty())
      return Population<F>();
    auto fFunction = [&fFront](const individual_t<F> &ind) {
      for (auto it = fFront.begin(); it != fFront.end();) {
        if ((*it)->IsDominating(*ind))
          return false;
        if ((*it)->IsDominated(*ind))
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

}

#endif
