#ifndef __NSGA__
#define __NSGA__

#include "gacounters/GACD.h"
#include "gacounters/GAComparator.h"
#include "gacounters/GACounter.h"
#include "gacounters/GANDRank.h"

#include "gaoperators/SBXCrossover.h"
#include "gaoperators/PolynomialMutation.h"
#include "gaoperators/TournamentSelection.h"

#include "generic/Population.h"
#include "generic/Algorithm.h"

namespace geantvmoop {

template <typename F> class NSGA : public Algorithm<NSGA<F>, F> {

private:
  Population<F> fPopulation;
  std::unordered_map<individual_t<F>, double> fIndividualCrowDist;
  std::unordered_map<individual_t<F>, int> fIndividualRank;

public:
  NSGA(F Function) : Algorithm<NSGA<F>, F>(Function) {}
  int fPopulationSize = 100;
  double fMutation = 0.2;

  void Initialize() {
    fPopulation = Population<F>{fPopulationSize};
    fIndividualCrowDist = GACD::CalculateCDPop(fPopulation);
    fIndividualRank = GANDRank::CalculateRank(fPopulation);
  }

  void Evolution() {
    GAComparator<F> initialcomp(&fIndividualRank, &fIndividualCrowDist);
    TournamentSelection<GAComparator<F> > tournament(initialcomp);
    Population<F> fPopPool =
        tournament.SelectionGAMultiple(fPopulation, fPopulationSize * 2);
    for (unsigned int j = 0; j < fPopPool.size() - 1; j += 2) {
      individual_t<F> offspring =
          SBXCrossover::CrossoverGA(fPopPool[j], fPopPool[j + 1]);
      if (Generator::GetInstance().RNGDouble() < fMutation)
        PolynomialMutation::MutateGA(offspring);
      fPopulation.push_back(offspring);
    }
// Temporary until i will not find fix for Clang 3.7 (namespace) or upgrade to
// 3.8
//#ifdef __linux__
    fIndividualRank = GANDRank::CalculateRank(fPopulation);
//#endif
    fIndividualCrowDist = GACD::CalculateCDPop(fPopulation);
    GAComparator<F> nsgacomparator(&fIndividualRank, &fIndividualCrowDist);
    std::sort(fPopulation.begin(), fPopulation.end(), nsgacomparator);
    Population<F> fNextPop;
    for (int t = 0; t < fPopulationSize; ++t)
      fNextPop.push_back(fPopulation[t]);
    fPopulation = fNextPop;
  }

  void Print(std::ostream &os) {
    auto fLastIndividual = fPopulation[fPopulation.size() - 1];
    os << "pareto front: " << fIndividualRank[fLastIndividual]
       << " | worst crowding: ";
    os << fIndividualCrowDist[fLastIndividual] << std::endl;
  }

  PF<F> GetPF() {
    PF<F> fFront;
    for (unsigned int i = 0; i < fPopulation.size(); ++i)
      fFront.AddIndToPF(fPopulation[i]);
    return fFront;
  }
};
}

#endif
