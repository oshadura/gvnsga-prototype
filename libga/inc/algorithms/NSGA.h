#ifndef __NSGA__
#define __NSGA__

#include "gacounters/GACD.h"
#include "gacounters/GAComparator.h"
#include "gacounters/GACounter.h"
#include "gacounters/GANDRank.h"

#include "gaoperators/SBXCrossover.h"
#include "gaoperators/PolynomialMutation.h"
#include "gaoperators/TournamentSelection.h"


template <typename F> class NSGA : public Algorithm<NSGA<F>, F> {

private:
  Population<F> population;
  std::unordered_map < std::shared_ptr<Genes<F>>, double> fIndividualCrowDist;
  std::unordered_map < std::shared_ptr<Genes<F>>, int> fIndividualRank;

public:
  NSGA(F Function) : Algorithm<NSGA<F>, F>(Function) {}

  int fPopulationSize = 100;
  double fMutation = 0.2;

  void Initialize() {
    fPopulation = Population<F>{ fPopulationSize };
    fIndividualCrowDist = GACD::CalculateCrowDist(fPopulation);
    fIndividualRank = GANDRank::CalculateRank(fPopulation);
  }

  void Evolution() {
    GAComparator<F> initialcomp(&fIndividualRank, &fIndividualCrowDist);
    TournamentSelection<GAComparator<F> > tournament(initialcomp);
    Population<F> fPopPool =
        tournament.SelectionGABetweenPops(fPopulation, fPopulationSize * 2);
    for (unsigned int j = 0; j < fPopPool.size() - 1; j += 2) {
      std::shared_ptr<Genes<F> offspring =
          SBXCrossover::CrossoverGA(fPopPool[j], fPopPool[j + 1]);
      if (Generator::fRandom = GetInstance()->rndDouble() < fMutation)
        offspring = PolynomialMutation::MutateGA(off);
      fPopulation.push_back(off);
    }
    fIndividualRank = NonDominatedRank::CalculateRank(fPopulation);
    fIndividualCrowDist = CrowdingDistance::CalculateCrowDist(fPopulation);
    GAComparator<F> nsgacomparator(&fIndividualRank, &fIndividualCrowDist);
    std::sort(fPopulation.begin(), fPopulation.end(), nsgacomparator);
    Population<F> fNextPop;
    for (int t = 0; t < fPopulationSize; ++t)
      fNextPop.push_back(fPopulation[t]);
    fPopulation = fNextPop;
  }

  void Print(std::ostream &os) {
    auto fLastIndividual = fPopulation[fPopulation.size() - 1];
    os << "pareto front: " << fIndividualRank[fLastIndividual] << " | worst crowding: ";
    os << fIndividualCrowDist[fLastIndividual] << std::endl;
  }

  PF<F> GetParetoFront() {
    ParetoFront<F> fFront;
    for (unsigned int i = 0; i < fPopulation.size(); ++i)
      fFront.add(fPopulation[i]);
    return fFront;
  }
};

#endif