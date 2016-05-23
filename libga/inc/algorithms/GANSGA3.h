//===--- GANSGA3.h - LibGA ----------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Nsga2.h
 * @brief Implementation of NSGA-III algorithms for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __NSGA3__
#define __NSGA3__

#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "generic/ReferencePoint.h"
#include "addstructures/GANDRank.h"
#include "addstructures/GACD.h"
#include "addstructures/GAComparator.h"
#include "tools/Random.h"
#include <iostream>

#include "gaoperators/GATournamentSelection.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GAPolMutation.h"
#include "addstructures/GAComparator.h"

namespace geantvmoop {

template <typename F> class GANSGA3 : public GAAlgorithm<GANSGA3<F>, F> {

private:
  Population<F> population;
  //std::unordered_map<individual_t<F>, double> fIndCrowDist;
  //std::unordered_map<individual_t<F>, int> fIndRank;
  std::vector<ReferencePoint> fReference;

public:
  GANSGA3(F problem) : GAAlgorithm<GANSGA3<F>, F>(problem) {}
  int fPopulationSize = 1000;
  double PMut = 0.2;

  void InitializeImpl() {
    population = Population<F>{fPopulationSize};
    //Stupid solution only for DTLZ1!
    GenerateRP(&fReference, 3, 4);
    //fIndCrowDist = GACD::CalculateIndicator(population);
    //fIndRank = GANDRank::CalculateIndicator(population);
  }

  void EvolutionImpl() {
    /*
    GAComparator<F> cmp(&fIndRank, &fIndCrowDist);
    GATournamentSelection<GAComparator<F>> selector(cmp);
    Population<F> matingPool =
        selector.MultipleSelection(population, fPopulationSize * 2);
    for (unsigned int j = 0; j < matingPool.size() - 1; j += 2) {
      individual_t<F> offspring =
          GASBXCrossover::Crossover(matingPool[j], matingPool[j + 1]);
      if (Random::getInstance()->rndDouble() < PMut)
        offspring = GAPolMutation::Mutation(offspring);
      population.push_back(offspring);
    }
    fIndRank = GANDRank::CalculateIndicator(population);
    fIndCrowDist = GACD::CalculateIndicator(population);
    GAComparator<F> comp(&fIndRank, &fIndCrowDist);
    std::sort(population.begin(), population.end(), comp);
    Population<F> next;
    for (int l = 0; l < fPopulationSize; ++l)
      next.push_back(population[l]);
    population = next;
    */
  }

  void PrintImpl(std::ostream &os) {
    /*
    auto last = population[population.size() - 1];
    os << "Pareto front: " << fIndRank[last] << " | worst crowding: ";
    os << fIndCrowDist[last] << std::endl;
    */
  }

  PF<F> GetParetoFrontImpl() {
    PF<F> fFront;
    for (unsigned int i = 0; i < population.size(); ++i)
      fFront.Add(population[i]);
    return fFront;
  }

private:
  Population<F> pop;

  PF<F> fFront;
};
}

#endif
