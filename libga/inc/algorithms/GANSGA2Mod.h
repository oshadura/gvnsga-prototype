//===--- GANSGA2.h - LibGA ----------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Nsga2.h
 * @brief Implementation of NSGA-II algorithms for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef MOO_NSGAIIMOD_H
#define MOO_NSGAIIMOD_H

#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "gaoperators/GATournamentSelection.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GAPolMutation.h"
#include "gaoperators/PCAinvPCA.h"
#include "addstructures/GAComparator.h"
#include "addstructures/GANDRank.h"
#include "addstructures/GACD.h"
#include "addstructures/GAComparator.h"
#include "tools/Random.h"
#include "output/CSVManager.h"
#include "output/HistogramManager.h"

#include <iostream>

namespace geantvmoop {

template <typename F> class GANSGA2Mod : public GAAlgorithm<GANSGA2Mod<F>, F> {

private:
  Population<F> population;
  std::unordered_map<individual_t<F>, double> fIndCrowDist;
  std::unordered_map<individual_t<F>, int> fIndRank;

public:
  GANSGA2Mod(F problem) : GAAlgorithm<GANSGA2Mod<F>, F>(problem) {}
  int fPopulationSize = 10;
  double PMut = 0.2;
  int fCurrentGeneration = 0;

  void InitializeImpl() {
    fCurrentGeneration = 1; // initializing generation
    population = Population<F>{ fPopulationSize };
    fIndCrowDist = GACD::CalculateIndicator(population);
    fIndRank = GANDRank::CalculateIndicator(population);
  }

  void EvolutionImpl() {
    GAComparator<F> cmp(&fIndRank, &fIndCrowDist);
    GATournamentSelection<GAComparator<F> > selector(cmp);
    Population<F> matingPool =
        selector.MultipleSelection(population, fPopulationSize * 2);
    for (unsigned int j = 0; j < matingPool.size() - 1; j += 2) {
      individual_t<F> offspring =
          GASBXCrossover::Crossover(matingPool[j], matingPool[j + 1]);
      if (Random::GetInstance().RandomDouble() < PMut)
        offspring = GAPolMutation::Mutation(offspring);
      population.push_back(offspring);
    }
    fIndRank = GANDRank::CalculateIndicator(population);
    fIndCrowDist = GACD::CalculateIndicator(population);
    GAComparator<F> comp(&fIndRank, &fIndCrowDist);
    std::sort(population.begin(), population.end(), comp);
    HistogramManager<F>::GetInstance().HistoFill(
        population, "population_nsga2_mod.root", fCurrentGeneration);
    std::cout << "---------------------------\n" << std::endl;
    for (auto entry : population) {
      //  for (int i = 0; i < entry.GetOutput().size(); ++i) {
      std::cout << entry->GetOutput()[0] << ", " << entry->GetOutput()[1]
                << ", " << entry->GetOutput()[3] << ", ";
      //  }
      std::cout << " | rank: " << fIndRank[entry] << " | crowding: ";
      std::cout << fIndCrowDist[entry] << std::endl;
    }
    std::cout << "---------------------------\n" << std::endl;
    Population<F> next;
    for (int l = 0; l < fPopulationSize; ++l)
      next.push_back(population[l]);
    std::cout << "--------------TRANFORMATION-------------\n" << std::endl;
    if (fCurrentGeneration % 2 == 0) {
      PCAinvPCA cleanupoperator;
      population = cleanupoperator.NR(next);
    } else {
      population = next;
    }
    std::cout << "Moving to next generation " << fCurrentGeneration
              << std::endl;
    CSVManager::GetInstance().CSVOutput("output.nsgalpca", population);
    ++fCurrentGeneration;
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "---------------------------\n" << std::endl;
  }

  void PrintImpl(std::ostream &os) {
    // Print everything!
    auto last = population[population.size() - 1];
    os << "| Pareto front: " << fIndRank[last] << " | Crowding distance: ";
    os << fIndCrowDist[last] << std::endl;
  }

  PF<F> GetParetoFrontImpl() {
    PF<F> fFront;
    for (unsigned int i = 0; i < population.size(); ++i)
      fFront.Add(population[i]);
    std::cout << fFront << std::endl;
    return fFront;
  }
};
}

#endif