//===--- GANSGA2UPCA.h - LibGA ----------------------------------------------*-
// C++
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
#pragma once

#ifndef MOO_NSGAIIUPCA_H
#define MOO_NSGAIIUPCA_H

#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

#define CLEAR "\033[2J" // clear screen escape code

#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "gaoperators/GATournamentSelection.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GAPolynomialMutation.h"
#include "gaoperators/GASimpleMutation.h"
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

template <typename F>
class GANSGA2UPCA : public GAAlgorithm<GANSGA2UPCA<F>, F> {

private:
  Population<F> population;
  std::unordered_map<individual_t<F>, double> fIndCrowDist;
  std::unordered_map<individual_t<F>, int> fIndRank;

public:
  GANSGA2UPCA(F problem) : GAAlgorithm<GANSGA2UPCA<F>, F>(problem) {}
  int fPopulationSize = 10;
  double PMut = 0.6;
  double PCross = 0.9;
  int fCurrentGeneration = 0;

  void InitializeImpl() {
    fCurrentGeneration = 1; // initializing generation
    population = Population<F>{fPopulationSize};
    fIndCrowDist = GACD::CalculateIndicator(population);
    fIndRank = GANDRank::CalculateIndicator(population);
  }

  void EvolutionImpl() {
    GAComparator<F> cmp(&fIndRank, &fIndCrowDist);
    GATournamentSelection<GAComparator<F>> selector(cmp);
    Population<F> matingPool =
        selector.MultipleSelection(population, fPopulationSize * 2);
    for (unsigned int j = 0; j < matingPool.size() - 1; j += 2) {
      individual_t<F> offspring =
          GASBXCrossover::Crossover(matingPool[j], matingPool[j + 1]);
      // std::cout << "Element for Crossover " << j << " and " << j + 1
      //          << std::endl;
      if (Random::GetInstance().RandomDouble() < PMut)
        offspring = GAPolynomialMutation::Mutation(offspring, PMut);
      // offspring = GASimpleMutation::Mutation(offspring, PMut);
      population.push_back(offspring);
      auto last = population.size() - 1;
      // std::cout << "Mutation had been happened with " << std::endl;
      // for (int j = 0; j < population.GetTGenes(last).size(); ++j) {
      //     std::cout << population.GetGeneValue(last, j) << "|";
      // }
      // std::cout << std::endl;
    }
    fIndRank = GANDRank::CalculateIndicator(population);
    fIndCrowDist = GACD::CalculateIndicator(population);
    GAComparator<F> comp(&fIndRank, &fIndCrowDist);
    std::sort(population.begin(), population.end(), comp);
    HistogramManager<F>::GetInstance().HistoFill(
        population, "population_nsga2_upcas.root", fCurrentGeneration);
    Population<F> next;
    for (int l = 0; l < fPopulationSize; ++l)
      next.push_back(population[l]);
    std::cout << "--------------TRANFORMATION IS GOING-------------\n"
              << std::endl;
    if (fCurrentGeneration > 10) {
      PCAinvPCA cleanupoperator;
      population = cleanupoperator.NR(next);
      // Modification to avoid 0 equal ranks and crowding distance
      fIndRank = GANDRank::CalculateIndicator(population);
      fIndCrowDist = GACD::CalculateIndicator(population);
      GAComparator<F> comp(&fIndRank, &fIndCrowDist);
      std::sort(population.begin(), population.end(), comp);
    } else {
      population = next;
    }
    std::cout << "-----------------------------------------------\n"
              << std::endl;
    std::cout << "---------------After transformation------------\n"
              << std::endl;
    for (int i = 0; i < population.size(); ++i) {
      std::cout << "Individual " << i << std::endl;
      for (int j = 0; j < population.GetTGenes(i).size(); ++j) {
        std::cout << GREEN << population.GetGeneValue(i, j) << "|";
      }
      std::cout << RESET << "\nFitness function value: " << std::endl;
      for (int k = 0; k < population.GetTFitness(i).size(); ++k) {
        std::cout << BLUE << population.GetObjectiveValue(i, k) << "|";
      }
      auto ind = population[i];
      std::cout << RESET << "\n| Rank: " << RED << fIndRank[ind] << RESET
                << " | Crowding distance value: ";
      std::cout << MAGENTA << fIndCrowDist[ind] << RESET << std::endl;
    }
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "Moving to next generation " << fCurrentGeneration
              << std::endl;
    CSVManager::GetInstance().CSVOutput("output.nsgaupca", population, fIndRank,
                                        fIndCrowDist);
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
