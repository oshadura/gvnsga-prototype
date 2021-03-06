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
#pragma once

#ifndef MOO_NSGAII_H
#define MOO_NSGAII_H

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
#include "addstructures/GAComparator.h"
#include "addstructures/GANDRank.h"
#include "addstructures/GACD.h"
#include "addstructures/GAComparator.h"
#include "tools/Random.h"
#include "output/CSVManager.h"
#include "output/HistogramManager.h"

#include <iostream>

namespace geantvmoop {

template <typename F> class GANSGA2 : public GAAlgorithm<GANSGA2<F>, F> {

private:
  Population<F> population;
  std::unordered_map<individual_t<F>, double> fIndCrowDist;
  std::unordered_map<individual_t<F>, int> fIndRank;

public:
  GANSGA2(F problem) : GAAlgorithm<GANSGA2<F>, F>(problem) {}
  int fPopulationSize = 10;
  double PMut = 0.4;
  double PCross = 0.9;
  int fCurrentGeneration = 0;

  void InitializeImpl() {
    fCurrentGeneration = 1; // initializing generation
    population = Population<F>{ fPopulationSize };
    fIndCrowDist = GACD::CalculateIndicator(population);
    fIndRank = GANDRank::CalculateIndicator(population);
    /*
    for (unsigned int i = 0; i < population.size(); ++i)
      std::cout << "| Pareto front: " << fIndRank[population[0]] << std::endl;
    */
  }

  void EvolutionImpl() {
    GAComparator<F> cmp(&fIndRank, &fIndCrowDist);
    GATournamentSelection<GAComparator<F> > selector(cmp);
    Population<F> matingPool =
        selector.MultipleSelection(population, fPopulationSize * 2);
    for (unsigned int j = 0; j < matingPool.size() - 1; j += 2) {
      // if (Random::GetInstance().RandomDouble() < PCross)
      individual_t<F> offspring = GASBXCrossover::Crossover(matingPool[j], matingPool[j + 1]);
      if (Random::GetInstance().RandomDouble() < PMut)
        offspring = GAPolynomialMutation::Mutation(offspring, PMut);
        //offspring = GASimpleMutation::Mutation(offspring, PMut);
      population.push_back(offspring);
    }
    fIndRank = GANDRank::CalculateIndicator(population);
    fIndCrowDist = GACD::CalculateIndicator(population);
    GAComparator<F> comp(&fIndRank, &fIndCrowDist);
    std::sort(population.begin(), population.end(), comp);
    //HistogramManager<F>::GetInstance().HistoFill(
    //    population, "population_nsga2.root", fCurrentGeneration);
    std::cout << "---------------------------\n" << std::endl;
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
    Population<F> next;
    for (int l = 0; l < fPopulationSize; ++l)
      next.push_back(population[l]);
    population = next;
    HistogramManager<F>::GetInstance().HistoFill(
        population, "population_nsga2.root", fCurrentGeneration);
    // std::cout << population << std::endl;
    CSVManager::GetInstance().CSVOutput("output.nsga", population, fIndRank, fIndCrowDist);
    ++fCurrentGeneration;
    std::cout << "Moving to next generation->" << fCurrentGeneration
              << std::endl;
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
