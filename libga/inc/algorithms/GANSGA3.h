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
#pragma once

#ifndef __NSGA3__
#define __NSGA3__

#include "addstructures/GACD.h"
#//include "addstructures/GAComparator.h"
//#include "addstructures/GANDRank.h"
#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "generic/ReferencePoint.h"
#include "tools/Random.h"
#include <iostream>

#include "addstructures/GAComparator.h"
#include "gaoperators/GAPolynomialMutation.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GATournamentSelection.h"

namespace geantvmoop {

template <typename F> class GANSGA3 : public GAAlgorithm<GANSGA3<F>, F> {

private:
  Population<F> population;
  std::vector<ReferencePoint> fReference;

public:
  GANSGA3(F problem) : GAAlgorithm<GANSGA3<F>, F>(problem) {}
  int fPopulationSize = 100;
  double PMut = 0.2;

  void InitializeImpl() {
    population = Population<F>{fPopulationSize};
    // more genertic..
    GenerateRP(&fReference, 3, 4);
  }

  void EvolutionImpl() {
    // Comparator based on reference points
    //GATournamentSelection<GAComparator<F> > selector(cmp);
    //Population<F> matingPool =
    //   selector.MultipleSelection(population, fPopulationSize * 2);
    //for (unsigned int j = 0; j < matingPool.size() - 1; j += 2) {
      //individual_t<F> offspring =
      //   GASBXCrossover::Crossover(matingPool[j], matingPool[j + 1]);
      //if (Random::GetInstance().RandomDouble() < PMut)
       // offspring = GAPolynomialMutation::Mutation(offspring, PMut);
      //population.push_back(offspring);
    //}
    //....
    Population<F> next;
    for (int l = 0; l < fPopulationSize; ++l)
      next.push_back(population[l]);
    population = next;

  }

  void PrintImpl(std::ostream &os) {
    for (int i = 0; i < population.size(); ++i) {
      std::cout << "Individual " << i << std::endl;
      for (int j = 0; j < population.GetTGenes(i).size(); ++j) {
        std::cout << population.GetGeneValue(i, j) << "|";
      }
      std::cout << "\nFitness function value: " << std::endl;
      for (int k = 0; k < population.GetTFitness(i).size(); ++k) {
        std::cout << population.GetObjectiveValue(i, k) << "|";
      }
    }
    std::cout << "---------------------------\n" << std::endl;
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
