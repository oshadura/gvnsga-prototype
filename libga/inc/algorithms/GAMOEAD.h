//===--- GAMOEAD.h - LibGA -----------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GAMoead.h
 * @brief Implementation of MOEAD algorithms for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef MOO_MOEAD_H
#define MOO_MOEAD_H

#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "generic/TWeight.h"
#include "gaoperators/GASelection.h"
#include "tools/Random.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GAPolynomialMutation.h"
#include "gaoperators/GASimpleSelection.h"
#include "output/HistogramManager.h"

namespace geantvmoop {

template <typename F, std::size_t SizePop> class GAMOEAD : public GAAlgorithm<GAMOEAD<F,SizePop>, F> {
private:
  std::vector<Weights> fWeights;
  std::vector<std::vector<int>> fNearest;
  std::vector<double> fRefPoint;
  Population<F, SizePop> pop;
  std::vector<double> fFitness;
  PF<F, SizePop> fFront;
  double PMut = 0.2;
  GASimpleSelection SimpleSelection;
  int fCounter = 0;
  int fCurrentGeneration = 0;

public:
  int fPopulationSize = 10;
  // T closest  weight vectors
  int T = 7;

  GAMOEAD(F problem) : GAAlgorithm<GAMOEAD<F,SizePop>, F>(problem) {}

  void InitializeImpl() {
    fCurrentGeneration = 1; // initializing generation
    fWeights = Weights::GetWeights(fPopulationSize);
    fPopulationSize = fWeights.size();
    std::cout << "Population weights: " << fPopulationSize << std::endl;
    if (T >= fPopulationSize)
      throw std::runtime_error("Please set T lower than population size!");
    for (auto w : fWeights)
      fNearest.push_back(w.GetNearestNeighbor(fWeights, T));
    pop = Population<F, SizePop>();
    fRefPoint = GetRP(pop);
    fFitness = std::vector<double>(fPopulationSize,
                                   std::numeric_limits<double>::max());
    for (int i = 0; i < pop.size(); ++i)
      fFitness[i] = GetFitness(fWeights[i], pop[i]->GetOutput());
  }

  template <typename T> int RandomVectorIndex(std::vector<T> v) {
    return Random::GetInstance().RandomInt(0, v.size());
  }

  void EvolutionImpl() {
    fCounter = 0;
    HistogramManager<F, SizePop>::GetInstance().HistoFill(pop, "population_moead.root",
                                                 fCurrentGeneration);
    for (int i = 0; i < pop.size(); ++i) {
      auto a = pop[fNearest[i][RandomVectorIndex(fNearest[i])]];
      auto b = pop[fNearest[i][RandomVectorIndex(fNearest[i])]];
      individual_t<F> offspring = GASBXCrossover::Crossover(a, b);
      if (Random::GetInstance().RandomDouble() < PMut)
        offspring = GAPolynomialMutation::Mutation(offspring, PMut);
      UpdateRP(fRefPoint, offspring);
      auto out = offspring->GetOutput();
      for (int IterNearest : fNearest[i]) {
        double value = GetFitness(fWeights[IterNearest], out);
        if (value < fFitness[IterNearest]) {
          pop[IterNearest] = offspring;
          fFitness[IterNearest] = value;
          ++fCounter;
        }
      }
      fFront.Add(offspring);
    }
    // std::cout << fFront << std::endl;
    // std::cout << pop << std::endl;
    ++fCurrentGeneration;
  }

  void PrintImpl(std::ostream &os) {
    os << "Generation counter: " << fCounter << std::endl;
  }

  PF<F, SizePop> GetParetoFrontImpl() { return fFront; }

  static void UpdateRP(std::vector<double> &ref, const individual_t<F> &ind) {
    int numOfObjectives = F::GetNObjectives();
    auto v = ind->GetOutput();
    for (int i = 0; i < numOfObjectives; ++i) {
      ref[i] = std::min(ref[i], v[i]);
    }
  }

  static std::vector<double> GetRP(const Population<F, SizePop> &pop) {
    int numOfObjectives = F::GetNObjectives();
    std::vector<double> fRef(numOfObjectives);
    for (int i = 0; i < numOfObjectives; ++i) {
      auto v = pop.GetObjective(i);
      fRef[i] = *(std::min_element(v.begin(), v.end()));
    }
    return fRef;
  }

  double GetFitness(const Weights &w, const std::vector<double> &output) {
    return GetTchebichew(w, output);
  }

  double GetTchebichew(const Weights &w, const std::vector<double> &output) {
    double maxDistance = 0;
    for (int i = 0; i < F::GetNObjectives(); ++i) {
      maxDistance = std::max(maxDistance, w[i] * (output[i] - fRefPoint[i]));
    }
    return maxDistance;
  }
};
}

#endif
