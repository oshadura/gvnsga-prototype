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

#ifndef MOO_MOEAD_H
#define MOO_MOEAD_H

#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "generic/TWeight.h"
#include "gaoperators/GASelection.h"
#include "tools/Random.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GAPolMutation.h"
#include "gaoperators/GASimpleSelection.h"

namespace geantvmoop {

template <typename F> class GAMOEAD : public GAAlgorithm<GAMOEAD<F>, F> {
private:
  std::vector<Weights> fWeights;
  std::vector<std::vector<int>> fNearest;
  std::vector<double> fRefPoint;
  Population<F> pop;
  std::vector<double> fFitness;
  PF<F> fFront;
  double PMut = 0.2;
  GASimpleSelection SimpleSelection;
  int fCounter = 0;

public:
  int fPopulationSize = 1000;
  int T = 5;

  GAMOEAD(F problem) : GAAlgorithm<GAMOEAD<F>, F>(problem) {}

  void InitializeImpl() {
    fWeights = Weights::GetWeights(fPopulationSize);
    fPopulationSize = fWeights.size();
    std::cout << fPopulationSize << std::endl;
    if (T >= fPopulationSize)
      throw std::runtime_error("Please set T lower than population size!");
    for (auto w : fWeights)
      fNearest.push_back(w.GetNearestNeighbor(fWeights, T));
    pop = Population<F>{fPopulationSize};
    fRefPoint = GetRP(pop);
    fFitness = std::vector<double>(fPopulationSize,
                                   std::numeric_limits<double>::max());
    for (int i = 0; i < pop.size(); ++i)
      fFitness[i] = GetFitness(fWeights[i], pop[i]->GetOutput());
  }

  template <typename T> int RandomVectorIndex(std::vector<T> v) {
    return Random::getInstance()->rndInt(0, v.size());
  }

  void EvolutionImpl() {
    fCounter = 0;
    for (int i = 0; i < pop.size(); ++i) {
      auto a = pop[fNearest[i][RandomVectorIndex(fNearest[i])]];
      auto b = pop[fNearest[i][RandomVectorIndex(fNearest[i])]];
      individual_t<F> offspring = GASBXCrossover::Crossover(a, b);
      if (Random::getInstance()->rndDouble() < PMut)
        offspring = GAPolMutation::Mutation(offspring);
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
      fFront.add(offspring);
    }
  }

  void PrintImpl(std::ostream &os) { os << fCounter << std::endl; }

  PF<F> GetParetoFrontImpl() { return fFront; }

  ////////////////////////////////////

  static void UpdateRP(std::vector<double> &ref, const individual_t<F> &ind) {
    int numOfObjectives = F::GetNObjectives();
    auto v = ind->GetOutput();
    for (int i = 0; i < numOfObjectives; ++i) {
      ref[i] = std::min(ref[i], v[i]);
    }
  }

  static std::vector<double> GetRP(const Population<F> &pop) {
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
