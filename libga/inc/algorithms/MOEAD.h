#ifndef __MOED__
#define __MOED__

#include "generic/Algorithm.h"
#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/PF.h"
#include "generic/GAWeights.h"
#include "gaoperators/SBXCrossover.h"
#include "gaoperators/PolynomialMutation.h"
#include "gaoperators/RandomSelection.h"

namespace geantvmoop {

template <typename F> class MOEAD : public Algorithm<MOEAD<F>, F> {

public:
  MOEAD(F problem) : Algorithm<MOEAD<F>, F>(problem) {}

  void Initialize() {
    fWeights = Weights::Weights(fPopulationSize);
    fPopulationSize = fWeights.size();
    for (auto weight : fWeights) {
      fNearest.push_back(weight.GetNearestNeighborByIndex(fWeights, T));
    }
    fPopulation = Population<F>{ fPopulationSize };
    fRefPoint = GetRP(fPopulationSize);
    fFitness = std::vector<double>(fPopulationSize,
                                   std::numeric_limits<double>::max());
    for (int i = 0; i < fPopulation.size(); ++i) {
      fFitness[i] = GetFitness(fWeights[i], fPopulation[i]->GetOutput());
    }
  }

  void Evolution() {
    for (int i = 0; i < fPopulation.size(); ++i) {
      auto fIndividual1 = fPopulation[fNearest[i][RandomIndex(fNearest[i])]];
      auto fIndividual2 = fPopulation[fNearest[i][RandomIndex(fNearest[i])]];
      individual_t<F> fNewIndividual =
          SBXCrossover::CrossoverGA(fIndividual1, fIndividual2);
      UpdateRP(fRefPoint, fNewIndividual);
      auto fNewFitness = fNewIndividual->GetOutput();
      for (auto &&itnearest : fNearest[i]) {
        double fValue = GetFitness(fWeights[itnearest], fNewFitness);
        if (fValue < fFitness[itnearest]) {
          fPopulation[itnearest] = fNewIndividual;
          fFitness[itnearest] = fValue;
        }
      }
      fFront.AddIndToPF(fNewIndividual);
    }
  }

  void Print(std::ostream &os) { ; }

  PF<F> GetPF() { return fFront; }

  template <typename T> int RandomIndex(std::vector<T> fVectorIndex) {
    return Generator::GetInstance().RNGInteger(0, fVectorIndex.size());
  }

  static void UpdateRP(std::vector<double> &fRef, const individual_t<F> &ind) {
    int fNObjectives = F::GetNObjectives();
    auto fFitVector = ind->GetOutput();
    for (int i = 0; i < fNObjectives; ++i) {
      fRef[i] = std::min(fRef[i], fFitVector[i]);
    }
  }

  static std::vector<double> GetRP(const Population<F> &pop) {
    int fNObjectives = F::GetNObjectives();
    std::vector<double> fRef(fNObjectives);
    for (int i = 0; i < fNObjectives; ++i) {
      auto fFitVector = pop.GetObjective(i);
      fRef[i] = *(std::min_element(fFitVector.begin(), fFitVector.end()));
    }
    return fRef;
  }

  // Output suppose to be fFitness
  double GetFitness(const Weights &fWeight,
                    const std::vector<double> &fFitVector) {
    return GetTchebichew(fWeight, fFitVector);
  }

  // Output suppose to be fFitness
  double GetTchebichew(const Weights &fWeight,
                       const std::vector<double> &fFitVector) {
    double fMaxDistance = 0;
    for (int i = 0; i < F::GetNObjectives(); ++i) {
      fMaxDistance =
          std::max(fMaxDistance, fWeight[i] * (fWeight[i] - fRefPoint[i]));
    }
    return fMaxDistance;
  }

private:
  std::vector<Weights> fWeights;
  std::vector<std::vector<int> > fNearest;
  // TBD: Can we add this in a Reference point class?
  std::vector<double> fRefPoint;
  Population<F> fPopulation;
  std::vector<double> fFitness;
  PF<F> fFront;
  RandomSelection selection;
  double fMutProbability = 0.3;

public:
  int fPopulationSize = 100;
  int T = 5;
};
}

#endif
