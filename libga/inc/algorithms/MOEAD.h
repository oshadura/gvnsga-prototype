#ifndef __MOED__
#define __MOED__

#include "generic/Algorithm.h"
#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/PF.h"
#include "MOEADWeights.h"
#include "gaoperators/SBXCrossover.h"
#include "gaoperators/PolynomialMutation.h"
#include "gaoperators/TournamentSelection.h"


template <typename F>
class MOEAD : public Algorithm<MOEAD<F>, F> {

public:
  MOEAD(F problem) : Algorithm<MOEAD<F>, F>(problem) {}

  void Initialize() {}

  template <typename T> int RandomIndex(std::vector<double> v) {}

  void Evolution() {}

  void Print(std::ostream &os) { os << fGen << std::endl; }

  PF<F> GetParetoFront() { return fFront; }

  /*

  static void UpdateReferencePoint(std::vector<double> &fRef,
                                   const Genes<double> &ind) {
    int fNObjectives = F::GetNObjectives();
    auto fFitVector = ind->GetOutput();
    for (int i = 0; i < fNObjectives; ++i) {
      fRef[i] = std::min(ref[i], fFitVector[i]);
    }
  }

  static std::vector<double> GetReferencePoint(const Population<F> &pop) {
    int fNObjectives = F::GetNObjectives();
    std::vector<double> fRef(fNObjectives);
    for (int i = 0; i < fNObjectives; ++i) {
      auto fFitVector = pop.GetFitness(i);
      fRef[i] = *(std::min_element(fFitVector.begin(), fFitVector.end()));
    }
    return ref;
  }

  // Output suppose to be fFitness
  double GetFitness(const Weights &fWeight, const std::vector<double> &fFitVector) {
    return GetTchebichew(fWeight, fFitVector);
  }

  // Output suppose to be fFitness
  double GetTchebichew(const Weights &fWeight, const std::vector<double> &fFitVector) {
    double fMaxDistance = 0;
    for (int i = 0; i < F::GetNumOfObjectives(); ++i) {
      fMaxDistance = std::max(fMaxDistance, fWeight[i] * (fWeight[i] - fRefPoint[i]));
    }
    return fMaxDistance;
  }
  
  */
private:
  std::vector<Weights> fWeights;

  std::vector<std::vector<int> > fNearest;

  // Can we add this in a Reference point class?
  std::vector<double> fRefPoint;

  Population<F> pop;

  PF<F> fFront;

  int fGen;

  //RandomSelection selection;

};

#endif
