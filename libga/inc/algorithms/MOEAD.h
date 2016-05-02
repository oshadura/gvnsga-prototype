#ifndef __MOED__
#define __MOED__

#include "generic/Algorithm.h"
#include "generic/TGenes.h"
#include "generic/Population.h"
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

  PF<double> GetParetoFront() { return fFront; }
  /////////////////////////////////////////////////////////
  static void UpdateReferencePoint(std::vector<double> &ref,
                                   const Genes<double> &ind) {
    int fNObjectives = F::GetNObjectives();
    /*
    auto v = ind->getOutput();
    for (int i = 0; i < fNObjectives; ++i) {
      ref[i] = std::min(ref[i], v[i]);
    }
    */
  }

  static std::vector<double> GetReferencePoint(const Population<double> &pop) {
    int fNObjectives = F::GetNObjectives();
    std::vector<double> ref(fNObjectives);
    // WRONG !
    /*
    for (int i = 0; i < fNObjectives; ++i) {
      auto v = pop.GetFitness(i);
      ref[i] = *(std::min_element(v.begin(), v.end()));
    }
    */
    return ref;
  }

  // Output suppose to be fFitness
  double GetFitness(const Weights &w, const std::vector<double> &output) {
    return GetTchebichew(w, output);
  }

  // Output suppose to be fFitness
  double GetTchebichew(const Weights &w, const std::vector<double> &output) {
    double maxDistance = 0;
    for (int i = 0; i < T::getNumOfObjectives(); ++i) {
      maxDistance = std::max(maxDistance, w[i] * (output[i] - fRefPoint[i]));
    }
    return maxDistance;
  }

private:
  std::vector<Weights> fWeights;

  std::vector<std::vector<int> > fNearest;

  // Can we add this in a Reference point class?
  std::vector<double> fRefPoint;

  Population<double> pop;

  PF<double> fFront;

  TournamentSelection fSelection;

  int fGen;

  //================= User defined parameters ================//
  Double_t fPCross;
  Double_t fEtaCross;
  Int_t fNCross;
  Int_t fNMut;
  Int_t fNGen;
  Int_t fSizePop;
  Int_t fNParam;
  std::vector<std::pair<Double_t, Double_t> > fInterval;
  Int_t fNCons;
  Int_t fNObjectives;
  Double_t fPMut;
  Double_t fEtaMut;
  Double_t fEpsilonC;
  Bool_t fCrowdingObj;
  //================= User defined parameters ================//
};

#endif
