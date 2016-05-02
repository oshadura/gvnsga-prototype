#ifndef __CMAES__
#define __CMAES__

#include "generic/Algorithm.h"
#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/PF.h"

template <typename F> class CMAES : public Algorithm<CMAES<F>, F> {

public:
  CMAES(F problem) : Algorithm<CMAES<F>, F>(problem) {}

  void Initialize() {}

  void Evolution() {}

  void Print(std::ostream &os) { os << fGen << std::endl; }

  PF<double> GetParetoFront() { return fFront; }

private:
  Population<double> pop;

  PF<double> fFront;

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