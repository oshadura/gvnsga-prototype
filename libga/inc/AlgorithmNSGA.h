#ifndef __ALGORITHMNSGA__
#define __ALGORITHMNSGA__

#include "TObject.h"
#include "Functions.h"
#include <vector>

template <class T> class Population;
template <class T> class Genes;
//class Functions; 

class AlgorithmNSGA {
public:
  AlgorithmNSGA();
  virtual ~AlgorithmNSGA();
  void Initialize() throw(ExceptionMessenger);
  void Selection(Population<Double_t> &oldpop, Population<Double_t> &newpop) throw(ExceptionMessenger);
  Genes<Double_t> &Tournament(Genes<Double_t> &ind1,
                              Genes<Double_t> &ind2) const;
  void Crossover(const Genes<Double_t> &parent1, const Genes<Double_t> &parent2,
                 Genes<Double_t> &child1, Genes<Double_t> &child2);
  void NextStep();
  void Evolution();
  ////////////////////////////////////////////////////
  Double_t GetPCross() const { return fPCross; }
  Double_t GetEtaCross() const { return fEtaCross; }
  //////////////////////////////////////////////////
  void SetPCross(Double_t pcross) { fPCross = pcross; }
  void SetEtaCross(Double_t etacross) { fEtaCross = etacross; }
  //////////////////////////////////////////////////
  Int_t GetGenTotalNumber() const { return fNGen; }
  void SetGenTotalNumber(Int_t gen) { fNGen = gen; }
  //====================================================//
  void SetFunction(Functions::functype f) { this->function = f; }
  void SetPopFunction(Functions::popfunctype f) { this->popfunction = f; }
  void SetCrowdingObj(Bool_t s) { this->fCrowdingObj = s; }
  void SetNParam(Int_t p) { this->fNParam = p; }
  void SetNCons(Double_t nc) { this->fNCons = nc; }
  void SetNObjectives(Double_t no) { this->fNObjectives = no; }
  void SetPopulationSize(Double_t ps) { this->fSizePop = ps; }
  void SetPMut(Double_t pmut) { this->fPMut = pmut; }
  void SetEtaMut(Double_t em) { this->fEtaMut = em; }
  void SetEpsilonC(Double_t ec) { this->fEpsilonC = ec; }
  void SetLimit(std::vector<std::pair<Double_t, Double_t>> lim) {
    this->fInterval = lim; 
  }
  void Report(std::ostream& os) const{
        os << "Population size = " << fSizePop
       << "\nNumber of generations = " << fNGen
       << "\nNumber of objective functions = " << fNObjectives
       << "\nNumber of constraints = " << fNCons
       << "\nNumber of variables = " << fNParam;

    if (fNParam != 0) {
        for (int i = 0; i < fNParam; ++i) {
            os << "\nLower limit of real variable " << (i+1)
               << " = " << fInterval[i].first;
            os << "\nUpper limit of real variable " << (i+1)
               << " = " << fInterval[i].second;
        }
        os << "\nProbability of crossover of real variable = " << fPCross;
        os << "\nProbability of mutation of real variable = " << fPMut;
        os << "\nDistribution index for crossover = " << fEtaCross;
        os << "\nDistribution index for mutation = " << fEtaMut;
    }
  }

private:
  Functions::functype function;
  Functions::popfunctype popfunction;
  Int_t fGen; // count
  //================= User defined parameters ================//
  Double_t fPCross;
  Double_t fEtaCross;
  Int_t fNCross;
  Int_t fNMut;
  Int_t fNGen;
  Int_t fSizePop;
  Int_t fNParam;
  std::vector<std::pair<Double_t, Double_t>> fInterval;
  Int_t fNCons;
  Int_t fNObjectives;
  Double_t fPMut;
  Double_t fEtaMut;
  Double_t fEpsilonC;
  Bool_t fCrowdingObj;
  //================= User defined parameters ================//

public:
  Population<Double_t> *fParentPop;
  Population<Double_t> *fChildPop;
  Population<Double_t> *fMixedPop;
  std::ofstream configuration;

  ClassDef(AlgorithmNSGA, 1)
};

#endif
