#ifndef __ALGORITHMNSGA__
#define __ALGORITHMNSGA__

#include "TObject.h"
#include "Functions.h"
#include <vector>

template <class T> class Population;
template <class T> class Genes;
class Functions;

class AlgorithmNSGA {
public:
  AlgorithmNSGA();
  virtual ~AlgorithmNSGA();
  ////////////////////////////////////////////////////
  void Initialize();
  void Selection(Population<Double_t> &oldpop, Population<Double_t> &newpop);
  Genes<Double_t> &Tournament(Genes<Double_t> &ind1,
                              Genes<Double_t> &ind2) const;
  void Crossover(const Genes<Double_t> &parent1, const Genes<Double_t> &parent2,
                 Genes<Double_t> &child1, Genes<Double_t> &child2);
  void NextStep();
  void Evolution();
  ////////////////////////////////////////////////////
  Double_t GetPCross() const { return fPCross; }
  Double_t GetEtaCross() const { return fEtaCross; }
  Double_t GetPMut() const { return fPMut; }
  Double_t GetEtaMut() const { return fEtaMut; }
  //////////////////////////////////////////////////
  void SetPCross(Double_t pcross) { fPCross = pcross; }
  void SetEtaCross(Double_t etacross) { fEtaCross = etacross; }
  void SetEtaMut(Double_t etamut) { fEtaMut = etamut; }
  void SetPMut(Double_t pmut) { fPMut = pmut; }
  //////////////////////////////////////////////////
  Int_t GetGenTotalNumber() const { return fNGen; }
  void SetGenTotalNumber(Int_t gen) { fNGen = gen; }
  /////////////////////////////////////////////////
  static AlgorithmNSGA *Instance();

private:
  Double_t fPCross;
  Double_t fEtaCross;
  Double_t fPMut;
  Double_t fEtaMut;
  Int_t fNCross;
  Int_t fNMut;
  Int_t fNGen;
  static AlgorithmNSGA *fgNSGA2;

public:
  // Ugly instantiation
  Population<Double_t> *fParentPop;
  Population<Double_t> *fChildPop;
  Population<Double_t> *fMixedPop;

  ClassDef(AlgorithmNSGA, 1)
};

#endif
