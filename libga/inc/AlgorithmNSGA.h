#ifndef __ALGORITHMNSGA__
#define __ALGORITHMNSGA__

#include "TObject.h"
#include "Functions.h"
#include <vector>

class Population;
class Genes;
class Functions;

class AlgorithmNSGA {
public:
  AlgorithmNSGA();
  virtual ~AlgorithmNSGA();
  ////////////////////////////////////////////////////
  void Initialize();
  void Selection(Population &oldpop, Population &newpop);
  Genes &Tournament(Genes &ind1, Genes &ind2) const;
  void Crossover(const Genes &parent1, const Genes &parent2, Genes &child1,
                 Genes &child2);
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
  Int_t GetGenTotalNumber() const {return fNGen;}
  /////////////////////////////////////////////////
  static AlgorithmNSGA *Instance();

private:
  Population *fParentPop;
  Population *fChildPop;
  Population *fMixedPop;
  Double_t fPCross;
  Double_t fEtaCross;
  Double_t fPMut;
  Double_t fEtaMut;
  Int_t fNCross;
  Int_t fNMut;
  Int_t fNGen;

  static AlgorithmNSGA* fgNSGA2;

  ClassDef(AlgorithmNSGA, 1)
};

#endif
