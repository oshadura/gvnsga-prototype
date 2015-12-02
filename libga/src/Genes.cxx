#include <vector>
#include <ostream>
#include <string>
#include <utility>

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "Functions.h"
#include "Genes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"

ClassImp(Genes)

Genes::Genes()
    : TObject(), fFitness(0), fNObjectives(0), fRank(0), fDominationCounter(0),
      fCrowdingDistance(0), fEvaluated(0), fDominated(0), ConstViol(0),
      fGenes(0), fEpsilonC(0) {}

Genes::Genes(std::vector<Double_t> &f)
    : TObject(), fFitness(0), fNObjectives(f.size()), fRank(0),
      fDominationCounter(0), fCrowdingDistance(0), fEvaluated(0), fDominated(0),
      ConstViol(0), fGenes(f), fEpsilonC(0) {
  fFitness.reserve(fNObjectives);
}

Genes &Genes::operator=(const Genes &gen) {
  // comparison operator
  if (this != &gen) {
    fGenes = gen.fGenes;
    fRank = gen.fRank;
    fDominationCounter = gen.fDominationCounter;
    fEvaluated = gen.fEvaluated;
    fFitness = gen.fFitness;
  }
  return *this;
}

void Genes::Set(/*Double_t fAllev, Double_t fBuffev, Double_t fThread,
                Double_t fPriority, Double_t fSteps, Double_t fVector*/) {
  TRandom3 rand;
  std::vector<Double_t> *Genes = &fGenes;
  // Random numbers without limits per each parameter
  Int_t nparam = Functions::Instance()->GetNParam();
  for (Int_t i = 0; i < nparam; ++i) {
    fGenes[i] = rand.Uniform(Functions::Instance()->GetIntervalLimit(i).first,
                             Functions::Instance()->GetIntervalLimit(i).second);
  }
}

/*
void Genes::SetIt(Int_t i){
  TRandom3 rand;
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it) {
    for (const std::pair<Double_t, Double_t> &limit :
(Functions::Instance()->GetInterval())) {
      Double_t value = rand.Uniform(limit.first, limit.second);
      position = std::distance(GetIndividuals().begin(), it);
      SetGenes(position,Genes::SetGene(i,value));
    }
  }
}
*/
void Genes::Clear(Option_t * /*option*/) { TObject::Clear(); }

Double_t Genes::CheckDominance(const Genes *ind2) {
  if (ConstViol < 0 && ind2->ConstViol < 0) {
    if (ConstViol > ind2->ConstViol)
      return 1; // ind1 less
    else if (ConstViol < ind2->ConstViol)
      return -1; // ind2 violates less
    else
      return 0; // they both violate equally
  } else if (ConstViol < 0 && ind2->ConstViol == 0) {
    // ind1 violates and ind2 doesn't => ind1 dominates
    return -1;
  } else if (ConstViol == 0 && ind2->ConstViol < 0) {
    // ind1 doesn't violate and ind2 does => ind1 dominates
    return 1;
  } else {
    Int_t fFlag1, fFlag2;
    for (Int_t i = 0; i < ind2->GetNObjectives(); ++i) {
      if (ind2->GetNObjectives() > 1) {
        if (GetFitness(i) < ind2->GetFitness(i)) {
          fFlag1 = 1;
        } else if (GetFitness(i) > ind2->GetFitness(i)) {
          fFlag2 = 1;
        }
      } else {
        if (GetFitness(i) < ind2->GetFitness(i) &&
            fabs(GetFitness(i) - ind2->GetFitness(i)) > GetEpsilonC()) {
          fFlag1 = 1;
        } else if (GetFitness(i) > ind2->GetFitness(i) &&
                   fabs(GetFitness(i) - ind2->GetFitness(i)) > GetEpsilonC()) {
          fFlag2 = 1;
        }
      }
    }
    if (fFlag1 == 1 && fFlag2 == 0) {
      return 1;
    } else if (fFlag1 == 0 && fFlag2 == 1) {
      return -1;
    } else {
      return 0;
    }
  }
}

// Polynomial mutation
Double_t Genes::Mutate() {
  TRandom rand;
  Double_t fRrnd, fDelta1, fDelta2, fMutPow, fDelta, fValue;
  Double_t y, yl, yu, xy;
  Int_t fNMut = 0;
  for (Int_t j = 0; j < Functions::Instance()->GetNParam(); ++j) {
    if (rand.Rndm() <= AlgorithmNSGA::Instance()->GetPMut()) {
      y = fGenes[j];
      yl = Functions::Instance()->GetIntervalLimit(j).first;
      yu = Functions::Instance()->GetIntervalLimit(j).second;
      fDelta1 = (y - yl) / (yu - yl);
      fDelta2 = (yu - y) / (yu - yl);
      fRrnd = rand.Rndm();
      fMutPow = 1.0 / (AlgorithmNSGA::Instance()->GetEtaMut() + 1.0);
      if (fRrnd <= 0.5) {
        xy = 1.0 - fDelta1;
        fValue = 2.0 * fRrnd +
                 (1.0 - 2.0 * fRrnd) *
                     (pow(xy, (AlgorithmNSGA::Instance()->GetEtaMut() + 1.0)));
        fDelta = pow(fValue, fMutPow) - 1.0;
      } else {
        xy = 1.0 - fDelta2;
        fValue = 2.0 * (1.0 - fRrnd) +
                 2.0 * (fRrnd - 0.5) *
                     (pow(xy, (AlgorithmNSGA::Instance()->GetEtaMut() + 1.0)));
        fDelta = 1.0 - (pow(fValue, fMutPow));
      }
      y = y + fDelta * (yu - yl);
      if (y < yl)
        y = yl;
      if (y > yu)
        y = yu;
      fGenes[j] = y;
      fNMut += 1;
    }
  }
  return fNMut;
}

void Genes::WriteGenesTree(Genes &ind, Population &pop, const char *file) {
  if (!file) {
    TFile *f = new TFile(file, "RECREATE");
    TTree *tree = new TTree("gvga", "Genetic Algorithm TTree");
    tree->Branch("Population", &pop);
    f->Write();
  } else {
    TFile *f = TFile::Open(file, "RECREATE");
    TTree *tree = (TTree *)f->Get("gvga");
    TTree *tr_c = new TTree("gvga_1", "Genetic Algorithm friend Tree");
    tree->AddFriend("gvga_1", f);
    f->Write();
  }
}

void Genes::UpdateGenesTree(Genes &ind1, Genes &ind2, Population &pop,
                            const char *file) {
  // Looks it is not possible update existing events, lets update the tree
  TFile *f = TFile::Open(file, "RECREATE");
  if (!f) {
    return;
  }
  TTree *tree = (TTree *)f->Get("gvga");
  TTree *output = new TTree("gvga_1", "Changing individual in a population");
  output->Branch("Population", &pop);
  tree->AddFriend("gvga_1", f);
  f->Write();
}

void Genes::ReadGenesTree(Genes &ind, Population &pop, const char *file) {

  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("gvga");
  tree->SetBranchAddress("Population", &pop);
  Int_t entries = (Int_t)(tree->GetEntries());
}

void Genes::EvaluateGene() {}
