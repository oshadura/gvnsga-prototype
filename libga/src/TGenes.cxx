#include <vector>
#include <ostream>
#include <string>
#include <utility>
#include <random>
#include <algorithm>
#include <stdexcept>

// Map for Genes[x]<->Limits[x] (?)
#include <map>

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "Functions.h"
#include "TGenes.h"
#include "Population.h"
#include "HistogramManager.h"

// ClassImp(Genes<T>)

template <class T>
Genes<T>::Genes() throw()
    : TObject(), fFitness(0), fRank(0), fDominationCounter(0),
      fCrowdingDistance(0), fEvaluated(false), fDominated(), ConstViol(0),
      fGenes(), fAllev(0), fBuffev(0), fThread(0), fPriority(0), fSteps(0),
      fVector(0), fTime(0), fMemory(0), setup(0), fConstraines(0) {}

template <class T>
Genes<T>::Genes(const Functions &config) throw(ExceptionMessenger)
    : TObject(), fFitness(0), fRank(0), fDominationCounter(0),
      fCrowdingDistance(0), fEvaluated(false), fDominated(), ConstViol(0),
      fGenes(), fAllev(0), fBuffev(0), fThread(0), fPriority(0), fSteps(0),
      fVector(0), fTime(0), fMemory(0), setup(&config), fConstraines(0) {
  fGenes.resize(setup->fNParam, 0);
  fFitness.resize(setup->fNObjectives, 0);
  fConstraines.resize(setup->fNCons, 0);
}

template <class T> Genes<T>::Genes(Genes &f) {}

template <class T> Genes<T> &Genes<T>::operator=(const Genes<T> &gen) {
  if (this != &gen) {
    fGenes = gen.fGenes;
    fRank = gen.fRank;
    fDominationCounter = gen.fDominationCounter;
    fEvaluated = gen.fEvaluated;
    fFitness = gen.fFitness;
    fCrowdingDistance = gen.fCrowdingDistance;
    fDominated = gen.fDominated;
    ConstViol = gen.ConstViol;
  }
  return *this;
}

template <class T> void Genes<T>::Set() throw(ExceptionMessenger) {
  /*
  if (!setup)
    throw ExceptionMessenger("Do something with setup function!");
  */
  TRandom3 rand;
  rand.SetSeed(time(NULL));
  for (Int_t i = 0; i < (setup->fNParam); ++i) {
    fGenes[i] = rand.Uniform(setup->GetIntervalLimit(i).first,
                             setup->GetIntervalLimit(i).second);
  }
}

template <class T>
void Genes<T>::Set(Functions &setup) throw(ExceptionMessenger) {
  /////////////////////////////
  /*
  TRandom3 rand;
  rand.SetSeed(time(NULL));
  for (Int_t i = 0; i < (setup.fNParam); ++i) {
    Double_t gene =
  rand.Uniform(setup.fInterval[i].first,setup.fInterval[i].second);
    fGenes.push_back(gene);
  }
  if (!setup)
    throw ExceptionMessenger("Do something with setup function!");
  */
  // Lets imagine that we have only one limit #0 for all parameters
  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  fGenes.resize(setup.fNParam);

  try {
    fGenes.resize(setup.fNParam);
  } catch (const std::length_error &le) {
    std::cerr << "Length error: " << le.what() << '\n';
  }

  std::uniform_real_distribution<T> dist(setup.fInterval[0].first,
                                         setup.fInterval[0].second);
  auto gen = std::bind(dist, mersenne_engine);
  std::generate(std::begin(fGenes), std::end(fGenes), gen);
  for (auto i : fGenes) {
    std::cout << "| " << i << " = element of gene |";
  }
}

template <class T> void Genes<T>::SetConstrain(Int_t i, T value) {
  fConstraines.emplace(fConstraines.begin() + i, value);
}

template <class T>
void Genes<T>::Evaluate(Functions &setup,
                        Genes<T> &ind) throw(ExceptionMessenger) {
  std::cout << "Again debug from Genes<T>::Evaluate():\n";
  printGenes(ind);
  (setup.evfunc)(ind);
  if (setup.fNCons) {
    ConstViol = 0;
  }
  fEvaluated = true;
}

template <class T> void Genes<T>::Clear(Option_t * /*option*/) {
  TObject::Clear();
  fGenes.clear();
  fRank = 0;
  fDominationCounter = 0.;
  fEvaluated = 0;
  fFitness.clear();
  fCrowdingDistance = 0;
  fDominated.clear();
  ConstViol = 0.;
}

template <class T>
T Genes<T>::CheckDominance(Functions *setup,
                           const Genes<T> *ind2) throw(ExceptionMessenger) {
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
    for (Int_t i = 0; i < setup->fNObjectives; ++i) {
      if (setup->fNObjectives > 1) {
        if (GetFitness(i) < ind2->GetFitness(i)) {
          fFlag1 = 1;
        } else if (GetFitness(i) > ind2->GetFitness(i)) {
          fFlag2 = 1;
        }
      } else {
        if (GetFitness(i) < ind2->GetFitness(i) &&
            fabs(GetFitness(i) - ind2->GetFitness(i)) > setup->fEpsilonC) {
          fFlag1 = 1;
        } else if (GetFitness(i) > ind2->GetFitness(i) &&
                   fabs(GetFitness(i) - ind2->GetFitness(i)) >
                       setup->fEpsilonC) {
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
template <class T> Int_t Genes<T>::Mutate() {
  TRandom rand;
  Double_t fRrnd, fDelta1, fDelta2, fMutPow, fDelta, fValue;
  Double_t y, yl, yu, xy;
  Int_t fNMut = 0;
  for (Int_t j = 0; j < setup->fNParam; ++j) {
    if (rand.Rndm() <= setup->fPMut) {
      y = fGenes[j];
      yl = setup->GetIntervalLimit(j).first;
      yu = setup->GetIntervalLimit(j).second;
      fDelta1 = (y - yl) / (yu - yl);
      fDelta2 = (yu - y) / (yu - yl);
      fRrnd = rand.Rndm();
      fMutPow = 1.0 / (setup->fEtaMut + 1.0);
      if (fRrnd <= 0.5) {
        xy = 1.0 - fDelta1;
        fValue = 2.0 * fRrnd +
                 (1.0 - 2.0 * fRrnd) * (pow(xy, (setup->fEtaMut + 1.0)));
        fDelta = pow(fValue, fMutPow) - 1.0;
      } else {
        xy = 1.0 - fDelta2;
        fValue = 2.0 * (1.0 - fRrnd) +
                 2.0 * (fRrnd - 0.5) * (pow(xy, (setup->fEtaMut + 1.0)));
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

template <class T>
void Genes<T>::WriteGenesTree(Genes<T> &ind, Population<T> &pop,
                              const char *file) {
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

template <class T>
void Genes<T>::UpdateGenesTree(Genes<T> &ind1, Genes<T> &ind2,
                               Population<T> &pop, const char *file) {
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
  f->Close();
}

template <class T>
void Genes<T>::ReadGenesTree(Genes<T> &ind, Population<T> &pop,
                             const char *file) {

  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("gvga");
  tree->SetBranchAddress("Population", &pop);
  Int_t entries = (Int_t)(tree->GetEntries());
}

// Ugly instantiating
template class Genes<Double_t>;
