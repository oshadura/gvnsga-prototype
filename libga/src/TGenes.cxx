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

templateClassImp(Genes)

    template <class T>
    Genes<T>::Genes() throw()
    : TObject(), fFitness(0), fRank(0), fDominationCounter(0),
      fCrowdingDistance(0), fEvaluated(false), fDominated(), ConstViol(0),
      fGenes(), fAllev(0), fBuffev(0), fThread(0), fPriority(0), fSteps(0),
      fVector(0), fTime(0), fMemory(0), setup(0), fConstraines() {}

template <class T>
Genes<T>::Genes(const Functions &config) throw(ExceptionMessenger)
    : TObject(), fFitness(), fRank(0), fDominationCounter(0),
      fCrowdingDistance(0), fEvaluated(false), fDominated(), ConstViol(0),
      fGenes(), fAllev(0), fBuffev(0), fThread(0), fPriority(0), fSteps(0),
      fVector(0), fTime(0), fMemory(0), setup(&config), fConstraines() {

  fGenes.resize(setup->fNParam, 0);
  fFitness.resize(setup->fNObjectives, 0);
  fConstraines.resize(setup->fNCons, 0);
}

template <class T> Genes<T> &Genes<T>::operator=(const Genes<T> &gen) {
  if (this != &gen) {
    fFitness = gen.fFitness;
    fRank = gen.fRank;
    fDominationCounter = gen.fDominationCounter;
    fCrowdingDistance = gen.fCrowdingDistance;
    fEvaluated = gen.fCrowdingDistance;
    fDominated = gen.fDominated;
    ConstViol = gen.ConstViol;
    fGenes = gen.fGenes;
    fAllev = gen.fAllev;
    fBuffev = gen.fBuffev;
    fThread = gen.fThread;
    fPriority = gen.fPriority;
    fSteps = gen.fSteps;
    fVector = gen.fVector;
    fTime = gen.fTime;
    fMemory = gen.fMemory;
    setup = gen.setup;
    fConstraines = gen.fConstraines;
  }
  return *this;
}

template <class T>
void Genes<T>::Set(Functions &setup, Genes<T> &ind) throw(ExceptionMessenger) {
  // Lets imagine that we have only one limit #0 for all parameters
  ind = Genes<T>(setup);
  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  ind.resize(setup.fNParam);
  std::uniform_real_distribution<T> dist(setup.fInterval[0].first,
                                         setup.fInterval[0].second);
  auto gen = std::bind(dist, mersenne_engine);
  std::generate(std::begin(ind), std::end(ind), gen);

  for (auto i : ind) {
    std::cout << "| " << i << " = element of gene |";
  }
}

#ifdef ENABLE_GEANTV
template <class T>
void Genes<T>::SetGeantV(Functions &setup,
                         Genes<T> &ind) throw(ExceptionMessenger) {
  // FIX GENERATORS
  // 1. Consider value that ([0] - 1) should be always smaller [0]
  // 2. Consider that [5] and [6] are generated in diferent way (new generator)
  ind = Genes<T>(setup);
  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  ind.resize(setup.fNParam);
  for (Int_t i = 0; i <= (setup.fNParam); ++i) {
    std::uniform_real_distribution<T> dist(setup.fInterval[i].first,
                                           setup.fInterval[i].second);
    auto gen = std::bind(dist, mersenne_engine);
    std::generate_n(std::begin(ind) + i, 1, gen);
  }

  for (auto i : ind) {
    std::cout << "| " << i << " = element of gene |";
  }
  /*
  TRandom rand;
  rand.SetSeed(5000);
  ind = Genes<T>(setup);
  ind.resize(setup.fNParam);
  rand.SetSeed(time(NULL));
  for (Int_t i = 0; i < (setup.fNParam); ++i) {
     fGenes[i] = rand.Uniform(setup .fInterval[i].first,
                             setup.fInterval[i].second);
     //fGenes.emplace(fGenes.begin() + i, value);
  }

  for (auto i : ind) {
    std::cout << "| " << i << " = element of gene |";
  }
  */
}
#endif

template <class T> void Genes<T>::SetConstrain(Int_t i, T value) {
  fConstraines.emplace(fConstraines.begin() + i, value);
}

//#ifdef ENABLE_GEANTV
/*
template <class T>
void Genes<T>::Evaluate(GeantPropagator *prop, Functions &setup,
                        Genes<T> &ind) throw(ExceptionMessenger) {
  std::cout << "-==============================================-" << std::endl;
  std::cout << "Again debug from Genes<T>::Evaluate():\n" << std::endl;
  printGenes(ind);
  (setup.evfunc)(prop, ind);
  if (setup.fNCons) {
    ConstViol = 0;
  }
  fEvaluated = true;
}
*/
// else
template <class T>
void Genes<T>::Evaluate(Functions &setup,
                        Genes<T> &ind) throw(ExceptionMessenger) {
  //std::cout << "-==============================================-" << std::endl;
  //std::cout << "Again debug from Genes<T>::Evaluate():\n" << std::endl;
  //printGenes(ind);
  (setup.evfunc)(ind);
  if (setup.fNCons) {
    ConstViol = 0;
  }
  fEvaluated = true;
}
//#endif

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
    Int_t fFlag1 = 0;
    Int_t fFlag2 = 0;
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
template <class T> Int_t Genes<T>::Mutate(const Functions *setup) {
  TRandom rand;
  Double_t fRnd, fDelta1, fDelta2, fMutPow, fDelta, fValue;
  Double_t y, LimitDown, LimitUp, xy;
  Int_t fNMut = 0;
  for (Int_t j = 0; j < setup->fNParam; ++j) {
    if (rand.Rndm() <= setup->fPMut) {
      y = fGenes[j];
      LimitDown = setup->GetIntervalLimit(j).first;
      LimitUp = setup->GetIntervalLimit(j).second;
      fDelta1 = (y - LimitDown) / (LimitUp - LimitDown);
      fDelta2 = (LimitUp - y) / (LimitUp - LimitDown);
      fRnd = rand.Rndm();
      fMutPow = 1.0 / (setup->fEtaMut + 1.0);
      if (fRnd <= 0.5) {
        xy = 1.0 - fDelta1;
        fValue =
            2.0 * fRnd + (1.0 - 2.0 * fRnd) * (pow(xy, (setup->fEtaMut + 1.0)));
        fDelta = pow(fValue, fMutPow) - 1.0;
      } else {
        xy = 1.0 - fDelta2;
        fValue = 2.0 * (1.0 - fRnd) +
                 2.0 * (fRnd - 0.5) * (pow(xy, (setup->fEtaMut + 1.0)));
        fDelta = 1.0 - (pow(fValue, fMutPow));
      }
      y = y + fDelta * (LimitUp - LimitDown);
      if (y < LimitDown)
        y = LimitDown;
      if (y > LimitUp)
        y = LimitUp;
      // std::cout << "Print new crossover part = " << y << std::endl;
      SetGene(j, y);
      fNMut += 1;
    }
  }
  return fNMut;
}

template <class T>
void Genes<T>::WriteGenesTree(Genes<T> &ind, Population<T> &pop,
                              const char *file) {
  R__LOAD_LIBRARY(libGa);
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
  R__LOAD_LIBRARY(libGa);
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
  R__LOAD_LIBRARY(libGa);
  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("gvga");
  tree->SetBranchAddress("Population", &pop);
  Int_t entries = (Int_t)(tree->GetEntries());
}

template class Genes<Double_t>;
