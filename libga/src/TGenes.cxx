#include <vector>
#include <ostream>
#include <string>
#include <utility>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <ctime>
// Map for Genes[x]<->Limits[x] (?)
#include <map>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

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
///////////////////
#include "Process.h"


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
  // Creating it but empty, working after with push_back
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
  ind.empty();
  ind.reserve(setup.fNParam);
  ind = Genes<T>(setup);
  std::random_device rnd_device;
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 mersenne_engine(seed);
  //std::mt19937 mersenne_engine(rnd_device());
  for (Int_t i = 0; i < (setup.fNParam); ++i) {
    std::uniform_real_distribution<T> dist(setup.fInterval[i].first,
                                           setup.fInterval[i].second);
    auto gen = std::bind(dist, std::ref(mersenne_engine));
    // std::generate(std::begin(ind), std::end(ind), gen);
    ind.SetGene(i, gen());
  }
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
  // 2. Consider that [5] and [6] are generated in diferent way (new
  // generator)
  ind.empty();
  ind.reserve(setup.fNParam);
  ind = Genes<T>(setup);
  std::random_device rnd_device;
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  //std::mt19937 mersenne_engine(rnd_device());
  std::mt19937 mersenne_engine(seed);
  for (Int_t i = 0; i < (setup.fNParam); ++i) {
    std::uniform_real_distribution<T> dist(setup.fInterval[i].first,
                                           setup.fInterval[i].second);
    auto gen = std::bind(dist, std::ref(mersenne_engine));
    // std::generate_n(std::begin(ind) + i, 1, gen());
    ind.SetGene(i, gen());

  }
  for (auto i : ind) {
    std::cout << "| " << i << " = element of gene |";
  }
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
  std::cout << "-==============================================-" <<
std::endl;
  std::cout << "Again debug from Genes<T>::Evaluate():\n" << std::endl;
  printGenes(ind);
  (setup.evfunc)(prop, ind);
  if (setup.fNCons) {
    ConstViol = 0;
  }
  fEvaluated = true;
}
*/


  /**
   * @brief [brief description]
   * @details 
  Create pipe
  Fork
  In parent:
    Close read end of pipe
    Write data to be evaluated down write end of pipe
    Close write end of pipe
    Wait for child to die
  In child
    Close write end of pipe
    Duplicate read end of pipe to stdin
    Close read end of pipe
    Exec the evaluate program
    Exit with an error if the exec returns
   */
template <class T>
void Genes<T>::Evaluate(Functions &setup,
                        Genes<T> &ind) throw(ExceptionMessenger) {

  // std::cout << "-==============================================-" <<
  // std::endl;
  // std::cout << "Again debug from Genes<T>::Evaluate():\n" << std::endl;
  // printGenes(ind);
  // Lets consider that we have inly 2 pipes
  size_t sizeofFitness = sizeof(fFitness) + sizeof(T)* fFitness.capacity();
  const int fNumberChildren = 1;
  int pipeGA[fNumberChildren + 1];
  //Array of pids to be dead
  pid_t fArrayDead[fNumberChildren]; 
  pid_t cpid;
  pipe(pipeGA);
  cpid = fork();
  for (int i = 0; i < fNumberChildren; ++i){
    if (cpid == 0) {
    std::cout << "Starting child.."<< std::endl;
    close(pipeGA[1]); // close the write-end of the pipe 
    (setup.evfunc)(ind);
    while (read(pipeGA[0], &fFitness, sizeofFitness) > 0){
      write(pipeGA[1], &fFitness, 1);
    }
    close(pipeGA[0]); // close the read-end of the pipe
    exit(EXIT_SUCCESS);
  }
  else{
    fArrayDead[i] = cpid;
    close(pipeGA[0]);
    write(pipeGA[1], &fFitness, sizeofFitness);
    close(pipeGA[1]); // close the read-end of the pipe
    for (int i = 0; i < fNumberChildren; ++i){
      std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.." << std::endl;
      waitpid(fArrayDead[i], NULL, 0);
      std::cout << "PID: " << fArrayDead[i] << " has shut down.." << std::endl;
    }
    //wait(NULL);
    std::cout << "WE ARE BACK TO MASTER JOB::"<< std::endl;
  }
  }
  if (setup.fNCons) {
    ConstViol = 0;
  }
  fEvaluated = true;
  // Cleaning array of previos pids
  fArrayDead[fNumberChildren] = NULL;
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
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> rand(0, 1);
  Double_t fRnd, fDelta1, fDelta2, fMutPow, fDelta, fValue;
  Double_t y, LimitDown, LimitUp, xy;
  Int_t fNMut = 0;
  for (Int_t j = 0; j < setup->fNParam; ++j) {
    if (rand(gen) <= setup->fPMut) {
      y = fGenes[j];
      LimitDown = setup->GetIntervalLimit(j).first;
      LimitUp = setup->GetIntervalLimit(j).second;
      fDelta1 = (y - LimitDown) / (LimitUp - LimitDown);
      fDelta2 = (LimitUp - y) / (LimitUp - LimitDown);
      fRnd = rand(gen);
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
      // std::cout << "New part of Gene<T> for Mutation() " << y << std::endl;
      if (y < LimitDown)
        y = LimitDown;
      if (y > LimitUp)
        y = LimitUp;
      SetGene(j, y);
      fNMut = fNMut + 1;
    }
  }
  // std::cout << "Print number of mutations in Gene: " << fNMut << std::endl;
  return fNMut;
}

template <class T>
void Genes<T>::WriteGenesTree(Genes<T> &ind, Population<T> &pop,
                              const char *file) {
  //R__LOAD_LIBRARY(libGa);
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
  //R__LOAD_LIBRARY(libGa);
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
  //R__LOAD_LIBRARY(libGa);
  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("gvga");
  tree->SetBranchAddress("Population", &pop);
  Int_t entries = (Int_t)(tree->GetEntries());
}

template class Genes<Double_t>;
