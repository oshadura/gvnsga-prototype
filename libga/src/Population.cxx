#include <vector>
#include <ostream>
#include <string>
#include <utility>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <algorithm> // std::copy

#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "TGenes.h"
#include "Functions.h"
#include "Population.h"
#include "HistogramManager.h"

// ClassImp(Population<T>)

  template <class T> Population<T>::Population(const Int_t fSizePop,
    const Int_t fNParam,
    const Int_t fNCons,
    const Int_t fNObjectives,
    const Double_t fEpsilonC,
    const Double_t fPMut,
    const Double_t fEtaMut,
    const std::vector<std::pair<Double_t, Double_t>> fInterval,
    const Functions::functype func) 
  throw (ExceptionMessenger) : fCrowdingObj(true), fPopFunction(NULL), setupPop(){

    setupPop.fNParam = fNParam;
    setupPop.fInterval = fInterval;
    setupPop.fNCons = fNCons;
    setupPop.fNObjectives = fNObjectives;
    setupPop.fEpsilonC = fEpsilonC;
    setupPop.fPMut = fPMut;
    setupPop.fEtaMut = fEtaMut;
    setupPop.evfunc = func;

    for (int i = 0; i < fSizePop; ++i)
    { 
      Genes<T> genes(&setupPop);
      //Genes<T>::printGenes(genes);
      fPopulation.push_back(genes);
    }

  }

/**
 * @brief Struct that is it crowding over objectives of
 variables though boolean comparing operator "<" on population for
 two individuals with different indexes and objectives/variables with index
 m
 */
template <typename T> struct Comparing {
  Comparing(Population<T> &p, Int_t indx) : pop(p), m(indx){};
  Population<T> &pop;
  Int_t m;
  Bool_t operator()(Int_t i, Int_t j) {
    return pop.fCrowdingObj
               ? pop.GetGenes(i).GetFitness(i) < pop.GetGenes(j).GetFitness(j)
               : pop.GetGenes(i)[i] < pop.GetGenes(j)[j];
  };
};

/*
template <class T> void Population<T>::Build() {
  TRandom3 rand;
  rand.SetSeed(time(NULL));
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it) {
    for (Int_t i = 1; i <= Functions::Instance()->GetNParam(); ++i) {
      for (const std::pair<Double_t, Double_t> &limit :
           (Functions::Instance()->GetInterval())) {
        Double_t value = rand.Uniform(limit.first, limit.second);
        Int_t position = std::distance(GetIndividuals().begin(), it);
        SetGenes(position, Genes<T>::SetGene(i, value));
      }
    }
  }
  WritePopulationTree(*this, "NSGA.root");
}
*/

template <class T> void Population<T>::Build() throw (ExceptionMessenger) {
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it) {
    //it->Genes<T>::Set();
    it->Genes<T>::Set(setupPop);
    std::cout << " Creating new individual.." << std::endl;
    }
  //WritePopulationTree(*this, "NSGA.root");
}

template <class T> void Population<T>::CrowdingDistanceAll() {
  for (Int_t i = 0; (ULong_t)i < fFront.size(); ++i)
    CrowdingDistanceFront(i);
}

template <class T> void Population<T>::CrowdingDistanceFront(Int_t i) {
  Genes<T> &F = fFront[i]; // Genes
  if (F.size() == 0)
    return;
  Int_t l = F.size();
  for (Int_t i = 0; i < l; ++i)
    GetGenes(F[i]).SetCrowdingDistance(0);
  Int_t limit = fCrowdingObj ? setupPop.fNObjectives
                             : setupPop.fNParam;
  for (Int_t m = 0; m < limit; ++m) {
    std::sort(F.begin(), F.end(), Comparing<T>(*this, m));
    GetGenes(F[0]).SetCrowdingDistance(INF);
    if (l > 1)
      GetGenes(F[l - i]).SetCrowdingDistance(INF);
    for (Int_t i = 1; i < l - 1; ++i) {
      if (GetGenes(F[i]).GetCrowdingDistance() != INF) {
        if (IsCrowdingObj() &&
            GetGenes(F[l - 1]).GetFitness(m) != GetGenes(F[0]).GetFitness(m)) {
          Double_t dist =
              (GetGenes(F[i + 1]).GetFitness(m) -
               GetGenes(F[i - 1]).GetFitness(m)) /
              (GetGenes(F[l - 1]).GetFitness(m) - GetGenes(F[0]).GetFitness(m));
          GetGenes(F[i]).SetCrowdingDistance(dist);
        } else if (!IsCrowdingObj() &&
                   GetGenes(F[l - 1])[m] != GetGenes(F[0])[m]) {
          Double_t dist = (GetGenes(F[i + 1])[m] - GetGenes(F[i - 1])[m]) /
                          (GetGenes(F[l - 1])[m] - GetGenes(F[0])[m]);
          GetGenes(F[i]).SetCrowdingDistance(dist);
        }
      }
    }
  }
}

template <class T> void Population<T>::FastNonDominantSorting() {
  fFront.resize(1);
  fFront[0].clear();
#pragma omp parallel for
  for (int i = 0; i < (ULong_t)fPopulation.size(); ++i) {
    std::vector<Double_t> fDom;
    Int_t fDomCount = 0;
    Genes<T> &p = fPopulation[i];
    for (Int_t j = 0; (ULong_t)j < fPopulation.size(); ++j) {
      Genes<T> &q = fPopulation[j];
      Int_t compare = p.Genes<T>::CheckDominance(&setupPop, &q);
      if (compare == 1) {
        fDom.push_back(j);
      } else if (compare == -1) {
        fDomCount += 1;
      }
    }
#pragma omp critical
    {
    p.SetDominatedCounter(fDomCount);
    p.GetDominated().clear();
    p.GetDominated() = fDom;
    if (p.GetDominatedCounter() == 0) {
      p.SetRank(1);
      fFront[0].push_back(i);
    }}
  }
  std::sort(fFront[0].begin(), fFront[0].end());
  int fi = 1;
  while (fFront[fi - 1].size() > 0) {
    Genes<T> &fronti = fFront[fi - 1];
    Genes<T> Q;
    for (Int_t i = 0; (ULong_t)i < fronti.size(); ++i) {
      Genes<T> &p = fPopulation[fronti[i]];
      for (Int_t j = 0; (ULong_t)j < p.GetDominated().size(); ++j) {
        Genes<T> &q = fPopulation[p.GetDominated(j)];
        q.SetDominatedCounter(-1); // -= 1;
        if (q.GetDominatedCounter() == 0) {
          q.SetRank(fi + 1);
          Q.push_back(p.GetDominated(j));
        }
      }
    }
    fi += 1;
    fFront.push_back(Q);
  }
}

template <class T>
void Population<T>::Merge(const Population<T> &population1,
                          const Population<T> &population2) {
  std::copy(population1.fPopulation.begin(), population1.fPopulation.end(),
            fPopulation.begin());
  std::copy(population2.fPopulation.begin(), population2.fPopulation.end(),
            fPopulation.begin() + population1.GetPopulationSize());
}

template <class T> void Population<T>::Clear(Option_t * /*option*/) {
  ;
} // Clear function

template <class T>
void Population<T>::WritePopulationTree(Population &pop, const char *file) {
  TStopwatch Timer;
  if (!file) {
    TFile *f = new TFile(file, "RECREATE");
    TTree *tree = new TTree("gvga", "Genetic Algorithm TTree");
    tree->Branch("Population", &pop);
    f->Write();
    tree->Fill();
    tree->Print();
  } else {
    TFile *f = TFile::Open(file, "RECREATE");
    TTree *tree = (TTree *)f->Get("gvga");
    TTree *tr_c = new TTree("gvga_1", "Genetic Algorithm friend Tree");
    tree->AddFriend("gvga_1", f);
    f->Write();
    tree->Fill();
    tree->Print();
  }
}

template <class T> Int_t Population<T>::Mutate(){
   Int_t tmp;
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it)
  {
    tmp += it-> Genes<T>::Mutate();
  }
  return tmp;
}

template <class T>
void Population<T>::UpdatePopulationTree(Population<T> &pop, const char *file) {
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

template <class T>
void Population<T>::ReadPopulationTree(Population<T> &pop, const char *file) {

  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("gvga");
  tree->SetBranchAddress("Population", &pop);
  Int_t entries = (Int_t)(tree->GetEntries());
}

// Testing work of stupid tree (could be deleted after)
template <class T>
Int_t Population<T>::PrintTree(const char *file, const char *name) {
  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = 0;
  f->GetObject(name, tree);
  if (tree) {
    tree->Print();
    return 0;
  } else {
    Error("PrintTree()", "Cannot find tree  %s", tree);
    return -1;
  }
}

template <class T> void Population<T>::Evaluate() {
#ifdef ENABLE_OPENMP
#pragma omp parrallel for
  for (int i = 0; i < GetPopulationSize(); ++i) {
    auto ind = GetGenes(i);
    //Genes<T>::Evaluate(ind);
    ind.Evaluate(&setupPop);
  }
#else
  /*
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it) {
    //Genes<T>::Evaluate(*it);
    it.Evaluate(&setupPop);

  }
  */
#endif
}

// Ugly instantiation
template class Population<Double_t>;
