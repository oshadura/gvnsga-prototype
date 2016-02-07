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
#include "TSystem.h"

#include "TGenes.h"
#include "Functions.h"
#include "Population.h"
#include "HistogramManager.h"

templateClassImp(Population)

    template <class T>
    Population<T>::Population(
        const Int_t fSizePop, const Int_t fNParam, const Int_t fNCons,
        const Int_t fNObjectives, const Double_t fEpsilonC,
        const Double_t fPMut, const Double_t fEtaMut,
        const std::vector<std::pair<Double_t, Double_t>> fInterval,
        const Functions::functype func) throw(ExceptionMessenger)
    : fCrowdingObj(true), fPopFunction(NULL), setupPop() {
  setupPop.fNParam = fNParam;
  setupPop.fInterval = fInterval;
  setupPop.fNCons = fNCons;
  setupPop.fNObjectives = fNObjectives;
  setupPop.fEpsilonC = fEpsilonC;
  setupPop.fPMut = fPMut;
  setupPop.fEtaMut = fEtaMut;
  setupPop.evfunc = func;
  for (int i = 0; i < fSizePop; ++i) {
    fPopulation.emplace_back(Genes<T>(setupPop));
  }
}

/**
 * @brief Struct that is it crowding over objectives of
 variables though boolean comparing operator "<" on population for
 two individuals with different indexes and objectives/variables with index
 m
 */
template <typename T> struct Comparing {
  Comparing(const Population<T> &p, Int_t index) : pop(p), m(index){};
  const Population<T> &pop;
  Int_t m;
  Bool_t operator()(Int_t i, Int_t j) {
    return pop.fCrowdingObj
               ? const_cast<Population<Double_t> &>(pop).GetGenes(i).GetFitness(
                     m) < const_cast<Population<Double_t> &>(pop)
                              .GetGenes(j)
                              .GetFitness(j)
               : const_cast<Population<Double_t> &>(pop).GetGenes(i).GetGene(
                     m) < const_cast<Population<Double_t> &>(pop)
                              .GetGenes(i)
                              .GetGene(j);
  };
};

template <class T> void Population<T>::Build() throw(ExceptionMessenger) {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    it->Genes<T>::SetGeantV(setupPop, *it);
    // fPopulation.emplace_back(&(*it).GetfGenes());
    std::cout << " Creating new individual.." << std::endl;
  }
  // WritePopulationTree(*this, "NSGA.root");
}

template <class T> void Population<T>::CrowdingDistanceAll() {
  for (Int_t i = 0; i < fFront.size(); ++i)
    CrowdingDistanceFront(i);
}

template <class T> void Population<T>::CrowdingDistanceFront(Int_t front) {
  std::vector<Int_t> &F = fFront[front];
  /*
  for(auto &i : F){
    std::cout << i << std::endl;
  }
  */
  if (F.size() == 0)
    return;
  for (Int_t i = 0; i < F.size(); ++i)
    GetGenes(F[i]).SetCrowdingDistance(0);
  Int_t limit = fCrowdingObj ? setupPop.fNObjectives : setupPop.fNParam;
  // std::cout << "Limit for calculating CrowDist (setupPop.fNObjectives :
  // setupPop.fNParam) = " << limit << std::endl;
  for (Int_t m = 0; m < limit; ++m) {
    std::sort(F.begin(), F.end(), Comparing<T>(*this, m));
    GetGenes(F[0]).SetCrowdingDistance(INF);
    if (F.size() > 1)
      GetGenes(F[F.size() - 1]).SetCrowdingDistance(INF);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Min in our front = " << GetGenes(F[0]).GetGene(0)
              << std::endl;
    std::cout << "Max in our front = " << GetGenes(F[F.size() - 1]).GetGene(0)
              << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
    for (Int_t i = 1; i < F.size() - 1; ++i) {
      if (GetGenes(F[i]).GetCrowdingDistance() != INF) {
        if (IsCrowdingObj() &&
            GetGenes(F[F.size() - 1]).GetFitness(m) !=
                GetGenes(F[0]).GetFitness(m)) {
          Double_t dist = (GetGenes(F[i + 1]).GetFitness(m) -
                           GetGenes(F[i - 1]).GetFitness(m)) /
                          (GetGenes(F[F.size() - 1]).GetFitness(m) -
                           GetGenes(F[0]).GetFitness(m));
          GetGenes(F[i]).SetCrowdingDistance(dist);
        } else if (!IsCrowdingObj() &&
                   GetGenes(F[F.size() - 1]).GetGene(m) !=
                       GetGenes(F[0]).GetGene(m)) {
          Double_t dist =
              (GetGenes(F[i + 1]).GetGene(m) - GetGenes(F[i - 1]).GetGene(m)) /
              (GetGenes(F[F.size() - 1]).GetGene(m) -
               GetGenes(F[0]).GetGene(m));
          GetGenes(F[i]).SetCrowdingDistance(dist);
        }
      }
    }
  }
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after crowding distance steps:" << std::endl;
    Genes<T>::printGenes(*it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
  }
}

template <class T> void Population<T>::FastNonDominantSorting() {
  // std::cout << fPopulation.size() <<std::endl;
  fFront.resize(1);
  fFront[0].clear();
#pragma omp parallel for
  // Dominance checking for each Gene
  for (int i = 0; i < fPopulation.size(); ++i) {
    std::vector<Int_t> fDom;
    Int_t fDomCount = 0;
    Genes<T> &p = fPopulation[i];
    for (Int_t j = 0; j < fPopulation.size(); ++j) {
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
      p.SetDominated(fDom);
      if (p.GetDominatedCounter() == 0) {
        p.SetRank(1);
        fFront[0].push_back(i);
      }
    }
  }
  std::sort(fFront[0].begin(), fFront[0].end());
  int fi = 1;
  while (fFront[fi - 1].size() > 0) {
    std::vector<Int_t> &fronti = fFront[fi - 1];
    std::vector<Int_t> Q;
    for (Int_t i = 0; i < fronti.size(); ++i) {
      Genes<T> &p = fPopulation[fronti[i]];
      for (Int_t j = 0; j < p.GetDominated().size(); ++j) {
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
  std::cout << "-==============================================-" << std::endl;
  std::cout << "Front for Population fFront<T> = [";
  for (int i = 0; i < fFront.size(); ++i) {
    for (auto &it : fFront[i]) {
      std::cout << it << ' ';
    }
  }
  std::cout << "]"
            << "\n" << std::endl;
}

template <class T>
void Population<T>::Merge(const Population<T> &population1,
                          const Population<T> &population2) {
  std::copy(population1.fPopulation.begin(), population1.fPopulation.end(),
            fPopulation.begin());
  std::copy(population2.fPopulation.begin(), population2.fPopulation.end(),
            fPopulation.begin() +
                const_cast<Population<T> &>(population1).GetPopulationSize());
}

template <class T> Int_t Population<T>::Mutate() {
  Int_t tmp;
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    const Functions *setupind = (*it).GetSetup();
    // std::cout << "-==============================================-"
    //          << std::endl;
    // std::cout << "So so, just to be sure -> number of objectives in "
    //             "Population::Mutation() "
    //          << (*it).GetSetup()->GetNObjectives() << std::endl;
    tmp += it->Genes<T>::Mutate(setupind);
    // this->printGenes(*it);
  }
  return tmp;
}

template <class T>
void Population<T>::WritePopulationTree(Population &pop, const char *file) {
  // if (!file) {
  gSystem->Load("libga/libGA.so");
  TFile *f = new TFile(file, "RECREATE");
  TTree *tree = new TTree("Population", "Genetic Algorithm TTree");
  tree->Branch("Population", "Population", &pop);
  for (int i = 0; i < pop.GetPopulationSize(); ++i) {
    for (auto it = pop.GetGenes(i).begin(); it != pop.GetGenes(i).end(); ++it) {
      tree->Branch("Genes", "Genes", &it);
    }
  }
  tree->Fill();
  tree->Print();
  // tree->Write();
  //}
  /*
  else {
    TFile *f = TFile::Open(file, "RECREATE");
    if (f->IsZombie()) {
    std::cout << "Error opening file" << std::endl;
    exit(-1);
    }
  */
  // Will be uncommented after installation of latest ROOT (> 15 September)
  /*
  else {
    TFile *ffriend = new TFile("NSGA-friend.root", "RECREATE");
    TTree *treecopy = new TTree("Population", "Genetic Algorithm TTree");
    treecopy->Branch("Population", "Population", &pop);
    TFile *f = new TFile("NSGA.root");
    TTree *tree = (TTree *)f->Get("tree");
    tree->AddFriend("treecopy", "NSGA-friend.root");
    tree->Fill();
    tree->Print();
  }
  */
}

template <class T>
void Population<T>::UpdatePopulationTree(Population<T> &pop, const char *file) {
  // Looks it is not possible update existing events, lets update the tree
  TFile *f = TFile::Open(file, "RECREATE");
  if (!f) {
    return;
  }
  TTree *tree = (TTree *)f->Get("Population");
  TTree *output =
      new TTree("Population", "Changing individual in a population");
  output->Branch("Population", &pop);
  tree->AddFriend("gvga_1", f);
  f->Write();
}

template <class T>
void Population<T>::ReadPopulationTree(Population<T> &pop, const char *file) {

  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("Population");
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

#ifdef ENABLE_GEANTV
template <class T> void Population<T>::Evaluate(GeantPropagator *prop) {
#ifdef ENABLE_OPENMP
  EvaluationOpenMP(prop);
#else
    Evaluation(prop);
#endif
}
#else
template <class T> void Population<T>::Evaluate() {
#ifdef ENABLE_OPENMP
  EvaluationOpenMP();
#else
    Evaluation();
#endif
}
#endif


#ifdef ENABLE_GEANTV
template <class T> void Population<T>::Evaluation(GeantPropagator* prop) {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    Genes<T>::Evaluate(prop, setupPop, *it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after sequence evaluation:" << std::endl;
    Genes<T>::printGenes(*it);
  }
}
#else
template <class T> void Population<T>::Evaluation() {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    Genes<T>::Evaluate(setupPop, *it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after sequence evaluation:" << std::endl;
    Genes<T>::printGenes(*it);
  }
}
#endif

#ifdef ENABLE_GEANTV
template <class T> void Population<T>::EvaluationOpenMP(GeantPropagator* prop) {
#pragma omp parallel for
  for (int i = 0; i < GetPopulationSize(); ++i) {
    fPopulation[i].Evaluate(prop, setupPop, fPopulation[i]);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after OPENMP evaluation:" << std::endl;
    Genes<T>::printGenes(fPopulation[i]);
  }
}
#else
template <class T> void Population<T>::EvaluationOpenMP() {
#pragma omp parallel for
  for (int i = 0; i < GetPopulationSize(); ++i) {
    fPopulation[i].Evaluate(setupPop, fPopulation[i]);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after OPENMP evaluation:" << std::endl;
    Genes<T>::printGenes(fPopulation[i]);
  }
}
#endif

// Ugly instantiation
template class Population<Double_t>;
