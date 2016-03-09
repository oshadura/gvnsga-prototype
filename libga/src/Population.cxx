#include <vector>
#include <ostream>
#include <string>
#include <utility>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <algorithm> // std::copy
#include <iostream>
#include <list>
#include <sstream>
#include <fstream>

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
    : fCrowdingObj(true), fPopFunction(NULL), setupPop(), fHisto(0) {
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

#ifdef ENABLE_GEANTV
template <class T> void Population<T>::Build() throw(ExceptionMessenger) {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    it->Genes<T>::SetGeantV(setupPop, *it);
    // fPopulation.emplace_back(&(*it).GetfGenes());
    std::cout << " Creating new individual.." << std::endl;
  }
 // WritePopulationTree(*this, "NSGA.root");
}
#else
template <class T> void Population<T>::Build() throw(ExceptionMessenger) {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    it->Genes<T>::Set(setupPop, *it);
    // fPopulation.emplace_back(&(*it).GetfGenes());
    std::cout << "Creating new individual.." << std::endl;
  }
 // WritePopulationTree(*this, "NSGA.root");
}
#endif

template <class T> void Population<T>::CrowdingDistanceAll() {
  // Doing for each front
  for (Int_t i = 0; i < fFront.size(); ++i)
    CrowdingDistanceFront(i);
}

template <class T> void Population<T>::CrowdingDistanceFront(Int_t front) {
  std::cout << "Crowding distance for Front<" << front << ">.." << std::endl;
  std::vector<Int_t> &F = fFront[front];
  /*
  for(auto &i : F){
    std::cout << i << std::endl;
  }
  */
  if (F.size() == 0)
    std::cout << BLUE << "FRONT size = 0" << RESET << std::endl;
  return;
  // Setting crowding distance as ZERO to all individuals
  for (Int_t i = 0; i < F.size(); ++i)
    GetGenes(F[i]).SetCrowdingDistance(0);
  // Pushing to go only for through fitness functions
  Int_t limit = fCrowdingObj ? setupPop.fNObjectives : setupPop.fNParam;
  // For each fitness function sort individuals according each front
  // Int_t limit = setupPop.fNObjectives;
  std::cout << "Limit for calculating CrowDist (setupPop.fNObjectives : "
               "setupPop.fNParam) = "
            << limit << std::endl;
  for (Int_t m = 0; m < limit; ++m) {
    std::sort(F.begin(), F.end(), Comparing<T>(*this, m));
    // Setting INF values for boundaries values (min, max)
    GetGenes(F[0]).SetCrowdingDistance(INF);
    if (F.size() > 1)
      GetGenes(F[F.size() - 1]).SetCrowdingDistance(INF);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Min in our front #0 = " << GetGenes(F[0]).GetGene(0)
              << std::endl;
    std::cout << "Max in our front #" << (F.size() - 1) << " = "
              << GetGenes(F[F.size() - 1]).GetGene(0) << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
    // Setting values of crowding distance for other elements
    for (Int_t i = 1; i < F.size() - 1; ++i) {
      if (GetGenes(F[i]).GetCrowdingDistance() != INF) {
        if (IsCrowdingObj() &&
            GetGenes(F[F.size() - 1]).GetFitness(m) !=
                GetGenes(F[0]).GetFitness(m)) {
          std::cout << "Crowding distance is defined over objectives..."
                    << std::endl;
          Double_t dist = (GetGenes(F[i + 1]).GetFitness(m) -
                           GetGenes(F[i - 1]).GetFitness(m)) /
                          (GetGenes(F[F.size() - 1]).GetFitness(m) -
                           GetGenes(F[0]).GetFitness(m));
          GetGenes(F[i]).SetCrowdingDistance(dist);
        } else if (!IsCrowdingObj() &&
                   GetGenes(F[F.size() - 1]).GetGene(m) !=
                       GetGenes(F[0]).GetGene(m)) {
          std::cout << "Crowding distance is defined over parameters.."
                    << std::endl;
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
    auto position = std::distance(fPopulation.begin(), it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout gene # " << position + 1
              << " after crowding distance steps:" << std::endl;
    Genes<T>::printGenes(*it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
  }
}

template <class T> void Population<T>::FastNonDominantSorting() {
  std::cout << "-==============================================-\n"
            << "Non-dominant sorting started - F[0] is filled.." << std::endl;
  fFront.resize(1);
  fFront[0].clear();
#pragma omp parallel for
  // Dominance checking for each Gene
  for (int i = 0; i < fPopulation.size(); ++i) {
    // std::cout << "\nSize of population = "<< fPopulation.size() << std::endl;
    // vector of individuals that dominated by p {selected individual}
    std::vector<Int_t> fDom;
    // number of individuals that dominate p {selected individual}
    Int_t fDomCount = 0;
    // Selected individual p (passing through loop)
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
      // std::cout << "fDomCount for " << i << " individual, counter value = "
      // << RED << fDomCount << RESET << std::endl;
      p.GetDominated().clear();
      p.SetDominated(fDom);
      /*
      for (auto i = fDom.begin(); i != fDom.end(); ++i){
        std::cout << CYAN << *i << ' ' << RESET;
      }
      std::cout << "\n" << std::endl;
      */
      if (p.GetDominatedCounter() == 0) {
        p.SetRank(1);
        fFront[0].push_back(i);
      }
    }
  }
  std::sort(fFront[0].begin(), fFront[0].end());
  int fi = 1;
  while (fFront[fi - 1].size() > 0) {
    std::cout << "Checking next fronts.." << std::endl;
    std::vector<Int_t> Q;
    std::vector<Int_t> &fronti = fFront[fi - 1];
    for (Int_t i = 0; i < fronti.size(); ++i) {
      Genes<T> &p = fPopulation[fronti[i]];
      for (Int_t j = 0; j < p.GetDominated().size(); ++j) {
        Genes<T> &q = fPopulation[p.GetDominated(j)];
        q.SetDominatedCounter((q.GetDominatedCounter() - 1));
        if (q.GetDominatedCounter() == 0) {
          q.SetRank(fi + 1);
          Q.push_back(p.GetDominated(j));
        }
      }
    }
    fi += 1;
    if (Q.size() > 0) {
      fFront.push_back(Q);
    } else {
      std::cout << "We have no more fronts to checkout.." << std::endl;
      return;
    }
  }
  std::cout << "-==============================================-" << std::endl;
  for (int i = 0; i < fFront.size(); ++i) {
    std::cout << RED << "fFront<" << i << "> = [";
    for (auto &it : fFront[i]) {
      std::cout << it << ' ';
    }
    std::cout << "]"
              << "\n" << RESET << std::endl;
  }
  std::cout << "-==============================================-" << std::endl;
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
  Int_t mutvalue;
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    const Functions *setupind = (*it).GetSetup();
    // std::cout << "-==============================================-" <<
    // std::endl;
    // std::cout << "So so, just to be sure -> number of objectives in "
    //             "Population::Mutation() "
    //          << (*it).GetSetup()->GetNObjectives() << std::endl;
    mutvalue += it->Genes<T>::Mutate(setupind);
    // this->printGenes(*it);
  }
  // std::cout << "Print number of mutations in Population: " << mutvalue <<
  // std::endl;
  return mutvalue;
}

template <class T>
void Population<T>::WritePopulationTree(Population &pop, const char *file) {
  fHisto = HistogramManager::Instance();
  //////////////////////////////////////
  if (!gSystem->AccessPathName(file, kFileExists)) {
    TFile *friendtree = new TFile("NSGApopulations.root", "RECREATE");
    TTree *treecopy = new TTree("GA", "Genetic Algorithm TTree");
    treecopy->Branch("Pop", "Pop", &pop);
    /*
    for (int i = 0; i < pop.GetPopulationSize(); ++i) {
      for (auto it = pop.GetGenes(i).begin(); it != pop.GetGenes(i).end();
           ++it) {
        treecopy->Branch("Genes", "Genes", &it);
      }
    }
    */
    ////////////////////////////////////
    TFile *f = TFile::Open("NSGA.root");
    if (f->IsZombie()) {
      std::cout << "Error opening file" << std::endl;
      exit(-1);
    }
    TTree *tree = (TTree *)f->Get("GA");
    tree->AddFriend("GA", "NSGApopulations.root");
    tree->Fill();
    ////////////////////////////////////
    // tree->Print();
    tree->Write();
    fHisto->HistoFill(pop, "NSGApopulations.root");
    f->cd();
    f->Close();
    friendtree->cd();
    friendtree->Close();
    gROOT->cd();
  } else {
    ///////////////////////////////////
    TFile *f = new TFile(file, "RECREATE");
    TTree *tree = new TTree("GA", "Genetic Algorithm TTree");
    tree->Branch("Pop", "Pop", &pop);
    /*
    for (int i = 0; i < pop.GetPopulationSize(); ++i) {
      for (auto it = pop.GetGenes(i).begin(); it != pop.GetGenes(i).end();
           ++it) {
        tree->Branch("Genes", "Genes", &it);
      }
    }
    */
    tree->Fill();
    tree->Write();
    fHisto->HistoFill(pop, const_cast<char *>(file));
    f->cd();
    f->Close();
    gROOT->cd();
    ///////////////////////////////////
  }
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
    printf("Cannot find tree  %s\n", tree);
    return -1;
  }
}

//#ifdef ENABLE_GEANTV
/*
template <class T> void Population<T>::Evaluate(GeantPropagator *prop) {
#ifdef ENABLE_OPENMP
  EvaluationOpenMP(prop);
#else
  Evaluation(prop);
#endif
}
*/
//#else
template <class T> void Population<T>::Evaluate() {
#ifdef ENABLE_OPENMP
  EvaluationOpenMP();
#else
  Evaluation();
#endif
}
//#endif

//#ifdef ENABLE_GEANTV
/*template <class T> void Population<T>::Evaluation(GeantPropagator *prop) {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    Genes<T>::Evaluate(prop, setupPop, *it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after sequence evaluation:" << std::endl;
    Genes<T>::printGenes(*it);
  }
}
#else
*/
template <class T> void Population<T>::Evaluation() {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    auto position = std::distance(fPopulation.begin(), it);
    Genes<T>::Evaluate(setupPop, *it);
    // std::cout << "-==============================================-"
    //          << std::endl;
    // std::cout << "Printout of gene " << position + 1 << " after sequence
    // evaluation:" << std::endl;
    // Genes<T>::printGenes(*it);
  }
}
//#endif

//#ifdef ENABLE_GEANTV
/*
template <class T> void Population<T>::EvaluationOpenMP(GeantPropagator *prop) {
#pragma omp parallel for
  for (int i = 0; i < GetPopulationSize(); ++i) {
    fPopulation[i].Evaluate(prop, setupPop, fPopulation[i]);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout after OPENMP evaluation:" << std::endl;
    Genes<T>::printGenes(fPopulation[i]);
  }
}
*/
//#else
template <class T> void Population<T>::EvaluationOpenMP() {
#pragma omp parallel for
  for (int i = 0; i < GetPopulationSize(); ++i) {
    fPopulation[i].Evaluate(setupPop, fPopulation[i]);
    // std::cout << "-==============================================-"
    //          << std::endl;
    // std::cout << "Printout after OPENMP evaluation:" << std::endl;
    // Genes<T>::printGenes(fPopulation[i]);
  }
}
//#endif

template <class T> void Population<T>::Print() {
  for (auto it = fPopulation.begin(); it != fPopulation.end(); ++it) {
    auto position = std::distance(fPopulation.begin(), it);
    std::cout << "-==============================================-"
              << std::endl;
    std::cout << "Printout of gene " << position + 1
              << " for population:" << std::endl;
    Genes<T>::printGenes(*it);
  }
}

template <class T>
void Population<T>::Store(const std::string &file, const Population<T> &pop) {
  std::ofstream *ofstream = new std::ofstream{file};
  // for (auto it = const_cast<Population<T> &>(pop).begin();it !=
  // const_cast<Population<T> &>(pop).end(); ++it) {
  for (int j = 0; j < pop.GetPopulationSetupSize(); ++j) {
    if (setupPop.fNObjectives > 0) {
      for (int i = 0; i < setupPop.fNObjectives; ++i) {
        const T &fitvalue =
            const_cast<Population<T> &>(pop).GetGenes(j).GetFitness(i);
        ofstream->write(reinterpret_cast<const char *>(&fitvalue),
                        sizeof(Double_t) * setupPop.fNObjectives);
      }
    }
    if (setupPop.fNCons > 0) {
      for (int i = 0; i < setupPop.fNCons; ++i) {
        const T &consvalue =
            const_cast<Population<T> &>(pop).GetGenes(j).GetConstrain(i);
        ofstream->write(reinterpret_cast<const char *>(&consvalue),
                        sizeof(Double_t) * setupPop.fNCons);
      }
    }
    if (setupPop.fNParam > 0) {
      for (int i = 0; i < setupPop.fNParam; ++i) {
        const T &value =
            const_cast<Population<T> &>(pop).GetGenes(j).GetGene(i);
        ofstream->write(reinterpret_cast<const char *>(&value),
                        sizeof(Double_t) * setupPop.fNParam);
      }
    }
    const T &constviol =
        const_cast<Population<T> &>(pop).GetGenes(j).GetConsViol();
    const T &cdist =
        const_cast<Population<T> &>(pop).GetGenes(j).GetCrowdingDistance();
    const Int_t &rank = const_cast<Population<T> &>(pop).GetGenes(j).GetRank();
    ofstream->write(reinterpret_cast<const char *>(&constviol),
                    sizeof(Double_t));
    ofstream->write(reinterpret_cast<const char *>(&rank), sizeof(Int_t));
    ofstream->write(reinterpret_cast<const char *>(&cdist), sizeof(Double_t));
  }
}

// Ugly instantiation
template class Population<Double_t>;
