#include <vector>
#include <ostream>
#include <string>
#include <algorithm>
#include <utility>

#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"

#include "Genes.h"
#include "Functions.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"

ClassImp(Population)
    /**
     * @brief Struct that is it crowding over objectives of
     variables though boolean comparing operator "<" on population for
     two individuals with different indexes and objectives/variables with index
     m
     */
    struct Comparing {
  Comparing(Population &p, Int_t indx) : pop(p), m(indx){};
  Population &pop;
  Int_t m;
  Bool_t operator()(Int_t i, Int_t j) {
    return pop.fCrowdingObj
               ? pop.GetGenes(i).GetFitness(i) < pop.GetGenes(j).GetFitness(j)
               : pop.GetGenes(i)[i] < pop.GetGenes(j)[j];
  };
};

void Population::Build() {
  TRandom3 rand;
  rand.SetSeed(time(NULL));
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it) {
    for (Int_t i = 1; i <= Functions::Instance()->GetNParam(); ++i) {
      for (const std::pair<Double_t, Double_t> &limit :
           (Functions::Instance()->GetInterval())) {
        Double_t value = rand.Uniform(limit.first, limit.second);
        Int_t position = std::distance(GetIndividuals().begin(), it);
        SetGenes(position, Genes::SetGene(i, value));
      }
    }
  }
  WritePopulationTree(*this, "NSGA.root");
}

void Population::CrowdingDistanceAll() {
  for (Int_t i = 0; (ULong_t)i < fFront.size(); ++i)
    CrowdingDistanceFront(i);
}

void Population::CrowdingDistanceFront(Int_t i) {
  Genes &F = fFront[i]; // Genes
  if (F.size() == 0)
    return;
  const Int_t l = F.size();
  for (Int_t i = 0; i < l; ++i)
    GetGenes(F[i]).SetCrowdingDistance(0);
  const Int_t limit =
      fCrowdingObj ? GetNObjectives() : Functions::Instance()->GetNParam();
  for (Int_t m = 0; m < limit; ++m) {
    std::sort(F.begin(), F.end(), Comparing(*this, m));
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

void Population::FastNonDominantSorting() {
  fFront.resize(1);
  fFront[0].clear();
  for (int i = 0; i < (ULong_t)fPopulation.size(); ++i) {
    std::vector<Double_t> fDom;
    Int_t fDomCount = 0;
    Genes &p = fPopulation[i];
    for (Int_t j = 0; (ULong_t)j < fPopulation.size(); ++j) {
      Genes &q = fPopulation[j];
      Int_t compare = p.Genes::CheckDominance(&q);
      if (compare == 1) {
        fDom.push_back(j);
      } else if (compare == -1) {
        fDomCount += 1;
      }
    }
    p.SetDominatedCounter(fDomCount);
    p.GetDominated().clear();
    p.GetDominated() = fDom;
    if (p.GetDominatedCounter() == 0) {
      p.SetRank(1);
      fFront[0].push_back(i);
    }
  }
  std::sort(fFront[0].begin(), fFront[0].end());
  int fi = 1;
  while (fFront[fi - 1].size() > 0) {
    Genes &fronti = fFront[fi - 1];
    std::vector<Double_t> Q;
    for (Int_t i = 0; (ULong_t)i < fronti.size(); ++i) {
      Genes &p = fPopulation[fronti[i]];
      for (Int_t j = 0; (ULong_t)j < p.GetDominated().size(); ++j) {
        Genes &q = fPopulation[p.GetDominated(j)];
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

void Population::Merge(const Population &population1,
                       const Population &population2) {
  std::copy(population1.fPopulation.begin(), population1.fPopulation.end(),
            fPopulation.begin());
  std::copy(population2.fPopulation.begin(), population2.fPopulation.end(),
            fPopulation.begin() + population1.GetPopulationSize());
}

void Population::Clear(Option_t * /*option*/) { ; } // Clear function

void Population::WritePopulationTree(Population &pop, const char *file) {
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

void Population::UpdatePopulationTree(Population &pop, const char *file) {
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

void Population::ReadPopulationTree(Population &pop, const char *file) {

  TFile *f = TFile::Open(file, "RECREATE");
  TTree *tree = (TTree *)f->Get("gvga");
  tree->SetBranchAddress("Population", &pop);
  Int_t entries = (Int_t)(tree->GetEntries());
}

// Testing work of stupid tree (could be deleted after)
Int_t Population::PrintTree(const char *file, const char *name) {
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

void Population::Evaluate() {
  for (auto it = GetIndividuals().begin(); it != GetIndividuals().end(); ++it)
  {
    //Functions::Instance()->SetFunctionGenes(Functions::Instance(), it); 
  }
}
