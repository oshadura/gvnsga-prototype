#include <vector>
#include <algorithm>
#include <cassert>

#include "TRandom.h"
#include "TRandom3.h"
//
#include "Functions.h"
#include "TGenes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"

AlgorithmNSGA *AlgorithmNSGA::fgNSGA2 = 0;

struct Sort {
  Population<Double_t> &pop;
  Sort(Population<Double_t> &population) : pop(population){};
  bool operator()(Int_t i, Int_t j) {
    const Genes<Double_t> &ind1 = pop.GetGenes(i);
    const Genes<Double_t> &ind2 = pop.GetGenes(j);
    if (ind1.GetRank() < ind2.GetRank())
      return true;
    else if (ind1.GetRank() == ind2.GetRank() &&
             ind1.GetCrowdingDistance() > ind2.GetCrowdingDistance())
      return true;
    return false;
  };
};

AlgorithmNSGA::AlgorithmNSGA()
    : function(0),popfunction(0), fPCross(0), fEtaCross(0), fPMut(0), fEtaMut(0), fNCross(0), fNMut(0),
      fNGen(0), fParentPop(), fChildPop(), fMixedPop() {
  fgNSGA2 = this;
}

AlgorithmNSGA::~AlgorithmNSGA() {
  if (fChildPop)
  {
    delete fChildPop;
    fChildPop = 0;
  }
  if (fParentPop)
  {
    delete fParentPop;
    fParentPop = 0;
  }
  if (fMixedPop){
    delete fMixedPop;
    fMixedPop = 0;
  }
}

AlgorithmNSGA *AlgorithmNSGA::Instance() {
  if (!fgNSGA2)
    AlgorithmNSGA *fgNSGA2 = new AlgorithmNSGA();
  return fgNSGA2;
}

void AlgorithmNSGA::Initialize() {
  std::cout << "Let's check NSGA2 configuration:"<< std::endl;


  Population<Double_t> *fChildPop = new Population<Double_t>();
  Population<Double_t> *fParentPop = new Population<Double_t>();
  Population<Double_t> *fMixedPop = new Population<Double_t>();
  //Missing check of input variables and creation of population with them 

  if(popfunction)
    fParentPop->SetPopFunction(popfunction);
    fChildPop->SetPopFunction(popfunction);
    fMixedPop->SetPopFunction(popfunction);

  fParentPop->Build();
  fParentPop->Evaluate();
  fParentPop->FastNonDominantSorting();
  fParentPop->CrowdingDistanceAll();
}

void AlgorithmNSGA::Selection(Population<Double_t> &oldpop,
                              Population<Double_t> &newpop) {
  static TRandom3 rand;
  const Int_t N = oldpop.GetPopulationSize();
  std::vector<Int_t> a1(N), a2(N);
  for (Int_t i = 0; i < N; ++i) {
    a1[i] = a2[i] = i;
  }
  for (Int_t i = 0; i < N; ++i) {
    std::swap(a1[rand.Uniform(i, N - 1)], a1[i]);
    std::swap(a2[rand.Uniform(i, N - 1)], a2[i]);
  }
  for (Int_t i = 0; i < N; i += 4) {
    Genes<Double_t> &p11 =
        Tournament(oldpop.GetGenes(a1[i]), oldpop.GetGenes(a1[i + 1]));
    Genes<Double_t> &p12 =
        Tournament(oldpop.GetGenes(a1[i + 2]), oldpop.GetGenes(a1[i + 3]));
    Crossover(p11, p12, newpop.GetGenes(i), newpop.GetGenes(i + 1));

    Genes<Double_t> &p21 =
        Tournament(oldpop.GetGenes(a2[i]), oldpop.GetGenes(a2[i + 1]));
    Genes<Double_t> &p22 =
        Tournament(oldpop.GetGenes(a2[i + 2]), oldpop.GetGenes(a2[i + 3]));
    Crossover(p21, p22, newpop.GetGenes(i + 2), newpop.GetGenes(i + 3));
  }
}

Genes<Double_t> &AlgorithmNSGA::Tournament(Genes<Double_t> &ind1,
                                           Genes<Double_t> &ind2) const {
  static TRandom rnd;
  Int_t fFlag = ind1.CheckDominance(&ind2);
  if (fFlag == 1) // Yes
    return ind1;
  else if (fFlag == -1) // Opposite
    return ind2;
  else if (ind1.GetCrowdingDistance() > ind2.GetCrowdingDistance())
    return ind1;
  else if (ind2.GetCrowdingDistance() > ind1.GetCrowdingDistance())
    return ind2;
  else if (rnd.Rndm() <= 0.5)
    return ind1;
  else
    return ind2;
}

void AlgorithmNSGA::Crossover(const Genes<Double_t> &parent1,
                              const Genes<Double_t> &parent2,
                              Genes<Double_t> &child1,
                              Genes<Double_t> &child2) {
  static TRandom rnd;
  Double_t y1, y2, yl, yu;
  Double_t c1, c2;
  Double_t alpha, beta, betaq;
  Int_t r = rnd.Rndm();
  if (r <= GetPCross()) {
    fNCross++;
    for (Int_t i = 0; i < Functions::Instance()->GetNParam(); i++) {
      if (fabs(parent1[i] - parent2[i]) > EPS) {
        if (parent1[i] < parent2[i]) {
          y1 = parent1[i];
          y2 = parent2[i];
        } else {
          y1 = parent2[i];
          y2 = parent1[i];
        }
        yl = Functions::Instance()->GetIntervalLimit(i).first;
        yu = Functions::Instance()->GetIntervalLimit(i).second;
        Int_t r = rnd.Rndm();
        beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
        alpha = 2.0 - pow(beta, -(GetEtaCross() + 1.0));
        if (r <= (1.0 / alpha)) { // This is a contracting crossover
          betaq = pow((r * alpha), (1.0 / (GetEtaCross() + 1.0)));
        } else { // This is an expanding crossover
          betaq = pow((1.0 / (2.0 - r * alpha)), (1.0 / (GetEtaCross() + 1.0)));
        }
        c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
        beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
        alpha = 2.0 - pow(beta, -(GetEtaCross() + 1.0));
        if (r <= (1.0 / alpha)) { // This is a contracting crossover
          betaq = pow((r * alpha), (1.0 / (GetEtaCross() + 1.0)));
        } else { // This is an expanding crossover
          betaq = pow((1.0 / (2.0 - r * alpha)), (1.0 / (GetEtaCross() + 1.0)));
        }
        c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
        c1 = fmin(fmax(c1, yl), yu);
        c2 = fmin(fmax(c2, yl), yu);
        if (rnd.Uniform() <= 0.5) {
          child1.SetGene(i, c2);
          child2.SetGene(i, c1);
        } else {
          child1.SetGene(i, c1);
          child2.SetGene(i, c2);
        }
      } else {
        child1.SetGene(i, parent1[i]);
        child2.SetGene(i, parent2[i]);
      }
    }
  } else {
    for (Int_t i = 0; i < Functions::Instance()->GetNParam(); i++) {
      child1.SetGene(i, parent1[i]);
      child2.SetGene(i, parent2[i]);
    }
  }
}

void AlgorithmNSGA::NextStep() {
  std::cout << "New generetion #" << (fParentPop->GetGenNumber()) + 1 << std::endl;
  Selection(*fParentPop, *fChildPop);
  fNMut = fChildPop->Mutate(); //not a std::pair (?)
  fChildPop->SetGenNumber((fChildPop->GetGenNumber()) + 1);
  fChildPop->Evaluate();
  // fNMut += fNMut;
  fMixedPop->Merge(*fParentPop, *fChildPop);
  fMixedPop->SetGenNumber((fMixedPop->GetGenNumber()) + 1);
  fMixedPop->FastNonDominantSorting();
  fParentPop->Clear();
  // until |Pt+1| + |Fi| <= N, until parent population is filled
  Int_t i = 0;
  while (fParentPop->GetPopulationSize() + (fMixedPop->GetFront(i)).size() <
         fMixedPop->GetPopulationSetupSize()) {
    Genes<Double_t> Fi = fMixedPop->GetFront(i);
    fMixedPop->CrowdingDistanceFront(i);            // calculate crowding in Fi
    for (Int_t j = 0; (Double_t)j < Fi.size(); ++j) // Pt+1 = Pt+1 U Fi
      fParentPop->GetIndividuals().push_back(fMixedPop->GetFront(j));
    i += 1;
  }
  fMixedPop->CrowdingDistanceFront(i); // calculate crowding in F
  std::sort(fMixedPop->GetFront(i).begin(), fMixedPop->GetFront(i).end(),
            Sort(*fMixedPop));
  const int extra =
      fParentPop->GetPopulationSetupSize() - fParentPop->GetPopulationSize();
  for (int j = 0; j < extra; ++j) // Pt+1 = Pt+1 U Fi[1:N-|Pt+1|]
    fParentPop->GetIndividuals().push_back(fMixedPop->GetFront(j));
  fParentPop->SetGenNumber((fParentPop->GetGenNumber()) + 1);
}

void AlgorithmNSGA::Evolution() {
  while (fNGen <= (Instance()->GetGenTotalNumber())) {
    NextStep();
  } // Check through population object
}

ClassImp(AlgorithmNSGA)
