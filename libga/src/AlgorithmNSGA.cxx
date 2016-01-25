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

// AlgorithmNSGA *AlgorithmNSGA::fgNSGA2 = 0;
ClassImp(AlgorithmNSGA)

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
    : function(0), popfunction(0), fPCross(0), fEtaCross(0), fNCross(0),
      fNMut(0), fNGen(0), fParentPop(0), fChildPop(0), fMixedPop(0), fSizePop(0),
      fNParam(0), fInterval(0), fNCons(0), fNObjectives(0), fPMut(0),
      fEtaMut(0), fEpsilonC(0), fCrowdingObj(true) {}

AlgorithmNSGA::~AlgorithmNSGA() {
  if (fChildPop) {
    delete fChildPop;
    fChildPop = 0;
  }
  if (fParentPop) {
    delete fParentPop;
    fParentPop = 0;
  }
  if (fMixedPop) {
    delete fMixedPop;
    fMixedPop = 0;
  }
}

void AlgorithmNSGA::Initialize() throw(ExceptionMessenger) {

  std::cout << "Let's check NSGA2 configuration..." << std::endl;

  if (fSizePop < 0)
    throw ExceptionMessenger("You didn't enter size of Population");
  if (fNParam < 0)
    throw ExceptionMessenger("And what is about number of parameters?");
  if (fNCons < 0)
    throw ExceptionMessenger(
        "Don't forget to put number constraints (at least 0)");
  if (fNObjectives <= 0)
    throw ExceptionMessenger(
        "What is about number of objectives? How you plan to measure things?");
  if (fEpsilonC <= 0)
    throw ExceptionMessenger("EpsilonC is not known for me!");
  if (fPMut <= 0)
    throw ExceptionMessenger(
        "Are you joking? Put probability of mutation as fast as possible..");
  if (fEtaMut < 0)
    throw ExceptionMessenger(
        "Are you joking? Put eta value of mutation as fast as possible.");
  // if (fInterval.size() !=  fNParam)
  //  throw ExceptionMessenger("Interval for generation individuals is not
  //  setuped at all!");
  if (function == 0)
    throw ExceptionMessenger(
        "Here I will not talk anymore with you: no function - no job");

  fChildPop =
      new Population<Double_t>(fSizePop, fNParam, fNCons, fNObjectives,
                               fEpsilonC, fPMut, fEtaMut, fInterval, function);
  fParentPop =
      new Population<Double_t>(fSizePop, fNParam, fNCons, fNObjectives,
                               fEpsilonC, fPMut, fEtaMut, fInterval, function);
  fMixedPop =
      new Population<Double_t>(fSizePop * 2, fNParam, fNCons, fNObjectives,
                               fEpsilonC, fPMut, fEtaMut, fInterval, function);
  // Missing check of input variables and creation of population with them
  // Report(configuration);
  std::cout<< "-==============================================-"<<std::endl;
  std::cout << "Population size = " << fSizePop
            << "\nNumber of generations = " << fNGen
            << "\nNumber of objective functions = " << fNObjectives
            << "\nNumber of constraints = " << fNCons
            << "\nNumber of variables = " << fNParam << std::endl;

  if (fNParam != 0) {
    for (int i = 0; i < fNParam; ++i) {
      std::cout << "\nLower limit of real variable " << (i + 1) << " = "
                << fInterval[i].first << std::endl;
      std::cout << "\nUpper limit of real variable " << (i + 1) << " = "
                << fInterval[i].second << std::endl;
    }
    std::cout << "\nProbability of crossover of real variable = " << fPCross
              << std::endl;
    std::cout << "\nProbability of mutation of real variable = " << fPMut
              << std::endl;
    std::cout << "\nDistribution index for crossover = " << fEtaCross
              << std::endl;
    std::cout << "\nDistribution index for mutation = " << fEtaMut << "\n"
              << std::endl;
  }
  if (popfunction) {
    fParentPop->SetPopFunction(popfunction);
    fChildPop->SetPopFunction(popfunction);
    fMixedPop->SetPopFunction(popfunction);
  }
  fGen = 1;
  std::cout<< "-==============================================-"<<std::endl;
  std::cout << "New generetion #" << fGen << std::endl;
  fParentPop->Build();
  fParentPop->Evaluate();
  fParentPop->FastNonDominantSorting();
  fParentPop->CrowdingDistanceAll();
  }

void AlgorithmNSGA::Selection(Population<Double_t> &oldpop,
                              Population<Double_t> &newpop) throw(ExceptionMessenger) {
  static TRandom3 rand;
  const Int_t PopSizeCheck = oldpop.GetPopulationSize();
  if((newpop.GetPopulationSize()) != PopSizeCheck)
    throw ExceptionMessenger("OMG! New population has wrong size");
  std::vector<Int_t> VecIndexGenes1(PopSizeCheck), VecIndexGenes2(PopSizeCheck);
  for (Int_t i = 0; i < PopSizeCheck; ++i) {
    VecIndexGenes1[i] = VecIndexGenes2[i] = i;
  }
  for (Int_t i = 0; i < PopSizeCheck; ++i) {
    std::swap(VecIndexGenes1[rand.Uniform(i, PopSizeCheck - 1)], VecIndexGenes1[i]);
    std::swap(VecIndexGenes2[rand.Uniform(i, PopSizeCheck - 1)], VecIndexGenes2[i]);
  }
  for (Int_t i = 0; i < PopSizeCheck; i += 4) {
    Genes<Double_t> &Combination11 =
        Tournament(oldpop.GetGenes(VecIndexGenes1[i]), oldpop.GetGenes(VecIndexGenes1[i + 1]));
    Genes<Double_t> &Combination12 =
        Tournament(oldpop.GetGenes(VecIndexGenes1[i + 2]), oldpop.GetGenes(VecIndexGenes1[i + 3]));
    Crossover(Combination11, Combination12, newpop.GetGenes(i), newpop.GetGenes(i + 1));
    Genes<Double_t> &Combination21 =
        Tournament(oldpop.GetGenes(VecIndexGenes2[i]), oldpop.GetGenes(VecIndexGenes2[i + 1]));
    Genes<Double_t> &Combination22 =
        Tournament(oldpop.GetGenes(VecIndexGenes2[i + 2]), oldpop.GetGenes(VecIndexGenes2[i + 3]));
    Crossover(Combination21, Combination22, newpop.GetGenes(i + 2), newpop.GetGenes(i + 3));
  }
} 

Genes<Double_t> &AlgorithmNSGA::Tournament(Genes<Double_t> &ind1,
                                           Genes<Double_t> &ind2) const {
  static TRandom rnd;
  const Functions* setupind = ind1.GetSetup();
  std::cout << "So so - number of objectives " << ind1.GetSetup() << std::endl;
  Int_t fFlag = ind1.CheckDominance(const_cast<Functions* >(setupind), &ind2);
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
  Double_t Component1, Component2, LimitDown, LimitUp;
  Double_t Crossover1, Crossover2;
  Double_t Alpha, Beta, Betaq;
  Int_t Rand = rnd.Rndm();
  if (Rand <= GetPCross()) {
    fNCross++;
    for (Int_t i = 0; i < fNParam; i++) {
      if (fabs(parent1[i] - parent2[i]) > EPS) {
        if (parent1[i] < parent2[i]) {
          Component1 = parent1[i];
          Component2 = parent2[i];
        } else {
          Component1 = parent2[i];
          Component2 = parent1[i];
        }
        LimitDown = fInterval[i].first;
        LimitUp = fInterval[i].second;
        Int_t Rand = rnd.Rndm();
        Beta = 1.0 + (2.0 * (Component1 - LimitDown) / (Component2 - Component1));
        Alpha = 2.0 - pow(Beta, -(GetEtaCross() + 1.0));
        if (Rand <= (1.0 / Alpha)) { // This is a contracting crossover
          Betaq = pow((Rand * Alpha), (1.0 / (GetEtaCross() + 1.0)));
        } else { // This is an expanding crossover
          Betaq = pow((1.0 / (2.0 - Rand * Alpha)), (1.0 / (GetEtaCross() + 1.0)));
        }
        Crossover1 = 0.5 * ((Component1 + Component2) - Betaq * (Component2 - Component1));
        Beta = 1.0 + (2.0 * (LimitUp - Component2) / (Component2 - Component1));
        Alpha = 2.0 - pow(Beta, -(GetEtaCross() + 1.0));
        if (Rand <= (1.0 /Alpha)) { // This is a contracting crossover
          Betaq = pow((Rand * Alpha), (1.0 / (GetEtaCross() + 1.0)));
        } else { // This is an expanding crossover
          Betaq = pow((1.0 / (2.0 - Rand * Alpha)), (1.0 / (GetEtaCross() + 1.0)));
        }
        Crossover2 = 0.5 * ((Component1 + Component2) + Betaq * (Component2 - Component1));
        Crossover1 = fmin(fmax(Crossover1, LimitDown), LimitUp);
        Crossover2 = fmin(fmax(Crossover2, LimitDown), LimitUp);
        if (rnd.Uniform() <= 0.5) {
          child1.SetGene(i, Crossover2);
          child2.SetGene(i, Crossover1);
        } else {
          std::cout << "Print Crossover part 1 = " << Crossover1 << std::endl;
          child1.SetGene(i, Crossover1);
          std::cout << "Print Crossover part 2 = " << Crossover2 << std::endl;
          child2.SetGene(i, Crossover2);
        }
      } else {
        child1.SetGene(i, parent1[i]);
        child2.SetGene(i, parent2[i]);
      }
    }
  } else {
    for (Int_t i = 0; i < fNParam; i++) {
      child1.SetGene(i, parent1[i]);
      child2.SetGene(i, parent2[i]);
    }
  }
}

void AlgorithmNSGA::NextStep() {
  std::cout<< "-==============================================-"<<std::endl;
  std::cout << "New generetion #" << fGen + 1 << std::endl;
  Selection(*fParentPop, *fChildPop);
  fNMut = fChildPop->Mutate(); // not a std::pair (?)
  fChildPop->GenCounter = fNGen + 1;
  fChildPop->Evaluate();
  // fNMut += fNMut;
  fMixedPop->Merge(*fParentPop, *fChildPop);
  fMixedPop->GenCounter = fGen + 1;
  fMixedPop->FastNonDominantSorting();
  fParentPop->Clear();
  // until |Pt+1| + |Fi| <= N, until parent population is filled
  Int_t i = 0;
  while (fParentPop->GetPopulationSize() + (fMixedPop->GetFront(i)).size() <
         fMixedPop->GetPopulationSetupSize()) {
    Genes<Double_t> Fi = fMixedPop->GetFront(i);
    fMixedPop->CrowdingDistanceFront(i);            // calculate crowding in Fi
    for (Int_t j = 0; (Double_t)j < Fi.size(); ++j) // Pt+1 = Pt+1 U Fi
      fParentPop->fPopulation.push_back(fMixedPop->GetFront(j));
    i += 1;
  }
  fMixedPop->CrowdingDistanceFront(i); // calculate crowding in F
  std::sort(fMixedPop->GetFront(i).begin(), fMixedPop->GetFront(i).end(),
            Sort(*fMixedPop));
  const int extra =
      fParentPop->GetPopulationSetupSize() - fParentPop->GetPopulationSize();
  for (int j = 0; j < extra; ++j) // Pt+1 = Pt+1 U Fi[1:N-|Pt+1|]
    fParentPop->fPopulation.push_back(fMixedPop->GetFront(j));
  fParentPop->GenCounter = fGen + 1;
}

void AlgorithmNSGA::Evolution() {
  while (fGen <= GetGenTotalNumber()) {
    NextStep();
  } // Check through population object
}
