#include <vector>
#include <algorithm>
#include <cassert>
#include <random>

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

    struct Sorting {
  const Population<Double_t> &pop;
  Sorting(Population<Double_t> &population) : pop(population){};
  bool operator()(Int_t i, Int_t j) {
    const Genes<Double_t> &ind1 =
        const_cast<Population<Double_t> &>(pop).GetGenes(i);
    const Genes<Double_t> &ind2 =
        const_cast<Population<Double_t> &>(pop).GetGenes(j);
    if (ind1.GetRank() < ind2.GetRank())
      return true;
    else if (ind1.GetRank() == ind2.GetRank() &&
             ind1.GetCrowdingDistance() > ind2.GetCrowdingDistance())
      return true;
    return false;
  };
};

#ifdef ENABLE_GEANTV
AlgorithmNSGA::AlgorithmNSGA()
    : function(0), popfunction(0), fPCross(0), fEtaCross(0), fNCross(0),
      fNMut(0), fNGen(0), fParentPop(0), fChildPop(0), fMixedPop(0),
      fSizePop(0), fNParam(0), fInterval(0), fNCons(0), fNObjectives(0),
      fPMut(0), fEtaMut(0), fEpsilonC(0), fCrowdingObj(true) /*,fProp(0)*/ {}
#else
AlgorithmNSGA::AlgorithmNSGA()
    : function(0), popfunction(0), fPCross(0), fEtaCross(0), fNCross(0),
      fNMut(0), fNGen(0), fParentPop(0), fChildPop(0), fMixedPop(0),
      fSizePop(0), fNParam(0), fInterval(0), fNCons(0), fNObjectives(0),
      fPMut(0), fEtaMut(0), fEpsilonC(0), fCrowdingObj(true) {}
#endif

AlgorithmNSGA::~AlgorithmNSGA() {
  if (fChildPop) {
    //delete fChildPop;
    fChildPop = 0;
  }
  if (fParentPop) {
    //delete fParentPop;
    fParentPop = 0;
  }
  if (fMixedPop) {
    //delete fMixedPop;
    fMixedPop = 0;
  }
}

void AlgorithmNSGA::Initialize() throw(ExceptionMessenger) {

  std::cout << "INITIALIZATION - Let's check NSGA2 configuration..."
            << std::endl;

  if (fSizePop < 4 && fSizePop % 4 != 0)
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
      new std::shared_ptr<Population<Double_t>>(fSizePop, fNParam, fNCons, fNObjectives,
                               fEpsilonC, fPMut, fEtaMut, fInterval, function);
  fParentPop =
      new std::shared_ptr<Population<Double_t>>(fSizePop, fNParam, fNCons, fNObjectives,
                               fEpsilonC, fPMut, fEtaMut, fInterval, function);
  fMixedPop =
      new std::shared_ptr<Population<Double_t>>(fSizePop * 2, fNParam, fNCons, fNObjectives,
                               fEpsilonC, fPMut, fEtaMut, fInterval, function);
  // Missing check of input variables and creation of population with them
  // Report(configuration);
  std::cout << "-==============================================-" << std::endl;
  std::cout << "Population size = " << fSizePop
            << "\nNumber of generations = " << fNGen
            << "\nNumber of objective functions = " << fNObjectives
            << "\nNumber of constraints = " << fNCons
            << "\nNumber of variables = " << fNParam;

  if (fNParam != 0) {
    for (int i = 0; i < fNParam; ++i) {
      std::cout << "\nLower limit of real variable " << (i + 1) << " = "
                << fInterval[i].first;
      std::cout << "\nUpper limit of real variable " << (i + 1) << " = "
                << fInterval[i].second;
    }
    std::cout << "\nProbability of crossover of real variable = " << fPCross;
    std::cout << "\nProbability of mutation of real variable = " << fPMut;
    std::cout << "\nDistribution index for crossover = " << fEtaCross;
    std::cout << "\nDistribution index for mutation = " << fEtaMut << "\n"
              << std::endl;
  }
  if (popfunction) {
    fParentPop->SetPopFunction(popfunction);
    fChildPop->SetPopFunction(popfunction);
    fMixedPop->SetPopFunction(popfunction);
  }

  fParentPop->fCrowdingObj = fCrowdingObj;
  fMixedPop->fCrowdingObj = fCrowdingObj;
  fChildPop->fCrowdingObj = fCrowdingObj;
  fGen = 1;
  std::cout << "=================New generation #" << fGen
            << "================" << std::endl;
  fParentPop->Build();
  //#ifdef ENABLE_GEANTV
  //  fParentPop->Evaluate(fProp);
  //#else
  fParentPop->Evaluate();
  //#endif
  std::cout << "Calculation of FastNDS for parent population" << std::endl;
  fParentPop->FastNonDominantSorting();
  std::cout << "-==============================================-" << std::endl;
  std::cout << "Calculation of CD for parent population" << std::endl;
  fParentPop->CrowdingDistanceAll();
  std::cout << "-==============================================-" << std::endl;
}

void AlgorithmNSGA::Selection(
    Population<Double_t> &oldpop,
    Population<Double_t> &newpop) throw(ExceptionMessenger) {
  std::random_device rd;
  std::mt19937 gen(rd());
  const Int_t PopSizeCheck = oldpop.GetPopulationSize();
  if ((newpop.GetPopulationSize()) != PopSizeCheck)
    throw ExceptionMessenger("OMG! New population has wrong size");
  std::vector<Int_t> Population1(PopSizeCheck), Population2(PopSizeCheck);
  for (Int_t i = 0; i < PopSizeCheck; ++i) {
    Population1[i] = Population2[i] = i;
  }
  for (Int_t i = 0; i < PopSizeCheck; ++i) {
    std::uniform_int_distribution<> dist(i, PopSizeCheck - 1);
    // std::swap(Population1[rand.Integer(PopSizeCheck - 1)], Population1[i]);
    // std::swap(Population2[rand.Integer(PopSizeCheck - 1)], Population2[i]);
    // Int_t rand1 = rand() % (PopSizeCheck - 1) + i;
    // std::swap(Population1[rand1], Population1[i]);
    // Int_t rand2 = rand() % (PopSizeCheck - 1) + i;
    // std::swap(Population2[rand2], Population2[i]);
    std::swap(Population1[dist(gen)], Population1[i]);
    std::swap(Population2[dist(gen)], Population2[i]);
  }
  for (Int_t i = 0; i < PopSizeCheck; i += 4) {
    Genes<Double_t> &Combination11 = Tournament(
        oldpop.GetGenes(Population1[i]), oldpop.GetGenes(Population1[i + 1]));
    Genes<Double_t> &Combination12 =
        Tournament(oldpop.GetGenes(Population1[i + 2]),
                   oldpop.GetGenes(Population1[i + 3]));
    // std::cout << "-================== Combination11 & Combination12
    // ============================-" << std::endl;
    // Combination11.Genes<Double_t>::printGenes(Combination11);
    // std::cout << "-==============================================-" <<
    // std::endl;
    // Combination12.Genes<Double_t>::printGenes(Combination12);
    Crossover(Combination11, Combination12, newpop.GetGenes(i),
              newpop.GetGenes(i + 1));
    Genes<Double_t> &Combination21 = Tournament(
        oldpop.GetGenes(Population2[i]), oldpop.GetGenes(Population2[i + 1]));
    Genes<Double_t> &Combination22 =
        Tournament(oldpop.GetGenes(Population2[i + 2]),
                   oldpop.GetGenes(Population2[i + 3]));
    // std::cout << "-================== Combination21 & Combination22
    // ============================-" << std::endl;
    // Combination21.Genes<Double_t>::printGenes(Combination21);
    // std::cout << "-==============================================-" <<
    // std::endl;
    // Combination22.Genes<Double_t>::printGenes(Combination22);
    // std::cout << "-==============================================-" <<
    // std::endl;
    Crossover(Combination21, Combination22, newpop.GetGenes(i + 2),
              newpop.GetGenes(i + 3));
  }
}

Genes<Double_t> &AlgorithmNSGA::Tournament(Genes<Double_t> &ind1,
                                           Genes<Double_t> &ind2) const {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0, 1);
  const std::shared_ptr<Functions> setupind = ind1.GetSetup();
  // std::cout << "So so, just to be sure - number of objectives on Tournament()
  // " << ind1.GetSetup()->GetNObjectives() << std::endl;
  Int_t fFlag = ind1.CheckDominance(const_cast<std::shared_ptr<Functions>>(setupind), &ind2);
  if (fFlag == 1) // Yes
    return ind1;
  else if (fFlag == -1) // Opposite
    return ind2;
  else if (ind1.GetCrowdingDistance() > ind2.GetCrowdingDistance())
    return ind1;
  else if (ind2.GetCrowdingDistance() > ind1.GetCrowdingDistance())
    return ind2;
  else if (dist(gen) <= 0.5)
    return ind1;
  else
    return ind2;
}

void AlgorithmNSGA::Crossover(const Genes<Double_t> &parent1,
                              const Genes<Double_t> &parent2,
                              Genes<Double_t> &child1,
                              Genes<Double_t> &child2) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0, 1);
  Double_t Component1, Component2, LimitDown, LimitUp;
  Double_t Crossover1, Crossover2;
  Double32_t Alpha, Beta, Betaq;

  Double_t Rand = dist(gen);
  if (Rand <= GetPCross()) {
    fNCross++;
    for (Int_t i = 0; i < fNParam; ++i) {
      if (fabs(parent1.GetGene(i) - parent2.GetGene(i)) > EPS) {
        if (parent1.GetGene(i) < parent2.GetGene(i)) {
          Component1 = parent1.GetGene(i);
          //std::cout << "Component1 = " << Component1 << std::endl;
          Component2 = parent2.GetGene(i);
          //std::cout << "Component2 = " << Component2 << std::endl;
        } else {
          Component1 = parent2.GetGene(i);
          //std::cout << "Component1 = " << Component1 << std::endl;
          Component2 = parent1.GetGene(i);
          //std::cout << "Component2 = " << Component2 << std::endl;
        }
        LimitDown = fInterval[i].first;
        LimitUp = fInterval[i].second;
        //std::cout << "-======================= Preparations before Crossover1=======================-"<< std::endl;
        Beta =
            1.0 + (2.0 * (Component1 - LimitDown) / (Component2 - Component1));
        //std::cout << "Beta = " << Beta << std::endl;
        Alpha = 2.0 - pow(Beta, -(GetEtaCross() + 1.0));
        //std::cout << "Alpha = " << Alpha << std::endl;
        if (Rand <= (1 / Alpha)) { // This is a contracting crossover
          Betaq = pow((Rand * Alpha), (1.0 / (GetEtaCross() + 1.0)));
          //std::cout << "Betaq = " << Betaq << std::endl;
        } else { // This is an expanding crossover
          Betaq =
              pow((1.0 / (2.0 - Rand * Alpha)), (1.0 / (GetEtaCross() + 1.0)));
          //std::cout << "Betaq = " << Betaq << std::endl;
        }
        Crossover1 = 0.5 * ((Component1 + Component2) -
                            Betaq * (Component2 - Component1));
        //std::cout << "-======================= Preparations before Crossover2=======================-"<< std::endl;
        Beta = 1.0 + (2.0 * (LimitUp - Component2) / (Component2 - Component1));
        //std::cout << "Beta = " << Beta << std::endl;
        Alpha = 2.0 - pow(Beta, -(GetEtaCross() + 1.0));
        //std::cout << "Alpha = " << Alpha << std::endl;
        if (Rand <= (1 / Alpha)) { // This is a contracting crossover
          Betaq = pow((Rand * Alpha), (1.0 / (GetEtaCross() + 1.0)));
          //std::cout << "Betaq = " << Betaq << std::endl;
        } else { // This is an expanding crossover
          Betaq =
              pow((1.0 / (2.0 - Rand * Alpha)), (1.0 / (GetEtaCross() + 1.0)));
          //std::cout << "Betaq = " << Betaq << std::endl;
        }
        Crossover2 = 0.5 * ((Component1 + Component2) +
                            Betaq * (Component2 - Component1));
        //////////////////////////////////////////////////////////////////////
        Crossover1 = fmin(fmax(Crossover1, LimitDown), LimitUp);
        //std::cout << "-======================= Time for a new Crossover()=======================-"<< std::endl;
        //std::cout << "DEBUG: Crossover1 = "<< Crossover1 << std::endl;
        Crossover2 = fmin(fmax(Crossover2, LimitDown), LimitUp);
        //std::cout << "DEBUG: Crossover2 = "<< Crossover2 << std::endl;
        if (Rand <= 0.5) {
          child1.SetGene(i, Crossover2);
          //std::cout << "============== DEBUG Child1: SetGene(i,Crossover2) on i-parameter place =================" << std::endl;
          //child1.Genes<Double_t>::printGenes(child1);
          child2.SetGene(i, Crossover1);
          //std::cout << "============== DEBUG Child2: SetGene(i,Crossover1) on i-parameter place =================" << std::endl;
          //child2.Genes<Double_t>::printGenes(child2);
        } else {
          child1.SetGene(i, Crossover1);
          //std::cout << "============== DEBUG Child1: SetGene(i,Crossover1) on i-parameter place =================" << std::endl;
          //child1.Genes<Double_t>::printGenes(child1);
          child2.SetGene(i, Crossover2);
          //std::cout << "============== DEBUG Child2: SetGene(i,Crossover2) on i-parameter place =================" << std::endl;
          //child2.Genes<Double_t>::printGenes(child2);
        }
      } else {
        child1.SetGene(i, parent1.GetGene(i));
        //std::cout << "============== DEBUG Child1: SetGene(i,Parent1-i-Genes) on i-parameter place ================="
        //          << std::endl;
        //child1.Genes<Double_t>::printGenes(child1);
        child2.SetGene(i, parent2.GetGene(i));
        //std::cout << "============== DEBUG Child2: SetGene(i,Parent2-i-Genes) on i-parameter place ================="
        //          << std::endl;
        //child2.Genes<Double_t>::printGenes(child2);
      }
    }
  } else {
    for (Int_t i = 0; i < fNParam; i++) {
      child1.SetGene(i, parent1.GetGene(i));
       //std::cout << "============== DEBUG Child1: SetGene(i,Parent1-i-Genes)on i-parameter place =================" << std::endl;
       //child1.Genes<Double_t>::printGenes(child1);
      child2.SetGene(i, parent2.GetGene(i));
       //std::cout << "============== DEBUG Child2: SetGene(i,Parent2-i-Genes) on i-parameter place =================" << std::endl;
       //child2.Genes<Double_t>::printGenes(child1);
    }
  }
  child1.SetEvaluated(false);
  child2.SetEvaluated(false);
}

void AlgorithmNSGA::NextStep() {
  std::cout << "-==============================================-"
            << "\n" << std::endl;
  std::cout << "New generation #" << fGen + 1 << std::endl;
  Selection(*fParentPop, *fChildPop);
  // Writing new Parent population in root file!
  //fParentPop->WritePopulationTree(*fParentPop, "NSGA.root");
  // Printing ChildPop output
  // fChildPop->Print();
  fNMut = fChildPop->Mutate();
  std::cout << "-==============================================-" << std::endl;
  std::cout << "Number of mutations in ChildPop: " << fNMut << "\n";
  fChildPop->GenCounter = fNGen + 1;
  //#ifdef ENABLE_GEANTV
  //  fChildPop->Evaluate(fProp);
  //#else
  fChildPop->Evaluate();
  //#endif
  // fNMut += fNMut;
  fMixedPop->Merge(*fParentPop, *fChildPop);
  fMixedPop->GenCounter = fGen + 1;
  std::cout << "Calculation of FastNDS for new mixed population" << std::endl;
  fMixedPop->FastNonDominantSorting();
  std::cout << "-==============================================-" << std::endl;
  fParentPop->fPopulation.clear();
  fParentPop->SetPopulationSize(0);
  Int_t i = 0;
  std::cout << "-==============================================-" << std::endl;
  for (int i = 0; i < fMixedPop->fFront.size(); ++i) {
    std::cout << CYAN << "MixedPop: fFront<" << i << "> = [";
    for (auto &it : fMixedPop->fFront[i]) {
      std::cout << it << ' ';
    }
    std::cout << "]"
              << "\n" << RESET << std::endl;
  }
  std::cout << "Next step: until |Pt+1| + |Fi| <= N" << std::endl;
  while (fParentPop->GetPopulationSize() + fMixedPop->GetFront(i).size() <
         fSizePop) {
    std::cout << "fParentPop->GetPopulationSize() = "
              << fParentPop->GetPopulationSize() << std::endl;
    std::cout << "fMixedPop->GetFront(i).size() = "
              << fMixedPop->GetFront(i).size() << std::endl;
    std::cout << "fSizePop = " << fSizePop << std::endl;
    std::cout << "-==============================================-"
              << std::endl;
    std::vector<Int_t> Fi = fMixedPop->GetFront(i);
    std::cout << "Calculation of CD for new mixed population" << std::endl;
    fMixedPop->CrowdingDistanceFront(i);
    std::cout << "-==============================================-"
              << std::endl;
    // calculate crowding in Fi
    std::cout << "Next step: Pt+1 = Pt+1 U Fi" << std::endl;
    for (Int_t j = 0; j < Fi.size(); ++j) {
      Genes<Double_t> &ind = fMixedPop->GetGenes(Fi[j]);
      // ind.printGenes(ind);
      fParentPop->push_back(ind);
    }
    i += 1;
  }
  fMixedPop->CrowdingDistanceFront(i);
  // calculate crowding in F
  std::sort(fMixedPop->GetFront(i).begin(), fMixedPop->GetFront(i).end(),
            Sorting(*fMixedPop));
  const int extra = fSizePop - (fParentPop->GetPopulationSize());
  std::cout << "Next step: Pt+1 = Pt+1 U Fi[1:N-|Pt+1|]" << std::endl;
  for (int j = 0; j < extra; ++j) {
    std::vector<Int_t> thisFront = fMixedPop->GetFront(i);
    Int_t i = thisFront[j];
    Genes<Double_t> &indadd = fMixedPop->GetGenes(i);
    fParentPop->push_back(indadd);
  }
  // std::cout << "-==============================================-" <<
  // std::endl;
  // fMixedPop->Print();
  // std::cout << "-==============================================-" <<
  // std::endl;
  // fParentPop->Print();
  // std::cout << "-==============================================-" <<
  // std::endl;
  // fChildPop->Print();
  // std::cout << "-==============================================-" <<
  // std::endl;
  fGen += 1;
  // Print new population
  fParentPop->Print();
  //Dump
  fParentPop->Store("ParentPop.backup", *fParentPop);
  //Write current populations in file
  std::ofstream populationcvs;
  populationcvs.open("GeantV.cvs", std::fstream::app);
  //fParentPop->CVSOutput(fParentPop->CreateCVS("Population.csv"), *fParentPop);
  fParentPop->CVSOutput(populationcvs, *fParentPop);
}

void AlgorithmNSGA::Evolution() {
  while (fGen < GetGenTotalNumber()) {
    std::cout << "Current generation is " << fGen << " that < "
              << "number of requested generations " << GetGenTotalNumber()
              << std::endl;
    NextStep();
  } // Check through population object
  if (fGen == GetGenTotalNumber())
    std::cout << "Finishing work.." << std::endl;
}
