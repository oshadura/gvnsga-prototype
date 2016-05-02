#ifndef __SBXCROSSOVER__
#define __SBXCROSSOVER__

#include "Crossover.h"

class Functions;

class SBXCrossover : public Crossover<SBXCrossover> {

public:
  /**
* @brief SBX Crossover
* @details TBD: add generator class
*
* @param
*
*/
  template <typename F, typename T>
  static void CrossoverGA(const Genes<T> &parent1,
                                       const Genes<T> &parent2,
                                       Genes<T> &child1, Genes<T> &child2) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    T Component1, Component2, LimitDown, LimitUp;
    T Crossover1, Crossover2;
    Double32_t Alpha, Beta, Betaq;
    T Rand = dist(gen);
    if (Rand <= GetPCross()) {
      fNCross++;
      for (Int_t i = 0; i < fNParam; ++i) {
        if (fabs(parent1.GetGene(i) - parent2.GetGene(i)) > EPS) {
          if (parent1.GetGene(i) < parent2.GetGene(i)) {
            Component1 = parent1.GetGene(i);
            // std::cout << "Component1 = " << Component1 << std::endl;
            Component2 = parent2.GetGene(i);
            // std::cout << "Component2 = " << Component2 << std::endl;
          } else {
            Component1 = parent2.GetGene(i);
            // std::cout << "Component1 = " << Component1 << std::endl;
            Component2 = parent1.GetGene(i);
            // std::cout << "Component2 = " << Component2 << std::endl;
          }
          LimitDown = fInterval[i].first;
          LimitUp = fInterval[i].second;
          // std::cout << "-======================= Preparations before
          // Crossover1=======================-"<< std::endl;
          Beta = 1.0 +
                 (2.0 * (Component1 - LimitDown) / (Component2 - Component1));
          // std::cout << "Beta = " << Beta << std::endl;
          Alpha = 2.0 - pow(Beta, -(GetEtaCross() + 1.0));
          // std::cout << "Alpha = " << Alpha << std::endl;
          if (Rand <= (1 / Alpha)) { // This is a contracting crossover
            Betaq = pow((Rand * Alpha), (1.0 / (GetEtaCross() + 1.0)));
            // std::cout << "Betaq = " << Betaq << std::endl;
          } else { // This is an expanding crossover
            Betaq = pow((1.0 / (2.0 - Rand * Alpha)),
                        (1.0 / (GetEtaCross() + 1.0)));
            // std::cout << "Betaq = " << Betaq << std::endl;
          }
          Crossover1 = 0.5 * ((Component1 + Component2) -
                              Betaq * (Component2 - Component1));
          // std::cout << "-======================= Preparations before
          // Crossover2=======================-"<< std::endl;
          Beta =
              1.0 + (2.0 * (LimitUp - Component2) / (Component2 - Component1));
          // std::cout << "Beta = " << Beta << std::endl;
          Alpha = 2.0 - pow(Beta, -(GetEtaCross() + 1.0));
          // std::cout << "Alpha = " << Alpha << std::endl;
          if (Rand <= (1 / Alpha)) { // This is a contracting crossover
            Betaq = pow((Rand * Alpha), (1.0 / (GetEtaCross() + 1.0)));
            // std::cout << "Betaq = " << Betaq << std::endl;
          } else { // This is an expanding crossover
            Betaq = pow((1.0 / (2.0 - Rand * Alpha)),
                        (1.0 / (GetEtaCross() + 1.0)));
            // std::cout << "Betaq = " << Betaq << std::endl;
          }
          Crossover2 = 0.5 * ((Component1 + Component2) +
                              Betaq * (Component2 - Component1));
          //////////////////////////////////////////////////////////////////////
          Crossover1 = fmin(fmax(Crossover1, LimitDown), LimitUp);
          // std::cout << "-======================= Time for a new
          // Crossover()=======================-"<< std::endl;
          // std::cout << "DEBUG: Crossover1 = "<< Crossover1 << std::endl;
          Crossover2 = fmin(fmax(Crossover2, LimitDown), LimitUp);
          // std::cout << "DEBUG: Crossover2 = "<< Crossover2 << std::endl;
          if (Rand <= 0.5) {
            child1.SetGene(i, Crossover2);
            // std::cout << "============== DEBUG Child1: SetGene(i,Crossover2)
            // on
            // i-parameter place =================" << std::endl;
            // child1.Genes<T> ::printGenes(child1);
            child2.SetGene(i, Crossover1);
            // std::cout << "============== DEBUG Child2: SetGene(i,Crossover1)
            // on
            // i-parameter place =================" << std::endl;
            // child2.Genes<T> ::printGenes(child2);
          } else {
            child1.SetGene(i, Crossover1);
            // std::cout << "============== DEBUG Child1: SetGene(i,Crossover1)
            // on
            // i-parameter place =================" << std::endl;
            // child1.Genes<T> ::printGenes(child1);
            child2.SetGene(i, Crossover2);
            // std::cout << "============== DEBUG Child2: SetGene(i,Crossover2)
            // on
            // i-parameter place =================" << std::endl;
            // child2.Genes<T> ::printGenes(child2);
          }
        } else {
          child1.SetGene(i, parent1.GetGene(i));
          // std::cout << "============== DEBUG Child1:
          // SetGene(i,Parent1-i-Genes)
          // on i-parameter place ================="
          //          << std::endl;
          // child1.Genes<T> ::printGenes(child1);
          child2.SetGene(i, parent2.GetGene(i));
          // std::cout << "============== DEBUG Child2:
          // SetGene(i,Parent2-i-Genes)
          // on i-parameter place ================="
          //          << std::endl;
          // child2.Genes<T> ::printGenes(child2);
        }
      }
    } else {
      for (Int_t i = 0; i < fNParam; i++) {
        child1.SetGene(i, parent1.GetGene(i));
        // std::cout << "============== DEBUG Child1:
        // SetGene(i,Parent1-i-Genes)on
        // i-parameter place =================" << std::endl;
        // child1.Genes<T> ::printGenes(child1);
        child2.SetGene(i, parent2.GetGene(i));
        // std::cout << "============== DEBUG Child2: SetGene(i,Parent2-i-Genes)
        // on i-parameter place =================" << std::endl;
        // child2.Genes<T> ::printGenes(child1);
      }
    }
    child1.SetEvaluated(false);
    child2.SetEvaluated(false);
  }
};

#endif
