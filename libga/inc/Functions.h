#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include "TGenes.h"
#include "invoke_cpp11.hpp"

#include <vector>
#include <limits>
#include <functional>

template <class T> class Genes;
template <class T> class Population;

#define EPS 1e-14
#define INF 1e+14

class Functions {

public:
  /**
   * @brief Simple constructor
   *
   * @param fNParam number of parameters
   * @param fNCons number of constrains
   * @param fInterval vectors of itervals for parameters
   * @param fConstrains Vector of constraines
   */
  Functions()
      : fNParam(0), fInterval(), fNCons(0), fNObjectives(0), fPMut(0),
        fEtaMut(0), fEpsilonC(EPS) {
  }

  /**
   * @brief Simple constructor with known number of parameters to be observed
   */
  Functions(Int_t nparam)
      : fNParam(nparam), fInterval(), fNCons(0), fNObjectives(0), fPMut(0),
        fEtaMut(0), fEpsilonC(EPS) {
  }
  /**
   * @brief [brief description]
   * @details [long description]
   *
   * @param func [description]
   */
  Functions(const Functions &func);
  /**
   * @brief [brief description]
   * @details [long description]
   */
  virtual ~Functions() {}
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  std::vector<std::pair<Double_t, Double_t>> GetInterval() const {
    return fInterval;
  }
  /**
   * @brief [brief description]
   * @details [long description]
   *
   * @param i [description]
   * @return [description]
   */
  std::pair<Double_t, Double_t> GetIntervalLimit(Int_t i) const {
    return fInterval.at(i);
  }
  /**
   * @brief [brief description]
   * @details [long description]
   */
  void SetInterval();
  /**
   * @brief [brief description]
   * @details [long description]
   * @param i [description]
   * @param fMin [description]
   * @param fMax [description]
   */
  void SetIntervalLimit(Int_t i, Double_t fMin, Double_t fMax);
  /**
   * @brief [brief description]
   * @details [long description]
   *
   * @param l [description]
   */
  void SetInterval(std::vector<std::pair<Double_t, Double_t>> l) {
    l = fInterval;
  }
  /**
  friend std::ostream& operator<<(std::ostream& os,
  std::vector<std::pair<Double_t,Double_t>> &limit){
  for(auto &x:limit){
    os << x.first << ":"<< x.second;
  }
    return os;
  }
  */

  void PrintLimit(std::vector<std::pair<Double_t, Double_t>> &limit) {
    std::cout << "Check what we create as a limit vector: [";
    for (auto &x : limit) {
      std::cout << x.first << ":" << x.second << " ";
    }
    std::cout << "]" << std::endl;
  }

public:
  typedef void (*functype)(Genes<Double_t> *);         // still dont know
  typedef void (*popfunctype)(Population<Double_t> &); // still dont know
  //////////////////////////////////////////////////////////////
  functype evfunc;
  std::vector<std::pair<Double_t, Double_t>> fInterval; // Interval
                                                        // settings for
                                                        // genes in
                                                        // cromosome ->
                                                        // inheritance (?)
                                                        // from
                                                        // function
  Int_t fNParam; // Number of parameters
  Int_t fNCons;
  Int_t fNObjectives; // Number of fitness values (objectives)
  Double_t fPMut;
  Double_t fEtaMut;
  Double_t fEpsilonC;

  ClassDef(Functions, 1)
};

#endif
