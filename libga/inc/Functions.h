#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include "TGenes.h"
#include "invoke_cpp11.hpp"

#include <vector>
#include <limits>
#include <functional>

template<class T> class Genes;
template<class T> class Population;

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
  Functions() : fNParam(0), fInterval(), fNCons(0), fNObjectives(0),fPMut(0), fEtaMut(0),fEpsilonC(EPS) {
    fgFunction = this;
  }

  /**
   * @brief Simple constructor with known number of parameters to be observed
   */
  Functions(Int_t nparam)
      : fNParam(nparam), fInterval(), fNCons(0), fNObjectives(0),fPMut(0), fEtaMut(0), fEpsilonC(EPS){
    fgFunction = this;
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
  virtual ~Functions() { fgFunction = 0; }
  ///////// Intervals definition //////////////
  std::vector<std::pair<Double_t, Double_t>> GetInterval() const {
    return fInterval;
  }
  std::pair<Double_t, Double_t> GetIntervalLimit(Int_t i) const {
    return fInterval.at(i);
  }
  void SetInterval();
  void SetIntervalLimit(Int_t i, Double_t fMin, Double_t fMax);
  void SetInterval(std::vector<std::pair<Double_t, Double_t>> l){ l = fInterval; }
  ///////// Parameter definition //////////////
  Int_t GetNParam() const { return fNParam; }
  void SetNParam(Int_t nparam) { nparam = fNParam; }
  //////// Constrains definition //////////////
  Int_t GetNCons() const { return fNCons; }
  void SetNCons(Int_t ncon) { ncon = fNCons; }
  /////////////////////////////////////////////
  Int_t GetNObjectives() const { return fNObjectives; }
  void SetNObjectives(Int_t nobj) { nobj = fNObjectives; }
  /////////////////////////////////////////////
  void SetEtaMut(Double_t etamut) { fEtaMut = etamut; }
  void SetPMut(Double_t pmut) { fPMut = pmut; }
  Double_t GetPMut() const { return fPMut; }
  Double_t GetEtaMut() const { return fEtaMut; }
  void SetEpsilonC(Double_t epsc) { fEpsilonC = epsc; }
  Double_t GetEpsilonC() const { return fEpsilonC; }
  ////////////////////////////////////////////
  static Functions *Instance();

private:
  Int_t fNParam; // Number of parameters
  std::vector<std::pair<Double_t, Double_t>> fInterval; // Interval
                                                        // settings for
                                                        // genes in
                                                        // cromosome ->
                                                        // inheritance (?)
                                                        // from
                                                        // function
  Int_t fNCons;          
  Int_t fNObjectives; // Number of fitness values (objectives)
  Double_t fPMut;
  Double_t fEtaMut;
  Double_t fEpsilonC;
  static Functions *fgFunction;

public:
  typedef void (*functype)(Genes<Double_t> *); // still dont know
  typedef void (*popfunctype)(Population<Double_t>&); // still dont know
  //////////////////////////////////////////////////////////////
  functype evfunc;

  ClassDef(Functions, 1)
};

#endif
