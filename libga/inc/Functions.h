#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include "TGenes.h"
#include "invoke_cpp11.hpp"

#include <vector>
#include <limits>
#include <functional>

template<class T> class Genes;
template<class T> class Population;

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
  Functions() : fNParam(0), fInterval(), fNCons(0), fNObjectives(0) {
    fgFunction = this;
  }

  /**
   * @brief Simple constructor with known number of parameters to be observed
   */
  Functions(Int_t nparam)
      : fNParam(nparam), fInterval(), fNCons(0), fNObjectives(0){
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

  ///////// Parameter definition //////////////
  Int_t GetNParam() const { return fNParam; }
  void SetNParam(Int_t nparam) { nparam = fNParam; }

  //////// Function definition ///////////////
  // Ugly instantiation
  //void SetFunctionGenes(void (*Function)(Genes<Double_t> &),Genes<Double_t> &ind);

  //////// Constrains definition //////////////
  Int_t GetNCons() const { return fNCons; }
  void SetNCons(Int_t ncon) { ncon = fNCons; }
  /////////////////////////////////////////////
  Int_t GetNObjectives() const { return fNObjectives; }
  void SetNObjectives(Int_t nobj) { nobj = fNObjectives; }
  /////////////////////////////////////////////
  static Functions *Instance();

private:
  Int_t fNParam; // Number of parameters
  mutable std::vector<std::pair<Double_t, Double_t>> fInterval; // Interval
                                                                // settings for
                                                                // genes in
                                                                // cromosome ->
                                                                // inheritance (?)
                                                                // from
                                                                // function
  Int_t fNCons;          
  Int_t fNObjectives; // Number of fitness values (objectives)
  static Functions *fgFunction;

public:
    //typedef std::function<void(Genes<Double_t> &)> fFunction;
  typedef void (*functype)(Genes<Double_t> *); // still dont know
  typedef void (*popfunctype)(Population<Double_t>&); // still dont know
  functype evfunc;

  ClassDef(Functions, 1)
};

#endif
