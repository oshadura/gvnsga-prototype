#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include "Genes.h"

#include <vector>
#include <limits>
#include <functional>

class Genes;

class Functions {

  typedef std::function<void(Genes &)> Function;
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
      : fNParam(0), fInterval(), fConstraines(),
        fNCons(0) {
          fgFunction = this;
        }

  /**
   * @brief Simple constructor with known number of parameters to be observed
   */
  Functions(Int_t nparam)
      : fNParam(nparam), fInterval(), fConstraines(),
        fNCons(0) {
          fgFunction = this;
        }
  /**
  Functions(Int_t nparam,
            const std::vector<std::pair<Double_t, Double_t>> limits)
      : fNParam(nparam), fInterval(limits), fConstraines(),
        PopulationFunction(), fNCons(0) {}
  
  Functions(Int_t nparam, Int_t nconst,
            const std::vector<std::pair<Double_t, Double_t>> limits,
            const std::vector<Double_t> constr)
      : fNParam(nparam), fInterval(limits), fConstraines(constr),
        fNCons(nconst), PopulationFunction() {}

  Functions(Int_t nparam, Int_t nconst,
            const std::vector<std::pair<Double_t, Double_t>> limits,
            const std::vector<Double_t> constr,
            std::function<void(Genes &)> fFunction)
      : fNParam(nparam), fInterval(limits), fConstraines(constr),
        PopulationFunction(fFunction), fNCons(nconst) {}
  */

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
  virtual ~Functions() {
    fgFunction = 0;
  }

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
  void SetNParams(Int_t nparam) { nparam = fNParam; }

  //////// Function definition ///////////////
  void SetFunction(void (*fFunction)());
  void SetFunctionGenes(void (*fFunction)(Genes &), Genes &Individual);

  //////// Constrains definition //////////////
  Int_t GetNCons() const { return fNCons; }
  void SetConstrain(Int_t i, Double_t value);
  std::vector<Double_t> GetConstraines() const { return fConstraines; }
  void SetNCons(Int_t ncon) { ncon = fNCons; }

  /////////////////////////////////////////////
  static Functions *Instance();

private:
  Int_t fNParam; // Number of parameterxs
  mutable std::vector<std::pair<Double_t, Double_t>> fInterval; // Interval
                                                                // settings for
                                                                // genes in
                                                                // cromosome ->
                                                                // inheritance
                                                                // from
                                                                // function
  mutable std::vector<Double_t> fConstraines; // Vector of constraines for NSGA
                                              // constrain based
  Int_t fNCons;                // Number of constrains
  static Functions *fgFunction;


  ClassDef(Functions, 1)
};

#endif
