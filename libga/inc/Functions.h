#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include "Genes.h"

#include <vector>
#include <limits>
#include <functional>

class Functions {

  typedef std::function<void(void)> Function;

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
      : fNParam(0), fInterval(), fConstraines(), PopulationFunction(),
        fNCons(0) {}

  /**
   * @brief Simple constructor with known number of parameters to be observed
   */
  Functions(Int_t nparam)
      : fNParam(nparam), fInterval(), fConstraines(), PopulationFunction(),
        fNCons(0) {}

  /**
   * @brief Simple constructor with known number of parameters and pre-setuped
   * pairs of vectors
   */
  Functions(Int_t nparam,
            const std::vector<std::pair<Double_t, Double_t> > limits)
      : fNParam(nparam), fInterval(limits), fConstraines(),
        PopulationFunction(), fNCons(0) {}
  /**
   * @brief Simple constructor
   */
  Functions(Int_t nparam, Int_t nconst,
            const std::vector<std::pair<Double_t, Double_t> > limits,
            const std::vector<Double_t> constr)
      : fNParam(nparam), fInterval(limits), fConstraines(constr),
        fNCons(nconst) {}

  Functions(Int_t nparam, Int_t nconst,
            const std::vector<std::pair<Double_t, Double_t> > limits,
            const std::vector<Double_t> constr,
            std::function<void(void)> fFunction)
      : fNParam(nparam), fInterval(limits), fConstraines(constr),
        PopulationFunction(fFunction), fNCons(nconst) {}

  Functions(const Functions &func);
  virtual ~Functions() {}

  ///////// Intervals definition //////////////
  std::vector<std::pair<Double_t, Double_t> > GetInterval() const {
    return fInterval;
  }
  std::pair<Double_t, Double_t> GetIntervalLimit(Int_t i) const {
    return fInterval.at(i);
  }
  Int_t GetNParam() const { return fNParam; }
  void SetInterval();
  void SetIntervalLimit(Int_t i, Double_t fMin, Double_t fMax);
  void SetFunction(void (*run)());
  void SetFunctionOpt(void (*run)(Int_t), Int_t i);

  //////// Constrains definition //////////////
  Int_t GetNCons() const { return fNCons; }
  void SetConstrain(Int_t i, Double_t value);
  std::vector<Double_t> GetConstraines() const { return fConstraines; }

  void SetNCons(Int_t ncon) { ncon = fNCons; }
  void SetNParams(Int_t nparam) { nparam = fNParam; }
  /////////////////////////////////////////////
  static Functions *Instance(); // ?????????

private:
  Int_t fNParam; // Number of parameterxs
  mutable std::vector<std::pair<Double_t, Double_t> > fInterval; // Interval
                                                                 // settings for
                                                                 // genes in
                                                                 // cromosome ->
                                                                 // inheritance
                                                                 // from
                                                                 // function
  Int_t fNBins; // Number of bins for statistics proposes
  mutable std::vector<Double_t> fConstraines; // Vector of constraines for NSGA
                                              // constrain based
  Function PopulationFunction; // type function to be passed from GeantV
  Int_t fNCons;                // Number of constrains

  ClassDef(Functions, 1)
};

#endif
