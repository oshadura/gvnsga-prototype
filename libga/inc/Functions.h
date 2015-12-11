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
   * 
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
  void SetInterval(std::vector<std::pair<Double_t, Double_t>> l){ l = fInterval; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  Int_t GetNParam() const { return fNParam; }
  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param nparam [description]
   */
  void SetNParam(Int_t nparam) { nparam = fNParam; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  Int_t GetNCons() const { return fNCons; }
  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param ncon [description]
   */
  void SetNCons(Int_t ncon) { ncon = fNCons; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  Int_t GetNObjectives() const { return fNObjectives; }
  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param nobj [description]
   */
  void SetNObjectives(Int_t nobj) { nobj = fNObjectives; }
  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param etamut [description]
   */
  void SetEtaMut(Double_t etamut) { fEtaMut = etamut; }
  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param pmut [description]
   */
  void SetPMut(Double_t pmut) { fPMut = pmut; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  Double_t GetPMut() const { return fPMut; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  Double_t GetEtaMut() const { return fEtaMut; }
  /**
   * @brief [brief description]
   * @details [long description]
   * 
   * @param epsc [description]
   */
  void SetEpsilonC(Double_t epsc) { fEpsilonC = epsc; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  Double_t GetEpsilonC() const { return fEpsilonC; }
  /**
   * @brief [brief description]
   * @details [long description]
   * @return [description]
   */
  static Functions *Instance();
  /**
  friend std::ostream& operator<<(std::ostream& os, std::vector<std::pair<Double_t,Double_t>> &limit){
  for(auto &x:limit){
    os << x.first << ":"<< x.second;
  }
    return os;
  }
  */

  void PrintLimit( std::vector<std::pair<Double_t,Double_t>> &limit){
  std::cout << "Check what we create as a limit vector: [";
  for(auto &x:limit){
    std::cout << x.first << ":"<< x.second<<" ";
  }
  std::cout << "]"<< std::endl;
  }

private:
  Int_t fNParam; // Number of parameters
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
  std::vector<std::pair<Double_t, Double_t>> fInterval; // Interval
                                                        // settings for
                                                        // genes in
                                                        // cromosome ->
                                                        // inheritance (?)
                                                        // from
                                                        // function

  ClassDef(Functions, 1)
};

#endif
