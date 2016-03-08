#include <cmath>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <vector>    // std::vector
#include <algorithm> // std::copy

//
//
//
// Constrained version...
//
//
//


#ifdef NUMERIC_LIB

#include "Population.h"
#include "Functions.h"
#include "HistogramManager.h"
#include "TGenes.h"
#include "AlgorithmNSGA.h"
#include "GeantVFitness.h"

#include <boost/math/constants/constants.hpp>

#include "Rtypes.h"

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

static const double pi = boost::math::constants::pi<Double_t>();

void DTLZ7(Genes<Double_t> &individual) {
  Int_t n = individual.GetSetup()->GetNParam();      // 30
  Int_t m = individual.GetSetup()->GetNObjectives(); // 3
  //
  Int_t k = n - m + 1;                               // 28
  //
  Double_t g = 0.0;
  //
  for (Int_t i = 0; i < m - 1; ++i) {
    individual.SetFitness(i, individual.GetGene(i));
  }
  //
  for(Int_t i = m - 1; i < n; ++i){
    g += individual.GetGene(i);
  }
  //
  g = 1 + 9*g/k;
  //
  Double_t h = m;
  //
  for (Int_t j = 0; j < m -1; ++j) {
    h -=individual.GetGene(j)/(1 + g)*(1 + sin(3*pi*individual.GetGene(j)));
  }
  individual.SetFitness(m -1, (1 + g)*h);
  return;
}

int main(int argc, char *argv[]) {
  // Function definition
  Functions *geantv = new Functions();
  // geantv->SetInterval(); // don't work because we initialize fNparam after...
  // STUPID SOLUTION // change on .emplace()
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));
  geantv->fInterval.push_back(std::make_pair(0, 1));

  // geantv->fInterval.resize(7);
  // geantv->std::fill(fInterval.begin(),fInterval.begin(),std::make_pair(0,
  // 1));

  std::cout << "-==============================================-" << std::endl;
  geantv->PrintLimit(geantv->fInterval);
  std::cout << "-==============================================-" << std::endl;
  // Algorithm  definition
  AlgorithmNSGA *nsga2 = new AlgorithmNSGA();
  nsga2->SetPCross(1.0);
  nsga2->SetPMut(0.033333);
  nsga2->SetGenTotalNumber(500);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(30);
  nsga2->SetNObjectives(3); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(100);
  nsga2->SetEtaMut(20);
  nsga2->SetEtaCross(15);
  nsga2->SetEpsilonC(0.01);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&DTLZ7);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
#else
int main(int argc, char *argv[]) {
  std::cout << "Numeric based test: disable Geant-V in cmake flags and re-run "
               "compilation"
            << std::endl;
}
#endif
