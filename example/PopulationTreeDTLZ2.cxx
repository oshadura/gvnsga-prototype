#include <cmath>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <vector>    // std::vector
#include <algorithm> // std::copy

#include "Population.h"
#include "Functions.h"
#include "HistogramManager.h"
#include "TGenes.h"
#include "AlgorithmNSGA.h"
#include "GeantVFitness.h"

#include <boost/math/constants/constants.hpp>

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

static const double halfpi = boost::math::constants::pi<Double_t>() / 2.0;

void DTLZ2(Genes<Double_t> &individual) {
  /*
  double y = 0.0;
  for (Int_t i = 0; i < individual.GetfGenes().size(); ++i) {
    y += pow(individual.GetGene(i) - 0.5, 2);
  }
  individual.SetFitness(0, y);
  */
  Int_t n = individual.GetSetup()->GetNParam(); //12
  Int_t m = individual.GetSetup()->GetNObjectives(); //3
  Int_t k = n - m + 1; //10

  Double_t g = 0.0;

  for (Int_t i = m - 1; i < n; ++i) {
    g += pow(individual.GetGene(i) - 0.5, 2);
  }

  //individual.GetFitnessVector().resize(m, 0);

  for (Int_t i = 0; i < m; ++i) {
    Double_t f = (1 + g);size_t j = 0;
    for (; i + m <=  m - 2; ++j){
      f *= cos(individual.GetGene(j) * halfpi);
    }
    if (m > 1) {
      f *=
          sin(individual.GetGene(j) * halfpi);
    }
    individual.SetFitness(i, f);
  }
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
  std::cout << "-==============================================-" << std::endl;
  geantv->PrintLimit(geantv->fInterval);
  std::cout << "-==============================================-" << std::endl;
  // Algorithm  definition
  AlgorithmNSGA *nsga2 = new AlgorithmNSGA();
  nsga2->SetPCross(0.5);
  nsga2->SetPMut(0.7);
  nsga2->SetGenTotalNumber(1);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(12);
  nsga2->SetNObjectives(3); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(4);
  nsga2->SetEtaMut(10);
  nsga2->SetEtaCross(10);
  nsga2->SetEpsilonC(0.7);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&DTLZ2);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
