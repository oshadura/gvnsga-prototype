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

#include "Rtypes.h"

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

static const double pi = boost::math::constants::pi<double>();

void DTLZ1(Genes<Double_t> &individual) {
  /*
  double y = 0.0;

  for(Int_t i = 0; i < individual.GetfGenes().size(); ++i) {
    y += pow(individual.GetGene(i) - 0.5, 2) - cos(20 * pi *
  (individual.GetGene(i) - 0.5));
  }
  individual.SetFitness(0, 100.0 * (y + individual.GetfGenes().size()));
  */

  Int_t n = individual.GetSetup()->GetNParam();      // 7
  Int_t m = individual.GetSetup()->GetNObjectives(); // 3
  Int_t k = n - m + 1;                               // 5

  Double_t g = 0.0;

  for (Int_t i = m - 1; i < n; ++i) {
    g += pow(individual.GetGene(i - 1) - 0.5, 2) -
         cos(20 * pi * (individual.GetGene(i - 1) - 0.5));
  }
  g = 100 * (k + g);

  //individual.GetFitnessVector().resize(m, 0);

  for (Int_t i = 0; i < m; ++i) {
    Double_t f = 0.5 * (1 + g);
    size_t j = 0;
    for (; m >= 2 + i && j <= m - 2 - i; ++j) {
      f *= individual.GetGene(j);
    }
    if (i > 0) {
      f *= (1 - individual.GetGene(j));
    }
    // std::cout << "Stupid iterator - tell me your value = " << i
    //          << " Function is " << f << std::endl;
    individual.SetFitness(i, f);
  }
  return;
}

int main(int argc, char *argv[]) {
  // Function definition
  Functions *geantv = new Functions();
  // geantv->SetInterval(); // don't work because we initialize fNparam after...
  // STUPID SOLUTION // cjange on .emplace()
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
  nsga2->SetPCross(1.0);
  nsga2->SetPMut(0.14285714);
  nsga2->SetGenTotalNumber(300);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(7);
  nsga2->SetNObjectives(3); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(4);
  nsga2->SetEtaMut(20);
  nsga2->SetEtaCross(15);
  nsga2->SetEpsilonC(0.7);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&DTLZ1);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
