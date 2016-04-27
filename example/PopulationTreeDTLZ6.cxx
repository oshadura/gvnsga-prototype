#include <cmath>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <vector>    // std::vector
#include <algorithm> // std::copy

#ifdef NUMERIC_LIB

#include "generic/Population.h"
#include "generic/Functions.h"
#include "output/HistogramManager.h"
#include "generic/TGenes.h"
#include "algorithms/AlgorithmNSGA.h"
#include "instrument/GeantVFitness.h"

#include <boost/math/constants/constants.hpp>

#include "Rtypes.h"

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

static const double halfpi = boost::math::constants::pi<Double_t>() / 2.0;

void DTLZ6(Genes<Double_t> &individual) {
  Int_t n = individual.GetSetup()->GetNParam();      // 16
  Int_t m = individual.GetSetup()->GetNObjectives(); // 3
  Int_t k = n - m + 1;                               // 20

  Double_t g = 0.0;

  for (Int_t i = m - 1; i < n; ++i) {
    g += pow(individual.GetGene(i), 0.1);
  }
  // individual.GetFitnessVector().resize(m, 0);
  std::vector<Double_t> theta(n);
  theta[0] = individual.GetGene(0) * halfpi;

  for (Int_t i = 1; i < theta.size(); ++i) {
    theta[i] = halfpi / (2 * (1 + g)) * (1 + 2 * g * individual.GetGene(i));
  }

  for (Int_t i = 0; i < m; ++i) {
    Double_t f = (1 + g);
    Int_t j = 0;
    for (; i + m <= m - 2; ++j) {
      f *= cos(theta[j]);
    }
    if (m > 0) {
      f *= sin(theta[j]);
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
  nsga2->SetPMut(0.0625);
  nsga2->SetGenTotalNumber(10);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(16);
  nsga2->SetNObjectives(3); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(100);
  nsga2->SetEtaMut(20);
  nsga2->SetEtaCross(15);
  nsga2->SetEpsilonC(0.01);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&DTLZ6);
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
