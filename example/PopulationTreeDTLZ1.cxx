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
  double y = 0.0;
  for(Int_t i = 0; i < individual.GetfGenes().size(); ++i) {
    y += pow(individual.GetGene(i) - 0.5, 2) - cos(20 * pi * (individual.GetGene(i) - 0.5));
  }
  individual.SetFitness(0, 100.0 * (y + individual.GetfGenes().size()));
  return;
}

int main(int argc, char *argv[]) {
  // Function definition
  Functions *geantv = new Functions();
  // geantv->SetInterval(); // don't work because we initialize fNparam after...
  // STUPID SOLUTION //
  geantv->fInterval.push_back(std::make_pair(1, 10));
  geantv->fInterval.push_back(std::make_pair(1, 10));
  geantv->fInterval.push_back(std::make_pair(1, 10));
  geantv->fInterval.push_back(std::make_pair(1, 10));
  geantv->fInterval.push_back(std::make_pair(1, 10));
  geantv->fInterval.push_back(std::make_pair(1, 10));
  geantv->PrintLimit(geantv->fInterval);
  // Algorithm  definition
  AlgorithmNSGA *nsga2 = new AlgorithmNSGA();
  nsga2->SetPCross(0.5);
  nsga2->SetPMut(0.7);
  nsga2->SetGenTotalNumber(2);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(6);
  nsga2->SetNObjectives(2); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(4);
  nsga2->SetEtaMut(10);
  nsga2->SetEtaCross(10);
  nsga2->SetEpsilonC(0.7);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&DTLZ1);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
