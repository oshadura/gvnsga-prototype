#include "GATest.h"
#include "algorithms/AlgorithmNSGA.h"
#include "problem/ProblemDTLZ1.h"

#include "generic/Population.h"
#include "generic/Functions.h"
#include "output/HistogramManager.h"
#include "generic/TGenes.h"
#include "algorithms/AlgorithmNSGA.h"
#include "instrumentation/GeantVFitness.h"

#include <boost/math/constants/constants.hpp>


class Nsga2 : public GATest {
public:
	ProblemDTLZ1 dtlz1;
};

/*
TEST_F(Nsga2, SolvingProblem) {
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

  // geantv->fInterval.resize(7);
  // geantv->std::fill(fInterval.begin(),fInterval.begin(),std::make_pair(0,
  // 1));

  std::cout << "-==============================================-" << std::endl;
  geantv->PrintLimit(geantv->fInterval);
  std::cout << "-==============================================-" << std::endl;
  // Algorithm  definition
  AlgorithmNSGA *nsga2 = new AlgorithmNSGA();
  nsga2->SetPCross(1.0);
  nsga2->SetPMut(0.14285714);
  nsga2->SetGenTotalNumber(100);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(7);
  nsga2->SetNObjectives(3); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(100);
  nsga2->SetEtaMut(20);
  nsga2->SetEtaCross(15);
  nsga2->SetEpsilonC(0.01);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&dtlz1.Evaluate);
  nsga2->Initialize();
  nsga2->Evolution();
  // !!! Result of test?
}
*/

/*
TEST_F(Nsga2, TestEvaluationOfOneIndividual){
	Genes<Double_t> individual;
	individual.resize(5);
	for (int i = 0; i < 5; ++i)
	{
		individual.push_back(1);
	}
	auto evaluation = dtlz1.Evaluate(individual);
	//Need still to manage..
	//EXPECT_NEAR(1, evaluation[0], 0.01);
	//EXPECT_NEAR(1, evaluation[1], 0.01);
	//EXPECT_NEAR(1, evaluation[2], 0.01);
	//EXPECT_NEAR(1, evaluation[3], 0.01);
	//EXPECT_NEAR(1, evaluation[4], 0.01);
	ASSERT_THROW(individual.GetFitnessVector(), std::runtime_error);
}
*/
