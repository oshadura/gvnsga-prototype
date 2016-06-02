#include "gtest/gtest.h"
#include "GATest.h"
#include "generic/Population.h"
#include "problem/DTLZ1.h"
#include "output/CSVManager.h"
#include "mva/PCA.h"
#include "mva/LPCA.h"
#include "mva/KPCA.h"

using namespace Eigen;

class AllPCA : public GATest {
public:
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::LPCA lpca;
  geantvmoop::KPCA kpca;
  geantvmoop::Population<geantvmoop::DTLZ1> pop{ 5 };
};

TEST_F(AllPCA, PCAConvertPopulationtoX) {
  lpca.UploadPopulation(pop);
  ASSERT_EQ(lpca.GetX().rows(), 5);
  ASSERT_EQ(lpca.GetX().cols(), 7);
}

TEST_F(AllPCA, LoadPopulationFromCSV) {
  geantvmoop::Population<geantvmoop::DTLZ1> population;
  geantvmoop::CSVManager::GetInstance().LoadCSV("datasimple", population);
  ASSERT_EQ(population.size(), 5);
}

TEST_F(AllPCA, PCALoadDataCSV) {
  lpca.LoadData("data");
  ASSERT_EQ(lpca.GetX().rows(), pop.size());
  ASSERT_EQ(lpca.GetX().cols(), 3);
}

TEST_F(AllPCA, PCAConvertXtoPopulation) {
  lpca.LoadData("dataDTLZ1");
  MatrixXd currentX = lpca.GetX();
  //std::cout << "Number of cols " << currentX.cols() << std::endl;
  //std::cout << "Number of rows " << currentX.rows() << std::endl;
  //std::cout << currentX << std::endl;
  geantvmoop::Population<geantvmoop::DTLZ1> population;
  lpca.UnloadPopulation(population, currentX);
  //std::cout << population << std::endl;
  ASSERT_EQ(lpca.GetX().rows(), population.size());
  ASSERT_EQ(lpca.GetX().cols(), 7);
}

TEST_F(AllPCA, LoadingDataLPCAByHands) {
  lpca.LoadData("data");
  MatrixXd test(3, 5);
  test(0, 0) = 2.0;
  test(1, 0) = 2.0;
  test(2, 0) = 5.0;
  test(0, 1) = 1.0;
  test(1, 1) = 1.0;
  test(2, 1) = 3.0;
  test(0, 2) = 1.0;
  test(1, 2) = 9.0;
  test(2, 2) = 8.0;
  test(0, 3) = 3.0;
  test(1, 3) = 7.0;
  test(2, 3) = 7.0;
  test(0, 4) = 1.0;
  test(1, 4) = 6.0;
  test(2, 4) = 3.0;
  MatrixXd Xtest;
  Xtest = lpca.GetX().transpose();
  ASSERT_EQ(Xtest, test);
}

TEST_F(AllPCA, RunLPCA) {
  lpca.LoadData("data.example");
  lpca.RunLPCA();
  lpca.Print();
  lpca.WriteTransformed("outputlpca");
}

TEST_F(AllPCA, RunLPCAReductionofComponents) {
  lpca.LoadData("data.example");
  lpca.RunLPCAWithReductionOfComponents();
  lpca.WriteTransformed("outputlpcatransform");
}

TEST_F(AllPCA, RunKPCA) {
  /*
  kpca.LoadData("data");
  kpca.RunKPCA();
  kpca.Print();
  kpca.WriteTransformed("outputkpca");
  */
}
