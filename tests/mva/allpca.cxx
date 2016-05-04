#include "gtest/gtest.h"
#include "GATest.h"
#include "mva/PCA.h"
#include "mva/LPCA.h"
#include "mva/KPCA.h"

using namespace Eigen;

class AllPCA : public GATest {
public:
  LPCA lpca;
  KPCA kpca;
};

/*
TEST_F(AllPCA, UploadPopulationLPCA){
  Population<double> pop;
  lpca.UploadPopulation(pop);
  ASSERT_EQ(lpca.GetX().rows(), pop.GetPopulationSize());
  ASSERT_EQ(lpca.GetX().cols(), pop.GetSetup()->GetNParam());
}
*/

TEST_F(AllPCA, LoadingDataLPCA) {
  lpca.LoadData("data");
  MatrixXd test(4, 5);
  test(0, 0) = 2.0;
  test(1, 0) = 2.0;
  test(2, 0) = 5.0;
  test(3, 0) = 1.0;
  test(0, 1) = 1.0;
  test(1, 1) = 1.0;
  test(2, 1) = 3.0;
  test(3, 1) = 3.0;
  test(0, 2) = 1.0;
  test(1, 2) = 9.0;
  test(2, 2) = 8.0;
  test(3, 2) = 4.0;
  test(0, 3) = 3.0;
  test(1, 3) = 7.0;
  test(2, 3) = 7.0;
  test(3, 3) = 6.0;
  test(0, 4) = 1.0;
  test(1, 4) = 6.0;
  test(2, 4) = 3.0;
  test(3, 4) = 3.0;
  ////////////////
  MatrixXd Xtest;
  Xtest = lpca.GetX().transpose();
  ASSERT_EQ(Xtest, test);
}

TEST_F(AllPCA, LoadingDataKPCA) {
  lpca.LoadData("data");
  MatrixXd test(4, 5);
  test(0, 0) = 2.0;
  test(1, 0) = 2.0;
  test(2, 0) = 5.0;
  test(3, 0) = 1.0;
  test(0, 1) = 1.0;
  test(1, 1) = 1.0;
  test(2, 1) = 3.0;
  test(3, 1) = 3.0;
  test(0, 2) = 1.0;
  test(1, 2) = 9.0;
  test(2, 2) = 8.0;
  test(3, 2) = 4.0;
  test(0, 3) = 3.0;
  test(1, 3) = 7.0;
  test(2, 3) = 7.0;
  test(3, 3) = 6.0;
  test(0, 4) = 1.0;
  test(1, 4) = 6.0;
  test(2, 4) = 3.0;
  test(3, 4) = 3.0;
  ////////////////
  MatrixXd Xtest;
  Xtest = lpca.GetX().transpose();
  ASSERT_EQ(Xtest, test);
}


TEST_F(AllPCA, RunLPCA) {
	lpca.LoadData("data");
	lpca.RunLPCA();
}