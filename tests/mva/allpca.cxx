#include "GATest.h"
#include "generic/Population.h"
#include "mva/KPCA.h"
#include "mva/LPCA.h"
#include "mva/LPCAWhite.h"
#include "mva/PCA.h"
#include "mva/RobustPCA.h"
#include "mva/UncenteredLPCA.h"
#include "mva/UncenteredTrickLPCA.h"
#include "mva/UncenteredWhiteLPCA.h"
#include "mva/UncenteredTrickSVD.h"
#include "output/CSVManager.h"
#include "problem/DTLZ1.h"
#include "gtest/gtest.h"

using namespace Eigen;

class AllPCA : public GATest {
public:
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::LPCA lpca;
  geantvmoop::LPCAWhite lpcawhite;
  geantvmoop::KPCA kpca;
  geantvmoop::RobustPCA robustpca;
  geantvmoop::UncenteredLPCA ulpca;
  geantvmoop::UncenteredWhiteLPCA uwlpca;
  geantvmoop::UncenteredTrickLPCA twlpca;
  geantvmoop::UncenteredTrickSVD svdpca;
  // Pointless allocation during tests
  geantvmoop::Population<geantvmoop::DTLZ1> pop{5};
};




TEST_F(AllPCA, SVD) {
  lpca.LoadData("data");
  MatrixXd X = lpca.GetX();
  std::cout << "Print X\n"<< X << std::endl;
  MatrixXd T = (X.adjoint() * X);
  std::cout << "Print T\n" << T << std::endl;
  JacobiSVD<MatrixXd> svdX(lpca.GetX(), ComputeThinU | ComputeThinV);
    auto sX = svdX.singularValues();
    std::cout << "X Its singular values are:" << std::endl
              << svdX.singularValues() << std::endl;
    auto lvecX = svdX.matrixU();
    std::cout
        << "X Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdX.matrixU() << std::endl;
    auto rvecX = svdX.matrixV();
    std::cout
        << "X Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdX.matrixV() << std::endl;
    ///////////////////////////////////////////////////////////////////////////////
    //Checking again...
    ///////////////////////////////////////////////////////////////////////////////
    MatrixXd diag = sX.asDiagonal();
    MatrixXd recheck = lvecX*diag*rvecX.transpose();
    MatrixXd recheckGramm = rvecX*diag*diag*rvecX.transpose();
    MatrixXd UGramm = lvecX.adjoint()*lvecX;
    MatrixXd Ssquare = diag*diag;
    MatrixXd check = lvecX.adjoint()*X*X.adjoint()*lvecX;
    std::cout << "Ssquare\n" << Ssquare << std::endl;
    std::cout << "Right check formula\n" << check << std::endl;
    std::cout << "Gramm matrix\n " << recheckGramm << std::endl;
    std::cout << "UtransposeU\n" << UGramm << std::endl;
    MatrixXd  diff = recheck - lpca.GetX();
    std::cout << "diff:\n" << diff << "\n";
    std::cout << "Print this thing: X = USV*:\n" << recheck << std::endl; 
}

TEST_F(AllPCA, LoadPopulationFromCSV) {
  // geantvmoop::Population<geantvmoop::DTLZ1> population;
  // geantvmoop::CSVManager::GetInstance().LoadCSV("datasimple", population);
  // ASSERT_EQ(population.size(), 5);
}

TEST_F(AllPCA, PCALoadDataCSV) {
  lpca.LoadData("data");
  ASSERT_EQ(lpca.GetX().rows(), pop.size());
  ASSERT_EQ(lpca.GetX().cols(), 4);
}

TEST_F(AllPCA, PCAConvertXtoPopulation) {
  lpca.LoadData("dataDTLZ1");
  MatrixXd currentX = lpca.GetX();
  geantvmoop::Population<geantvmoop::DTLZ1> population;
  lpca.UnloadPopulation(population, currentX);
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

////////////////////// Methods ///////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


TEST_F(AllPCA, RunLPCA) {
  lpca.LoadData("data");
  lpca.RunLPCA();
  lpca.Print();
  lpca.WriteTransformed("outputlpca");
}

TEST_F(AllPCA, RunLPCAReductionofComponents) {
  lpca.LoadData("data");
  lpca.RunLPCAWithReductionOfComponents();
  lpca.WriteTransformed("outputlpcatransform");
}

TEST_F(AllPCA, RunLPCAReductionofComponentsNoScale) {
  lpca.LoadData("data");
  lpca.RunLPCAWithReductionOfComponentsNoScale();
  lpca.WriteTransformed("outputlpcatransformnoscale");
}

TEST_F(AllPCA, RunLPCAWhiteReductionofComponents) {
  lpcawhite.LoadData("data");
  lpcawhite.RunLPCAWhiteWithReductionOfComponents();
  lpcawhite.WriteTransformed("outputlpcawtransform");
}

TEST_F(AllPCA, RunLPCAWhiteReductionofComponentsNoScale) {
  lpcawhite.LoadData("data");
  lpcawhite.RunLPCAWhiteWithReductionOfComponentsNoScale();
  lpcawhite.WriteTransformed("outputlpcawtransformnoscale");
}
TEST_F(AllPCA, RunUncenteredSVD) {
  svdpca.LoadData("datatest");
  //svdpca.RunUncenteredSVD();
  svdpca.RunUncenteredTrickSVDWithReductionOfComponents();
  //svdpca.WriteTransformed("outputsvdpca");
}
TEST_F(AllPCA, RunUncenteredLPCA) {
  ulpca.LoadData("datatest");
  ulpca.RunUncenteredLPCA();
  ulpca.Print();
  ulpca.WriteTransformed("outputulpca");
}

TEST_F(AllPCA, RunUncenteredLPCAReductionofComponents) {
  ulpca.LoadData("data");
  ulpca.RunUncenteredLPCAWithReductionOfComponents();
  ulpca.WriteTransformed("outputulpcatransform");
}


TEST_F(AllPCA, RunUncenteredLPCAScale) {
  ulpca.LoadData("data");
  ulpca.RunUncenteredLPCAWithReductionOfComponentsScale();
  ulpca.Print();
  ulpca.WriteTransformed("outputtlpca");
}

TEST_F(AllPCA, RunUncenteredWhiteLPCA) {
  uwlpca.LoadData("data");
  uwlpca.RunUncenteredWhiteLPCA();
  uwlpca.Print();
  uwlpca.WriteTransformed("outputulpca");
}

TEST_F(AllPCA, RunUncenteredWhiteLPCAReductionofComponents) {
  uwlpca.LoadData("data");
  uwlpca.RunUncenteredWhiteLPCAWithReductionOfComponents();
  uwlpca.WriteTransformed("outputulpcatransform");
}

TEST_F(AllPCA, RunUncenteredTrickLPCA) {
  twlpca.LoadData("data");
  twlpca.RunUncenteredTrickLPCA();
  twlpca.Print();
  twlpca.WriteTransformed("outputtlpca");
}

TEST_F(AllPCA, RunUncenteredTrickLPCAReductionofComponents) {
  twlpca.LoadData("data");
  twlpca.RunUncenteredTrickLPCAWithReductionOfComponents();
  twlpca.WriteTransformed("outputtlpcatransform");
}

// Too slow for bigdata set...
TEST_F(AllPCA, RunKPCA) {
  kpca.LoadData("data");
  kpca.RunKPCA();
  kpca.Print();
  kpca.WriteTransformed("outputkpca");
}

TEST_F(AllPCA, RunRobustPCA){
  robustpca.LoadData("data");
  robustpca.RobustPCAInexact();
  robustpca.Print();
  //robustpca.WriteTransformed("outputrobustpca");
}
