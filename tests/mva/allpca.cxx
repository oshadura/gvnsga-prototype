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

TEST_F(AllPCA, LoadingData) {
	lpca.LoadData("data");
	Matrix3f test;
	test << 2.0,2.0,5.0,1.0,
			1.0,1.0,3.0,3.0,
			1.0,9.0,8.0,4.0,
			3.0,7.0,7.0,6.0,
			1.0,6.0,3.0,3.0;
	EXPECT_EQ(lpca.GetX(), test);
}	