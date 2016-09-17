#include "GATest.h"
#include "generic/Population.h"
#include "mva/SpectralClustering.h"
#include "output/CSVManager.h"
#include "problem/DTLZ1.h"
#include "gtest/gtest.h"

using namespace Eigen;

class Clustering : public GATest {
public:
  Spectral spcl;
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::Population<geantvmoop::DTLZ1> pop{5};
};

TEST_F(Clustering, RunLPCAReductionofComponents) {
  spcl.read_data("data", ",");
  spcl.cluster();
}