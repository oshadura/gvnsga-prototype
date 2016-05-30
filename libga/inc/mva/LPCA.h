#ifndef __LPCA__
#define __LPCA__

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"

#include "PCA.h"

namespace geantvmoop {

using namespace Eigen;

class LPCA : public PCA<LPCA> {
public:
  LPCA() : normalise(0) {}

  explicit LPCA(MatrixXd &d) : normalise(0) { X = d; }

  virtual ~LPCA() {}

  void LoadData(const char *data, char sep = ',');

  void SetNormalise(const int i) {
    normalise = i;
  };
  MatrixXd &GetTransformed() { return TransformedCentered; }

  MatrixXd &GetX() { return X; }

  void RunLPCA();

  void RunLPCAWithReductionOfComponents();

  void RunRevertLPCAWithReductionOfComponents();

  void RunRevertLPCA();

  void Print();

  void WriteTransformed(std::string);

  void WriteEigenvectors(std::string);

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {}

private:
  MatrixXd X, Xcentered, C, K, eigenvectors, Transformed, TransformedCentered;
  VectorXd eigenvalues, cumulative;
  unsigned int normalise;

public:
  template <typename F> void UploadPopulation(Population<F> &pop) {
    for (int i = 0; i < pop.size(); ++i) {
      auto individual = pop.GetTGenes(i);
      for (int j = 0; j < individual.size(); ++j) {
        auto gene = individual[j];
        if (X.rows() < i + 1) {
          X.conservativeResize(i + 1, X.cols());
        }
        if (X.cols() < j + 1) {
          X.conservativeResize(X.rows(), j + 1);
        }
        X(i, j) = gene.GetGAValue();
      }
    }
    std::string sep = "\n----------------------------------------\n";
    std::cout << X << sep;
    Xcentered.resize(X.rows(), X.cols());
  }
};
}

#endif
