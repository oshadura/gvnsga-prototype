#ifndef __LPCA__
#define __LPCA__

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <generic/TGenes.h>
#include <generic/Functions.h>

#include "PCA.h"

namespace geantvmoop {

using namespace Eigen;

class LPCA : public PCA<LPCA> {
public:
  LPCA() : normalise(0) {}
  explicit LPCA(MatrixXd &d) : normalise(0) { X = d; }
  virtual ~LPCA() {}
  void LoadData(const char *data, char sep = ',');
  template <typename F> void UploadPopulation(Population<F> &pop);
  template <typename F> void LoadUpdatedPopulation(Population<F> &pop);
  void SetNormalise(const int i) {
    normalise = i;
  };
  MatrixXd &GetTransformed() { return transformed; }
  MatrixXd &GetX() { return X; }
  void RunLPCA();
  void Print();
  void WriteTransformed(std::string);
  void WriteEigenvectors(std::string);

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {}

private:
  MatrixXd X, Xcentered, C, K, eigenvectors, transformed;
  VectorXd eigenvalues, cumulative;
  unsigned int normalise;
};
}

#endif