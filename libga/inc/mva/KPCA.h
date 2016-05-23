#ifndef __KPCA__
#define __KPCA__

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "PCA.h"

namespace geantvmoop {

using namespace Eigen;

class KPCA : public PCA<KPCA> {

public:
  KPCA()
      : components(2), kernel_type(1), normalise(0), gamma(0.001),
        constant(1.0), order(2.0) {}
  explicit KPCA(MatrixXd &d)
      : components(2), kernel_type(1), normalise(0), gamma(0.001),
        constant(1.0), order(2.0) {
    X = d;
  }
  virtual ~KPCA() {}
  void LoadData(const char *data, char sep = ',');
  template <typename F> void UploadPopulation(Population<F> &pop);
  template <typename F> void LoadUpdatedPopulation(Population<F> &pop);
  void SetComponents(const int i) { components = i; };
  void SetKernel(const int i) { kernel_type = i; };
  void SetNormalise(const int i) { normalise = i; };
  void SetGamma(const double i) { gamma = i; };
  void SetConstant(const double i) { constant = i; };
  void SetOrder(const double i) { order = i; };
  MatrixXd &GetTransformed() { return transformed; }
  void RunKPCA();
  void Print();
  void WriteTransformed(std::string);
  void WriteEigenvectors(std::string);

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {}

private:
  double Kernel(const VectorXd &a, const VectorXd &b);
  MatrixXd X, Xcentered, C, K, eigenvectors, transformed;
  VectorXd eigenvalues, cumulative;
  unsigned int components, kernel_type, normalise;
  double gamma, constant, order;
};
}

#endif
