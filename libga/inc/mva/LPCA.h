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

public:
  template <typename F> void UploadPopulation(Population<F> &pop) {
    std::vector<double> fParameters;
    for (int i = 0; i < pop.size(); ++i) {
      auto individual = pop.GetTGenes(i);
      for (int i = 0; i < individual.size(); ++i) {
        auto gene = individual[i];
        fParameters.push_back(gene.GetGAValue());
        X.row(i) = VectorXd::Map(&fParameters.at(i), sizeof(double));
      }
    }
    std::string sep = "\n----------------------------------------\n";
    std::cout << X << sep;
  }
};
}

#endif