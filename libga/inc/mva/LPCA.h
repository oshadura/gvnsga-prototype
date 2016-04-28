#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "PCA.h"

using namespace Eigen;

class LPCA : public PCA<LPCA> {
public:
  LPCA() : normalise(0) {}
  explicit LPCA(MatrixXd &d) : normalise(0) { X = d; }
  void LoadData(const char *data, char sep = ',');
  void UploadPopulation(Population<double> &pop);
  void SetNormalise(const int i) {
    normalise = i;
  };
  MatrixXd &GetTransformed() { return transformed; }
  MatrixXd &GetX(){return X;}
  void RunPca();
  void Print();
  void WriteTransformed(std::string);
  void WriteEigenvectors(std::string);

private:
  MatrixXd X, Xcentered, C, K, eigenvectors, transformed;
  VectorXd eigenvalues, cumulative;
  unsigned int normalise;
};
