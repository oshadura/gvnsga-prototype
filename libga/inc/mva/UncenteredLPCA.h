#pragma once

#ifndef __UncenteredLPCA__
#define __UncenteredLPCA__

#include <fstream>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <unsupported/Eigen/MatrixFunctions>

#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"

#include "PCA.h"

// Uncentered version!

namespace geantvmoop {

using namespace Eigen;

class UncenteredLPCA : public PCA<UncenteredLPCA> {

private:
  MatrixXd X, C, K, eigenvectors, Transformed,
      covariance, dev, mean, devnew, meannew;
  VectorXd eigenvalues, cumulative, stddev, colmean, stddevnew, colmeannew;
  unsigned int normalise;

public:
  //We need to normalize for
  // \hat{U}^{\,t }\cdot  \hat{U}=  \hat I, \quad {U}^{\,t }_{i,i'}   {U}_{i',j}=\delta_{i,j}.
  UncenteredLPCA() : normalise(1) {}

  explicit UncenteredLPCA(MatrixXd &d) : normalise(1) { X = d; }

  virtual ~UncenteredLPCA() {}

  void SetNormalise(const int i) {
    normalise = i;
  };
  MatrixXd &GetTransformed() { return Transformed; }

  MatrixXd &GetX() { return X; }

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {
    Population<F> result;
    UploadPopulation(pop);
    RunUncenteredLPCAWithReductionOfComponents();
    UnloadPopulation(result, X);
    return result;
  }

  void LoadData(const char *data, char sep = ',') {
    unsigned int row = 0;
    std::ifstream reader;
    reader.open(data);
    if (reader.is_open()) {
      std::string line, token;
      while (std::getline(reader, line)) {
        std::stringstream tmp(line);
        unsigned int col = 0;
        while (std::getline(tmp, token, sep)) {
          if (X.rows() < row + 1) {
            X.conservativeResize(row + 1, X.cols());
          }
          if (X.cols() < col + 1) {
            X.conservativeResize(X.rows(), col + 1);
          }
          X(row, col) = std::atof(token.c_str());
          col++;
        }
        row++;
      }
      reader.close();
    } else {
      std::cout << "Failed to read file..." << data << std::endl;
    }
  }

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
  }

  template <typename F>
  void UnloadPopulation(Population<F> &newpop, MatrixXd &data) {
    typename F::Input ind;
    std::vector<individual_t<F> > poplist;
    std::string sep = "\n----------------------------------------\n";
    for (int i = 0; i < data.rows(); ++i) {
      for (int j = 0; j < data.cols(); ++j) {
        std::cout << "Gene to be added in a population[" << i << "," << j << "] is " << data(i, j) << std::endl;
        ind.push_back(data(i, j));
      }
      std::cout << "New gene added." << std::endl;
      TGenes<F> newind = ind;
      poplist.push_back(std::make_shared<geantvmoop::TGenes<F> >(newind));
      ind.clear();
    }
    newpop = Population<F>(poplist);
  }

  void RunUncenteredLPCA() {
    C = (X.adjoint() * X) / double(X.rows());
    EigenSolver<MatrixXd> edecomp(C);
    // Eigen values
    eigenvalues = edecomp.eigenvalues().real();
    // Eigen vectors
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
    // Eigen pairs [eigenvalue, eigenvector]
    std::vector<std::pair<double, VectorXd> > eigen_pairs;
    double c = 0.0;

    for (unsigned int i = 0; i < eigenvectors.cols(); i++) {
      if (normalise) {
        double norm = eigenvectors.col(i).norm();
        eigenvectors.col(i) /= norm;
      }
      eigen_pairs.push_back(
          std::make_pair(eigenvalues(i), eigenvectors.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(eigen_pairs.begin(), eigen_pairs.end(),
              [](const std::pair<double, VectorXd> a,
                 const std::pair<double, VectorXd> b)
                  ->bool { return (a.first > b.first); });

    for (unsigned int i = 0; i < eigen_pairs.size(); i++) {
      eigenvalues(i) = eigen_pairs[i].first;
      c += eigenvalues(i);
      cumulative(i) = c;
      eigenvectors.col(i) = eigen_pairs[i].second;
    }
    Transformed = X * eigenvectors;
  }

  void RunUncenteredLPCAWithReductionOfComponents() {
    double totalvar = 0;
    int i = 0;
    C = (X.adjoint() * X) / double(X.rows());
    EigenSolver<MatrixXd> edecomp(C);
    // Eigen values
    eigenvalues = edecomp.eigenvalues().real();
    // Eigen vectors
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
    // Eigen pairs [eigenvalue, eigenvector]
    std::vector<std::pair<double, VectorXd> > eigen_pairs;
    double c = 0.0;
    for (unsigned int i = 0; i < eigenvectors.cols(); i++) {
      if (normalise) {
        double norm = eigenvectors.col(i).norm();
        eigenvectors.col(i) /= norm;
      }
      eigen_pairs.push_back(
          std::make_pair(eigenvalues(i), eigenvectors.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(eigen_pairs.begin(), eigen_pairs.end(),
              [](const std::pair<double, VectorXd> a,
                 const std::pair<double, VectorXd> b)
                  ->bool { return (a.first > b.first); });
    // Printing current state information before  PC cutoff
    std::cout << "Printing original information after PCA" << std::endl;
    Transformed = X * eigenvectors;
    // Varince based selection (< 85 %)
    while (totalvar <= 0.85) {
      eigenvalues(i) = eigen_pairs[i].first;
      c += eigenvalues(i);
      cumulative(i) = c;
      eigenvectors.col(i) = eigen_pairs[i].second;
      totalvar = totalvar + (eigenvalues(i) / eigenvalues.sum());
      ++i;
    }
    Print();
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "REVERSE PCA: " << std::endl;
    eigenvectors.conservativeResize(eigenvectors.rows(), i);
    Transformed.conservativeResize(Transformed.rows(), i);
    // TransformedCentered.conservativeResize(Transformed.rows(), i);
    std::cout << "Reduced eigenvectors:\n" << eigenvectors << std::endl;
    std::cout << "Reduced tranformed matrix \n" << Transformed << std::endl;
    std::cout << "Total number of components to be used in transformed matrix: "
              << i << std::endl;
    // Transformed matrix
    MatrixXd NewDataMatrix, NewDataMatrixTransposed;
    NewDataMatrix = eigenvectors * Transformed.transpose();
    NewDataMatrixTransposed = NewDataMatrix.transpose();
    std::cout << "Transformed data matrix:\n" << NewDataMatrixTransposed
              << std::endl;
    X = NewDataMatrixTransposed;
  }


  void Print() {
    std::cout << "Input data:\n" << X << std::endl;
    std::cout << "Mean of columns:\n" << colmean << std::endl;
    std::cout << "Covariance matrix:\n" << C << std::endl;
    std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;
    std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
    std::cout << "Sorted eigenvalues:\n " << std::endl;
    for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
      if (eigenvalues(i) > 0) {
        std::cout << "PC " << i + 1 << ": Eigenvalue:\n " << eigenvalues(i);
        printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",
               eigenvalues(i) / eigenvalues.sum(),
               cumulative(i) / eigenvalues.sum());
      }
    }
    std::cout << std::endl;
    std::cout << "Sorted eigenvectors:\n" << eigenvectors << std::endl;
    std::cout << "Transformed data:\n" << Transformed << std::endl;
  }

  void WriteTransformed(std::string file) {
    std::ofstream outfile(file);
    for (unsigned int i = 0; i < Transformed.rows(); i++) {
      for (unsigned int j = 0; j < Transformed.cols(); j++) {
        outfile << Transformed(i, j);
        if (j != Transformed.cols() - 1)
          outfile << ",";
      }
      outfile << std::endl;
    }
    outfile.close();
    std::cout << "Written file " << file << std::endl;
  }

  void WriteEigenvectors(std::string file) {
    std::ofstream outfile(file);
    for (unsigned int i = 0; i < eigenvectors.rows(); i++) {
      for (unsigned int j = 0; j < eigenvectors.cols(); j++) {
        outfile << eigenvectors(i, j);
        if (j != eigenvectors.cols() - 1)
          outfile << ",";
      }
      outfile << std::endl;
    }
    outfile.close();
    std::cout << "Written file " << file << std::endl;
  }
};
}

#endif
