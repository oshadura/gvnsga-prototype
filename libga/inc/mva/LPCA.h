#pragma once

#ifndef __LPCA__
#define __LPCA__

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

namespace geantvmoop {

using namespace Eigen;

class LPCA : public PCA<LPCA> {

private:
  MatrixXd X, Xcentered, C, K, eigenvectors, Transformed, TransformedCentered,
      covariance, dev, mean, devnew, meannew;
  VectorXd eigenvalues, cumulative, stddev, colmean, stddevnew, colmeannew;
  unsigned int normalise;

public:
  LPCA() : normalise(0) {}

  explicit LPCA(MatrixXd &d) : normalise(0) { X = d; }

  virtual ~LPCA() {}

  void SetNormalise(const int i) {
    normalise = i;
  };
  MatrixXd &GetTransformed() { return TransformedCentered; }

  MatrixXd &GetX() { return X; }

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {
    Population<F> result;
    UploadPopulation(pop);
    RunLPCAWithReductionOfComponents();
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
      Xcentered.resize(X.rows(), X.cols());
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
    Xcentered.resize(X.rows(), X.cols());
  }

  template <typename F>
  void UnloadPopulation(Population<F> &newpop, MatrixXd &data) {
    // check if they are both the same size!
    // if (data.cols() != newpop.size())
    //  return;
    typename F::Input ind;
    std::vector<individual_t<F> > poplist;
    std::string sep = "\n----------------------------------------\n";
    for (int i = 0; i < data.rows(); ++i) {
      for (int j = 0; j < data.cols(); ++j) {
        // std::cout << "Gene to be added in a population[" << i << "," << j
        //          << "] is " << data(i, j) << std::endl;
        ind.push_back(data(i, j));
        // ind.SetGAValue(data(i, j));
      }
      // std::cout << "New gene added." << std::endl;
      TGenes<F> newind = ind;
      poplist.push_back(std::make_shared<geantvmoop::TGenes<F> >(newind));
      ind.clear();
    }
    newpop = Population<F>(poplist);
  }

  void RunLPCA() {
    // Centered matrix
    Xcentered = X.rowwise() - X.colwise().mean();
    C = (Xcentered.adjoint() * Xcentered) / double(X.rows());
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
    // Transformed matrix
    TransformedCentered = Xcentered * eigenvectors;
  }

  void RunLPCAWithReductionOfComponents() {
    // std dev
    double totalvar = 0;
    int i = 0;
    // Covariance matrix
    covariance = (Xcentered.adjoint() * Xcentered) / double(X.rows() - 1);
    // Mean of matrix X
    colmean = X.colwise().sum() / X.rows();
    // Vector of std deviation of colums
    stddev =
        (X.rowwise() - colmean.transpose()).array().pow(2).colwise().sum() /
        X.rows();
    mean = MatrixXd::Zero(X.rows(), X.cols());
    dev = MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); ++i) {
      mean.row(i) = colmean.transpose();
      dev.row(i) = stddev.transpose();
    }
    // Centered matrix
    Xcentered = (X.rowwise() - X.colwise().mean());    
    // Sqrt of sigma
    stddev = stddev.cwiseSqrt();
    Xcentered = Xcentered.array() / dev.array();
    C = (Xcentered.adjoint() * Xcentered) / double(X.rows());
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
    // TransformedCentered = Xcentered * eigenvectors;
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
    std::cout << "Sqrt of std dev vector of matrix X:\n" << stddev << std::endl;
    std::cout << "Sqrt of std dev matrix of matrix X:\n" << dev << std::endl;
    std::cout << "Column nean vector of matrix X:\n" << colmean << std::endl;
    std::cout << "Mean matrix:\n" << mean << std::endl;
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
    colmeannew = NewDataMatrixTransposed.colwise().sum() /
                 NewDataMatrixTransposed.rows();
    // Vector of std deviation of colums
    stddevnew = (NewDataMatrixTransposed.rowwise() - colmeannew.transpose())
                    .array()
                    .pow(2)
                    .colwise()
                    .sum() /
                NewDataMatrixTransposed.rows();
    std::cout << "New mean vector of matrix X':\n" << colmeannew << std::endl;
    std::cout << "New std dev of matrix X':\n" << stddevnew << std::endl;
    meannew = MatrixXd::Zero(X.rows(), X.cols());
    devnew = MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); ++i) {
      devnew.row(i) = stddevnew.transpose();
      meannew.row(i) = colmeannew.transpose();
    }
    NewDataMatrixTransposed = NewDataMatrixTransposed - meannew;
    NewDataMatrixTransposed = NewDataMatrixTransposed.array()/devnew.array();
    std::cout << "New Transformed data matrix X':\n" << NewDataMatrixTransposed << std::endl;
    /////////////////BACK///////////////////
    X =  dev.array() * NewDataMatrixTransposed.array();
    X = X + mean;
     std::cout << "New Transformed data matrix with reverse = X:\n" << X << std::endl;
  }


  void Print() {
    std::cout << "Input data:\n" << X << std::endl;
    std::cout << "Mean of columns:\n" << colmean << std::endl;
    std::cout << "Centered data:\n" << Xcentered << std::endl;
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
    for (unsigned int i = 0; i < TransformedCentered.rows(); i++) {
      for (unsigned int j = 0; j < TransformedCentered.cols(); j++) {
        outfile << TransformedCentered(i, j);
        if (j != TransformedCentered.cols() - 1)
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
