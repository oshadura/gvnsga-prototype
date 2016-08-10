#pragma once

#ifndef __UncenteredLPCA__
#define __UncenteredLPCA__

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>

#include "generic/GADouble.h"
#include "generic/GAVector.h"
#include "generic/Population.h"
#include "generic/TGenes.h"

#include "PCA.h"

// Uncentered version!

namespace geantvmoop {

using namespace Eigen;

class UncenteredLPCA : public PCA<UncenteredLPCA> {

private:
  MatrixXd X, C, K, eigenvectors, Transformed, covariance, dev, mean, devnew,
      meannew;
  VectorXd eigenvalues, cumulative, stddev, colmean, stddevnew, colmeannew,
      mean_column_centered;
  unsigned int normalise;

public:
  // We need to normalize for
  // \hat{U}^{\,t }\cdot  \hat{U}=  \hat I, \quad {U}^{\,t }_{i,i'}
  // {U}_{i',j}=\delta_{i,j}.
  UncenteredLPCA() : normalise(1) {}

  explicit UncenteredLPCA(MatrixXd &d) : normalise(1) { X = d; }

  virtual ~UncenteredLPCA() {}

  void SetNormalise(const int i) { normalise = i; };
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
    int row = 0;
    std::ifstream reader;
    reader.open(data);
    if (reader.is_open()) {
      std::string line, token;
      while (std::getline(reader, line)) {
        std::stringstream tmp(line);
        int col = 0;
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
    for (std::size_t i = 0; i < pop.size(); ++i) {
      auto individual = pop.GetTGenes(i);
      for (std::size_t j = 0; j < individual.size(); ++j) {
        auto gene = individual[j];
        if (X.rows() < (int)i + 1) {
          X.conservativeResize(i + 1, X.cols());
        }
        if (X.cols() < (int)j + 1) {
          X.conservativeResize(X.rows(), j + 1);
        }
        // Stupid thing, but it works..
        if (std::isinf(gene.GetGAValue())) {
          gene.SetGAValue(1);
        }
        if (std::isnan(gene.GetGAValue())) {
          gene.SetGAValue(0);
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
    std::vector<individual_t<F>> poplist;
    std::string sep = "\n----------------------------------------\n";
    for (int i = 0; i < data.rows(); ++i) {
      for (int j = 0; j < data.cols(); ++j) {
        std::cout << "Gene to be added in a population[" << i << "," << j
                  << "] is " << data(i, j) << std::endl;
        ind.push_back(data(i, j));
      }
      std::cout << "New gene added." << std::endl;
      TGenes<F> newind = ind;
      poplist.push_back(std::make_shared<geantvmoop::TGenes<F>>(newind));
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
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;

    for (unsigned int i = 0; i < eigenvectors.cols(); i++) {
      if (normalise) {
        double norm = eigenvectors.col(i).norm();
        eigenvectors.col(i) /= norm;
      }
      fEigenValues.push_back(
          std::make_pair(eigenvalues(i), eigenvectors.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(fEigenValues.begin(), fEigenValues.end(),
              [](const std::pair<double, VectorXd> &a,
                 const std::pair<double, VectorXd> &b) {
                if (a.first > b.first)
                  return true;
                if (a.first == b.first)
                  return a.first > b.first;
                return false;
              });
    for (unsigned int i = 0; i < fEigenValues.size(); i++) {
      eigenvalues(i) = fEigenValues[i].first;
      c += eigenvalues(i);
      cumulative(i) = c;
      eigenvectors.col(i) = fEigenValues[i].second;
    }
    Transformed = X * eigenvectors;
    // Checkout if we are right
    MatrixXd NewDataMatrix, NewDataMatrixTransposed;
    NewDataMatrix = eigenvectors * Transformed.transpose();
    NewDataMatrixTransposed = NewDataMatrix.transpose();
    std::cout << "Transformed back data matrix:\n"
              << NewDataMatrixTransposed << std::endl;
    X = NewDataMatrixTransposed.array();
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
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;
    for (unsigned int i = 0; i < eigenvectors.cols(); i++) {
      if (normalise) {
        double norm = eigenvectors.col(i).norm();
        eigenvectors.col(i) /= norm;
      }
      fEigenValues.push_back(
          std::make_pair(eigenvalues(i), eigenvectors.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(fEigenValues.begin(), fEigenValues.end(),
              [](const std::pair<double, VectorXd> &a,
                 const std::pair<double, VectorXd> &b) {
                if (a.first > b.first)
                  return true;
                if (a.first == b.first)
                  return a.first > b.first;
                return false;
              });
    // Printing current state information before  PC cutoff
    std::cout << "Printing original information after PCA" << std::endl;
    Transformed = X * eigenvectors;
    // Varince based selection (< 85 %)
    while (totalvar <= 0.95) {
      eigenvalues(i) = fEigenValues[i].first;
      c += eigenvalues(i);
      cumulative(i) = c;
      eigenvectors.col(i) = fEigenValues[i].second;
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
    std::cout << "Transformed data matrix:\n"
              << NewDataMatrixTransposed << std::endl;
    X = NewDataMatrixTransposed.array().abs();
  }

  void RunUncenteredLPCAWithReductionOfComponentsScale() {
    double totalvar = 0;
    int i = 0;
    colmean = X.colwise().sum() / X.rows();
    //====== Std deviation (Scaling) =================//
    stddev =
        (X.rowwise() - colmean.transpose()).array().pow(2).colwise().sum() /
        X.rows();
    mean = MatrixXd::Zero(X.rows(), X.cols());
    dev = MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); ++i) {
      mean.row(i) = colmean.transpose();
      dev.row(i) = stddev.transpose();
    }
    // Sqrt of sigma
    dev = dev.cwiseSqrt();
    X = X.array() / dev.array();
    //#ifdef DEBUG
    //=============== Output print===================//
    std::cout << "Std dev vector of matrix X:\n" << stddev << std::endl;
    std::cout << "Sqrt of std dev matrix of matrix X:\n" << dev << std::endl;
    std::cout << "Column nean vector of matrix X:\n" << colmean << std::endl;
    std::cout << "Mean matrix:\n" << mean << std::endl;
    std::cout << "---------------------------\n" << std::endl;
    C = (X.adjoint() * X) / double(X.rows());
    EigenSolver<MatrixXd> edecomp(C);
    // Eigen values
    eigenvalues = edecomp.eigenvalues().real();
    // Eigen vectors
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
    // Eigen pairs [eigenvalue, eigenvector]
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;
    for (unsigned int i = 0; i < eigenvectors.cols(); i++) {
      if (normalise) {
        double norm = eigenvectors.col(i).norm();
        eigenvectors.col(i) /= norm;
      }
      fEigenValues.push_back(
          std::make_pair(eigenvalues(i), eigenvectors.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(fEigenValues.begin(), fEigenValues.end(),
              [](const std::pair<double, VectorXd> &a,
                 const std::pair<double, VectorXd> &b) {
                if (a.first > b.first)
                  return true;
                if (a.first == b.first)
                  return a.first > b.first;
                return false;
              });
    // Printing current state information before  PC cutoff
    std::cout << "Printing original information after PCA" << std::endl;
    Transformed = X * eigenvectors;
    // Varince based selection (< 85 %)
    while (totalvar <= 0.95) {
      eigenvalues(i) = fEigenValues[i].first;
      c += eigenvalues(i);
      cumulative(i) = c;
      eigenvectors.col(i) = fEigenValues[i].second;
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
    mean_column_centered = NewDataMatrixTransposed.colwise().sum() /
                           NewDataMatrixTransposed.rows();
    stddevnew =
        (NewDataMatrixTransposed.rowwise() - mean_column_centered.transpose())
            .array()
            .pow(2)
            .colwise()
            .sum() /
        NewDataMatrixTransposed.rows();
    //#ifdef DEBUG
    std::cout << "New std dev of matrix X':\n" << stddevnew << std::endl;
    stddevnew = stddevnew.cwiseSqrt();
    //#endif
    meannew = MatrixXd::Zero(X.rows(), X.cols());
    devnew = MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); ++i) {
      devnew.row(i) = stddevnew.transpose();
      meannew.row(i) = mean_column_centered.transpose();
    }
    std::cout << "New sqrt of std dev of matrix X':\n" << devnew << std::endl;
    std::cout << "Transformed data matrix:\n"
              << NewDataMatrixTransposed << std::endl;
    X = NewDataMatrixTransposed.array().abs() * devnew.array();
  }

  void HistoFill() {}

  void Print() {
    std::cout << "Input data:\n" << X << std::endl;
    //#ifdef DEBUG
    std::cout << "Mean of columns:\n" << colmean << std::endl;
    //#endif
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
