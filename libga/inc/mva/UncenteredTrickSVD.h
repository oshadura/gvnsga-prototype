#pragma once

#ifndef __UncenteredTrickSVD__
#define __UncenteredTrickSVD__

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

class UncenteredTrickSVD : public PCA<UncenteredTrickSVD> {

private:
  MatrixXd X, Xtrick, C, K, eigenvectors, Transformed, covariance, dev, mean,
      devnew, meannew;
  VectorXd eigenvalues, cumulative, stddev, colmean, stddevnew, colmeannew;
  unsigned int normalise;

public:
  // We need to normalize for
  // \hat{U}^{\,t }\cdot  \hat{U}=  \hat I, \quad {U}^{\,t }_{i,i'}
  // {U}_{i',j}=\delta_{i,j}.
  UncenteredTrickSVD() : normalise(1) {}

  explicit UncenteredTrickSVD(MatrixXd &d) : normalise(1) { X = d; }

  virtual ~UncenteredTrickSVD() {}

  void SetNormalise(const int i) { normalise = i; };
  
  MatrixXd &GetTransformed() { return Transformed; }

  MatrixXd &GetX() { return X; }

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {
    Population<F> result;
    UploadPopulation(pop);
    RunUncenteredTrickSVDWithReductionOfComponents();
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
      Xtrick = X;
      X = X.rightCols(X.cols()); // was 2
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
        if (std::isinf(gene.GetGAValue())) {
          gene.SetGAValue(1);
        }
        if (std::isnan(gene.GetGAValue())) {
          gene.SetGAValue(0);
        }
        X(i, j) = gene.GetGAValue();
      }
    }
    Xtrick = X;
    X = X.rightCols(X.cols()); // was 2
    std::string sep = "\n----------------------------------------\n";
    MatrixXd Y = Xtrick;
    JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV);
    std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
    std::cout << "Its left singular vectors are the columns of the thin U matrix:"
         << std::endl << svd.matrixU() << std::endl;
    std::cout << "Its right singular vectors are the columns of the thin V matrix:"
         << std::endl << svd.matrixV() << std::endl;
  }

  template <typename F>
  void UnloadPopulation(Population<F> &newpop, MatrixXd &data) {
    typename F::Input ind;
    MatrixXd population;
    population.conservativeResize(Xtrick.rows(), Xtrick.cols());
    std::cout << data << "\n";
    //std::cout << Xtrick.leftCols(2) << "\n"; // was 2
    //population << Xtrick.leftCols(2), data; // was 2
    population << data;
    std::cout << "Finally..\n" << population << std::endl;
    std::vector<individual_t<F>> poplist;
    std::string sep = "\n----------------------------------------\n";
    for (int i = 0; i < population.rows(); ++i) {
      for (int j = 0; j < population.cols(); ++j) {
        //std::cout << "Gene to be added in a population[" << i << "," << j
        //          << "] is " << population(i, j) << std::endl;
        ind.push_back(population(i, j));
      }
      //std::cout << "New gene added." << std::endl;
      TGenes<F> newind = ind;
      poplist.push_back(std::make_shared<geantvmoop::TGenes<F>>(newind));
      ind.clear();
    }
    newpop = Population<F>(poplist);
  }

  void RunUncenteredTrickSVD() {
    JacobiSVD<MatrixXd> svd(X, ComputeFullU | ComputeFullV);
    auto s = svd.singularValues();
    std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout << "Its left singular vectors are the columns of the thin U matrix:"
         << std::endl << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout << "Its right singular vectors are the columns of the thin V matrix:"
         << std::endl << svd.matrixV() << std::endl;
    // C = (X.adjoint() * X) / double(X.rows());
    // EigenSolver<MatrixXd> edecomp(C);
    // Eigen values
    // eigenvalues = edecomp.eigenvalues().real();
    // Eigen vectors
    // eigenvectors = edecomp.eigenvectors().real();
    // cumulative.resize(eigenvalues.rows());
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
    std::cout << "Check::Transformed back data matrix:\n"
              << NewDataMatrixTransposed << std::endl;  
}

  void RunUncenteredTrickSVDWithReductionOfComponents() {
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    auto s = svd.singularValues();
    std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout << "Its left singular vectors are the columns of the thin U matrix:"
         << std::endl << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout << "Its right singular vectors are the columns of the thin V matrix:"
         << std::endl << svd.matrixV() << std::endl;
    double totalvar = 0;
    int i = 0;
    //C = (X.adjoint() * X) / double(X.rows());
    //EigenSolver<MatrixXd> edecomp(C);
    // Eigen values
    //eigenvalues = edecomp.eigenvalues().real();
    // Eigen vectors
    //eigenvectors = edecomp.eigenvectors().real();
    //cumulative.resize(eigenvalues.rows());
    // Eigen pairs [eigenvalue, eigenvector]
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;
    for (unsigned int i = 0; i < rvec.cols(); i++) {
      if (normalise) {
        double norm = rvec.col(i).norm();
        rvec.col(i) /= norm;
      }
      fEigenValues.push_back(
          std::make_pair(s(i), rvec.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(
        fEigenValues.begin(), fEigenValues.end(),
        [](std::pair<double, VectorXd> &a, std::pair<double, VectorXd> &b) {
          if (a.first > b.first)
            return true;
          if (a.first == b.first)
            return a.first > b.first;
          return false;
        });
    for (unsigned int i = 0; i < rvec.cols(); i++) {
      s(i) = fEigenValues[i].first;
      rvec.col(i) = fEigenValues[i].second;
    }
    // Printing current state information before  PC cutoff
    std::cout << "Printing original information after PCA" << std::endl;
    //Transformed = X * eigenvectors;
    // Varince based selection (< 95 %)
    while (totalvar <= 0.95) {
      s(i) = fEigenValues[i].first;
      c += s(i);
      //cumulative(i) = c;
      rvec.col(i) = fEigenValues[i].second;
      totalvar = totalvar + (s(i) / s.sum());
      ++i;
    }
    //Print();
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "REVERSE SVD: " << std::endl;
    for(unsigned int j = 0; j < X.rows(); i++){
     for(unsigned int i = 0; j < X.cols(); j++){ 
       MatrixXd eig = s.asDiagonal();
       Xtrick += std::sqrt(X.rows())*(std::sqrt(s(i))*rvec.col(i)*eig.row(i));
       
     }
    }
    //eigenvectors.conservativeResize(eigenvectors.rows(), i);
    //Transformed.conservativeResize(Transformed.rows(), i);
    // TransformedCentered.conservativeResize(Transformed.rows(), i);
    //std::cout << "Reduced eigenvectors:\n" << eigenvectors << std::endl;
    //std::cout << "Reduced tranformed matrix \n" << Transformed << std::endl;
    //std::cout << "Total number of components to be used in transformed matrix: "
    //          << i << std::endl;
    ////Transformed matrix
    //MatrixXd NewDataMatrix, NewDataMatrixTransposed;
    //NewDataMatrix = eigenvectors * Transformed.transpose();
    //NewDataMatrixTransposed = NewDataMatrix.transpose();
    //std::cout << "Transformed data matrix:\n"
    //          << NewDataMatrixTransposed << std::endl;
    //X = NewDataMatrixTransposed .array().abs();
    //std::cout << "Done." << std::endl;
    X = Xtrick;
  }
  
/*
  void Print() {
    //std::cout << "Input data:\n" << Xtrick << std::endl;
#ifdef DEBUG
    //std::cout << "Mean of columns:\n" << colmean << std::endl;
#endif
    //std::cout << "Covariance matrix:\n" << C << std::endl;
    std::cout << "Eigenvalues:\n" << s << std::endl;
    std::cout << "Eigenvectors:\n" << rvec << std::endl;
    std::cout << "Sorted eigenvalues:\n " << std::endl;
    for (unsigned int i = 0; i < s.rows(); i++) {
      if (s(i) > 0) {
        std::cout << "PC " << i + 1 << ": Eigenvalue:\n " << s(i);
        printf("\t(%3.3f of variance)\n", s(i) / s.sum());
      }
    }
    std::cout << std::endl;
    std::cout << "Sorted eigenvectors:\n" << S << std::endl;
    //std::cout << "Transformed data:\n" << Transformed << std::endl;
  }
*/
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
