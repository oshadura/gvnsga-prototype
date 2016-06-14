#pragma once

#ifndef __KPCA__
#define __KPCA__

#include <fstream>
#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "PCA.h"

namespace geantvmoop {

using namespace Eigen;

class KPCA : public PCA<KPCA> {

public:
  KPCA()
      : components(3), kernel_type(1), normalise(0), gamma(0.001),
        constant(1.0), order(2.0) {}
  explicit KPCA(MatrixXd &d)
      : components(3), kernel_type(1), normalise(0), gamma(0.001),
        constant(1.0), order(2.0) {
    X = d;
  }
  virtual ~KPCA() {}

  void SetComponents(const int i) {
    components = i;
  };

  void SetKernel(const int i) {
    kernel_type = i;
  };

  void SetNormalise(const int i) {
    normalise = i;
  };

  void SetGamma(const double i) {
    gamma = i;
  };

  void SetConstant(const double i) {
    constant = i;
  };

  void SetOrder(const double i) {
    order = i;
  };

  MatrixXd &GetTransformed() { return transformed; }

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {
    Population<F> result;
    UploadPopulation(pop);
    RunKPCAWithReductionOfComponents();
    UnloadPopulation(result, X);
    return result;
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

  void LoadData(const char *data, char sep = ',') {
    // Read data
    unsigned int row = 0;
    std::ifstream file(data);
    if (file.is_open()) {
      std::string line, token;
      while (std::getline(file, line)) {
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
      file.close();
      Xcentered.resize(X.rows(), X.cols());
    } else {
      std::cout << "Failed to read file " << data << std::endl;
    }
  }

  void RunKPCAWithReductionOfComponents() {
    // Fill kernel matrix
    K.resize(X.rows(), X.rows());
    for (unsigned int i = 0; i < X.rows(); i++) {
      for (unsigned int j = i; j < X.rows(); j++) {
        K(i, j) = K(j, i) = Kernel(X.row(i), X.row(j));
        // printf("k(%i,%i) = %f\n", i, j, K(i, j));
      }
    }
    std::cout << "Matrix X: \n" << K << std::endl;
    EigenSolver<MatrixXd> edecomp(K);
    eigenvalues = edecomp.eigenvalues().real();
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
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
    // http://stackoverflow.com/questions/5122804/sorting-with-lambda
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
    transformed.resize(X.rows(), components);
    for (unsigned int i = 0; i < X.rows(); i++) {
      for (unsigned int j = 0; j < components; j++) {
        for (int k = 0; k < K.rows(); k++) {
          transformed(i, j) += K(i, k) * eigenvectors(k, j);
        }
      }
    }
    std::cout << "Sorted eigenvalues:" << std::endl;
    for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
      if (eigenvalues(i) > 0) {
        std::cout << "PC " << i + 1 << ": Eigenvalue: " << eigenvalues(i);
        printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",
               eigenvalues(i) / eigenvalues.sum(),
               cumulative(i) / eigenvalues.sum());
      }
    }
    std::cout << std::endl;
    // std::cout << "Sorted eigenvectors:" << eigenvectors std::endl;
    // std::cout << "Transformed data:" << transformed << std::endl;
  }

  void RunKPCA() {
    // Fill kernel matrix
    K.resize(X.rows(), X.rows());
    for (unsigned int i = 0; i < X.rows(); i++) {
      for (unsigned int j = i; j < X.rows(); j++) {
        K(i, j) = K(j, i) = Kernel(X.row(i), X.row(j));
        //printf("k(%i,%i) = %f\n", i, j, K(i, j));
      }
    }
    std::cout << "Matrix X: \n" << K << std::endl;
    EigenSolver<MatrixXd> edecomp(K);
    eigenvalues = edecomp.eigenvalues().real();
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
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
    transformed.resize(X.rows(), components);
    for (unsigned int i = 0; i < X.rows(); i++) {
      for (unsigned int j = 0; j < components; j++) {
        for (int k = 0; k < K.rows(); k++) {
          transformed(i, j) += K(i, k) * eigenvectors(k, j);
        }
      }
    }
    std::cout << "Sorted eigenvalues:" << std::endl;
    for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
      if (eigenvalues(i) > 0) {
        std::cout << "PC " << i + 1 << ": Eigenvalue: " << eigenvalues(i);
        printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",
               eigenvalues(i) / eigenvalues.sum(),
               cumulative(i) / eigenvalues.sum());
      }
    }
  }

  void Print() {
    std::cout << "Input data:\n " << X << std::endl;
    std::cout << "Centered data: \n" << Xcentered << std::endl;
    std::cout << "Covariance matrix: \n" << C << std::endl;
    std::cout << "Eigenvalues:\n " << eigenvalues << std::endl;
    std::cout << "Eigenvectors:\n " << eigenvectors << std::endl;
    std::cout << "Sorted eigenvalues: \n" << std::endl;
    for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
      if (eigenvalues(i) > 0) {
        std::cout << "PC " << i + 1 << ": Eigenvalue: \n" << eigenvalues(i);
        printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",
               eigenvalues(i) / eigenvalues.sum(),
               cumulative(i) / eigenvalues.sum());
      }
    }
    std::cout << "Sorted eigenvectors:\n" << eigenvectors << std::endl;
    std::cout << "Transformed data:\n" << X *eigenvectors << std::endl;
    //std::cout << "Transformed centred data:\n" << transformed << std::endl;
  }

  void WriteTransformed(std::string file) {
    std::ofstream outfile(file);
    for (unsigned int i = 0; i < transformed.rows(); i++) {
      for (unsigned int j = 0; j < transformed.cols(); j++) {
        outfile << transformed(i, j);
        if (j != transformed.cols() - 1)
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

private:
  double Kernel(const VectorXd &a, const VectorXd &b) {
    switch (kernel_type) {
    case 2:
      return (std::pow(a.dot(b) + constant, order));
    default:
      return (std::exp(-gamma * ((a - b).squaredNorm())));
    }
  }

private:
  MatrixXd X, Xcentered, C, K, eigenvectors, transformed;
  VectorXd eigenvalues, cumulative;
  unsigned int components, kernel_type, normalise;
  double gamma, constant, order;
};
}

#endif
