#pragma once

#ifndef __ROBUSTPCA__
#define __ROBUSTPCA__

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <unsupported/Eigen/MatrixFunctions>

#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"

#include "PCA.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iostream>
#include <cmath>

namespace geantvmoop {

using namespace Eigen;

class RobustPCA : public PCA<RobustPCA> {

private:
  MatrixXd D, A, E;

public:
  RobustPCA() {}

  virtual ~RobustPCA() {}

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
          if (D.rows() < row + 1) {
            D.conservativeResize(row + 1, D.cols());
          }
          if (D.cols() < col + 1) {
            D.conservativeResize(D.rows(), col + 1);
          }
          D(row, col) = std::atof(token.c_str());
          col++;
        }
        row++;
      }
      reader.close();
      //std::cout << "Filed was proccessed.." << std::endl;
      MatrixXd A = MatrixXd::Zero(D.rows(), D.cols());
      MatrixXd E = MatrixXd::Zero(D.rows(), D.cols());
      //D = D.array() * D.array();
    } else {
      std::cout << "Failed to read file..." << data << std::endl;
    }
  }

  template <typename F> void UploadPopulation(Population<F> &pop) {
    for (int i = 0; i < pop.size(); ++i) {
      auto individual = pop.GetTGenes(i);
      for (int j = 0; j < individual.size(); ++j) {
        auto gene = individual[j];
        if (D.rows() < i + 1) {
          D.conservativeResize(i + 1, D.cols());
        }
        if (D.cols() < j + 1) {
          D.conservativeResize(D.rows(), j + 1);
        }
        D(i, j) = gene.GetGAValue();
      }
    }
    MatrixXd A = MatrixXd::Zero(D.rows(), D.cols());
    MatrixXd E = MatrixXd::Zero(D.rows(), D.cols());
    //D = D.array() * D.array();
    std::string sep = "\n----------------------------------------\n";
    std::cout << D << sep;
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

  template <typename F> Population<F> MVAImpl(Population<F> &pop) {
    Population<F> result;
    UploadPopulation(pop);
    RobustPCAInexact();
    UnloadPopulation(result, A);
    return result;
  }
  // Double suppose to be F::Input
  int LargerThan(const VectorXd &v, double value) {
    int count = 0;
    for (int i = 0; i < v.size(); ++i) {
      if (v[i] > value)
        ++count;
    }
    return count;
  }

  void Print() {
    std::cout << "The original matirix; D = \n" << D << std::endl;
    std::cout << "Estimated row rank matrix: A = \n" << A << std::endl;
    std::cout << "Estimated sparse matrix: E = \n" << E << std::endl;
    std::cout << "Reconstructed matrix: A + E = \n" << A + E << std::endl;
    std::cout << "Reconstruction Error = \n" << (D - (A + E)).norm() << std::endl;
  }

  /**
   * Robust Principal Component Analysis (RPCA) using the inexact augumented
   * Lagrance multiplier.
   * @param D observation matrix (D = A + E)
   * @param A  row-rank matrix
   * @param E  sparse matrix
   */
  void RobustPCAInexact() {

    const int M = D.rows();
    const int N = D.cols();

    // supplementary variable
    A = MatrixXd::Zero(M, N);
    E = MatrixXd::Zero(M, N);
    MatrixXd Y = D;
    ArrayXXd zero = ArrayXXd::Zero(M, N);

    // Parameters
    const double lambda = 1.0 / sqrt(std::max(M, N));
    const double rho = 1.5;

    JacobiSVD<MatrixXd> svd_only_singlar_values(Y);

    const double norm_two =
        svd_only_singlar_values.singularValues()(0); // can be tuned
    const double norm_inf = Y.array().abs().maxCoeff() / lambda;
    const double dual_norm = std::max(norm_two, norm_inf);
    const double d_norm = D.norm();

    Y /= dual_norm;

    double mu = 1.25 / norm_two;

    const double mu_bar = mu * 1.0e+7;

    bool converged = false;

    int max_iter = 1000;

    double error_tolerance = 1.0e-7;

    int iter = 0;

    int total_svd = 0;

    int sv = 10;
    
    while (!converged) {
      // update sparse matrix E
      ArrayXXd temp_T = D - A + (1.0 / mu) * Y;
      E = (temp_T - lambda / mu).max(zero) + (temp_T + lambda / mu).min(zero);

      // force non-negative
      E = E.array().max(zero);

      // SVD
      JacobiSVD<MatrixXd> svd(D - E + 1.0 / mu * Y,
                              ComputeFullU | ComputeFullV);
      MatrixXd U = svd.matrixU();
      MatrixXd V = svd.matrixV();
      VectorXd singularValues = svd.singularValues();

      // truncate dimention
      int svp = LargerThan(singularValues, 1 / mu);
      if (svp < sv) {
        sv = std::min(svp + 1, N);
      } else {
        sv = std::min(svp + static_cast<int>(0.05 * N + 0.5), N);
      }

      // update A
      MatrixXd S_th =
          (singularValues.head(svp).array() - 1.0 / mu).matrix().asDiagonal();
      A = U.leftCols(svp) * S_th * V.leftCols(svp).transpose();

      // force non-negative
      A = A.array().max(zero);

      total_svd += 1;
      MatrixXd Z = D - A - E;
      Y = Y + mu * Z;
      mu = std::min(mu * rho, mu_bar);

      // objective function
      double objective = Z.norm() / d_norm;

      if (objective < error_tolerance) {
        converged = true;
      }

      if (++iter >= max_iter) {
        break;
      }
    }
  }
};
}

#endif