#pragma once

#ifndef __SVDREPRESENTATION__
#define __SVDREPRESENTATION__

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

class SVDRepresentation : public PCA<SVDRepresentation> {

private:
  MatrixXd D;

public:
  SVDRepresentation() {}

  virtual ~SVDRepresentation() {}

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
    SVDOutput();
    UnloadPopulation(result, D);
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
    std::cout << "SVD representation of population: \n" << D << std::endl;
  }

  /**
   * SVD representation of population
   * For Checking how Pareto Fronty is dependable from singular values..
   */
  void SVDOutput() {
    //Print();
    const int M = D.rows();
    const int N = D.cols();
    const double lambda = 1.0 / sqrt(std::max(M, N));
    MatrixXd Y = D;
    JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV);
    std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
    std::cout << "Its left singular vectors are the columns of the thin U matrix:"
         << std::endl << svd.matrixU() << std::endl;
    std::cout << "Its right singular vectors are the columns of the thin V matrix:"
         << std::endl << svd.matrixV() << std::endl;
    //Vector3f rhs(1, 0, 0);
    //cout << "Now consider this rhs vector:" << endl << rhs << endl;
    //cout << "A least-squares solution of m*x = rhs is:" << endl
    //     << svd.solve(rhs) << endl;
    const double norm_two = svd.singularValues()(0); // can be tuned
    const double norm_inf = Y.array().abs().maxCoeff() / lambda;
    const double dual_norm = std::max(norm_two, norm_inf);
    const double d_norm = D.norm();
  }
};
}

#endif