//===--- SpectralCLustering.h - LibGA
//---------------------------------------------*-
// C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file SpectralClustering.h
 * @brief
 */
//===----------------------------------------------------------------------===//

#pragma once

#ifndef __SPECTRALCLUSTERING__
#define __SPECTRALCLUSTERING__

#include "generic/Population.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <vector>

using namespace Eigen;

namespace geantvmoop {

class Spectral {

public:
  Spectral()
      : centers(2), kernel_type(2), normalise(1), max_iters(1000), gamma(0.001),
        constant(1.0), order(2.0) {}
  explicit Spectral(MatrixXd &d)
      : centers(2), kernel_type(2), normalise(1), max_iters(1000), gamma(0.001),
        constant(1.0), order(2.0) {
    X = d;
  }
  void set_centers(const unsigned int i) { centers = i; };
  void set_kernel(const unsigned int i) { kernel_type = i; };
  void set_normalise(const unsigned int i) { normalise = i; };
  void set_gamma(const double i) { gamma = i; };
  void set_constant(const double i) { constant = i; };
  void set_order(const double i) { order = i; };
  void set_max_iters(const unsigned int i) { max_iters = i; };
  const std::vector<int> &get_assignments() const { return assignments; };

private:
  MatrixXd X, K, eigenvectors;
  VectorXd eigenvalues, cumulative;
  unsigned int centers, kernel_type, normalise, max_iters;
  double gamma, constant, order;
  std::vector<int> assignments;

public:
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

  template <typename F,  std::size_t SizePop> void UploadPopulation(Population<F, SizePop> &pop) {
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
        X(i, j) = gene.GetGAValue();
      }
    }
    std::string sep = "\n----------------------------------------\n";
    std::cout << X << sep;
  }

  template <typename F,  std::size_t SizePop>
  void UnloadPopulation(Population<F, SizePop> &newpop, MatrixXd &data) {
    // check if they are both the same size!
    if (data.cols() != newpop.size())
      return;
    typename F::Input ind;
    std::vector<individual_t<F>> poplist;
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
      poplist.push_back(std::make_shared<geantvmoop::TGenes<F>>(newind));
      ind.clear();
    }
    newpop = Population<F, SizePop>(poplist);
  }

  int WriteData(const char *data, char sep) {
    unsigned int row = 0;
    std::ofstream file(data);
    if (file.is_open()) {
      for (unsigned int i = 0; i < X.rows(); i++) {
        for (unsigned int j = 0; j < X.cols(); j++) {
          file << X(i, j) << sep;
        }
        file << assignments[row] << std::endl;
        row++;
      }
      file.close();
      return (1);
    } else {
      return (0);
    }
  }

  double Kernel(const VectorXd &a, const VectorXd &b) {
    switch (kernel_type) {
    case 2:
      return (std::pow(a.dot(b) + constant, order));
    default:
      return (std::exp(-gamma * ((a - b).squaredNorm())));
    }
  }

  void GenerateKernelMatrix() {
    K.resize(X.rows(), X.rows());
    for (unsigned int i = 0; i < X.rows(); i++) {
      for (unsigned int j = i; j < X.rows(); j++) {
        K(i, j) = K(j, i) = Kernel(X.row(i), X.row(j));
        // if(i == 0) cout << K(i,j) << " ";
      }
    }
    VectorXd d = K.rowwise().sum();
    for (unsigned int i = 0; i < d.rows(); i++) {
      d(i) = 1.0 / sqrt(d(i));
    }
    auto F = d.asDiagonal();
    MatrixXd l = (K * F);
    for (unsigned int i = 0; i < l.rows(); i++) {
      for (unsigned int j = 0; j < l.cols(); j++) {
        l(i, j) = l(i, j) * d(i);
      }
    }
    K = l;
  }

  void EigenDecomposition() {
    EigenSolver<MatrixXd> edecomp(K);
    eigenvalues = edecomp.eigenvalues().real();
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
    std::vector<std::pair<double, VectorXd>> eigen_pairs;
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
                 const std::pair<double, VectorXd> b) -> bool {
                return (a.first > b.first);
              });
    if (centers > eigen_pairs.size())
      centers = eigen_pairs.size();
    for (unsigned int i = 0; i < eigen_pairs.size(); i++) {
      eigenvalues(i) = eigen_pairs[i].first;
      c += eigenvalues(i);
      cumulative(i) = c;
      eigenvectors.col(i) = eigen_pairs[i].second;
    }
    std::cout << "Sorted eigenvalues:" << std::endl;
    for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
      if (eigenvalues(i) > 0) {
        std::cout << "PC " << i + 1 << ": Eigenvalue: " << eigenvalues(i);
        printf("\t(%3.3f of variance, cumulative = % 3.3f)\n",
               eigenvalues(i) / eigenvalues.sum(),
               cumulative(i) / eigenvalues.sum());
        std::cout << eigenvectors.col(i) << std::endl;
      }
    }
    std::cout << std::endl;
    MatrixXd tmp = eigenvectors;
    // Select top K eigenvectors where K = centers
    eigenvectors = tmp.block(0, 0, tmp.rows(), centers);
  }

  void Cluster() {
    std::cout << "Generating Kernel matrix...." << std::endl;
    GenerateKernelMatrix();
    std::cout << "Managing eigen decomposition.." << std::endl;
    EigenDecomposition();
    std::cout << "Calculating k-means... " << std::endl;
    Kmeans();
  }

  void Kmeans() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rand_index(0, eigenvectors.rows() - 1);
    MatrixXd centroids = MatrixXd::Zero(centers, eigenvectors.cols());
    MatrixXd old_centroids;
    std::vector<int> rands;
    while (rands.size() < centers) {
      int r = rand_index(gen);
      bool tag = false;
      for (unsigned int j = 0; j < rands.size(); ++j) {
        if (r == rands[j]) {
          tag = true;
          break;
        }
      }
      if (!tag) {
        centroids.row(rands.size()) = eigenvectors.row(r);
        rands.push_back(r);
      }
    }
    MatrixXd id = MatrixXd::Identity(centers, centers);
    VectorXd minvals(eigenvectors.rows());
    // Matrix to map vectors to centroids
    MatrixXd post(eigenvectors.rows(), centers);
    int r, c;
    double old_e = 0;
    for (unsigned int n = 0; n < max_iters; n++) {
      old_centroids = centroids;
      // Calculate distances
      MatrixXd d2(eigenvectors.rows(), centers);
      for (unsigned int j = 0; j < centers; j++) {
        for (int k = 0; k < eigenvectors.rows(); k++) {
          d2(k, j) = (eigenvectors.row(k) - centroids.row(j)).squaredNorm();
        }
      }
      // Assign to nearest centroid
      for (unsigned int k = 0; k < eigenvectors.rows(); k++) {
        // Get index of centroid
        minvals[k] = d2.row(k).minCoeff(&r, &c);
        // Set centroid
        post.row(k) = id.row(c);
      }
      // Adjust centeroids
      VectorXd num_points = post.colwise().sum();
      for (unsigned int j = 0; j < centers; j++) {
        if (num_points(j) > 0) {
          MatrixXd s = MatrixXd::Zero(1, eigenvectors.cols());
          for (unsigned int k = 0; k < eigenvectors.rows(); k++) {
            if (post(k, j) == 1) {
              s += eigenvectors.row(k);
            }
          }
          centroids.row(j) = s / num_points[j];
        }
      }
      // Calculate error - total squared distance from centroids
      double e = minvals.sum();
      double ediff = std::fabs(old_e - e);
      double cdiff = (centroids - old_centroids).cwiseAbs().maxCoeff();
      printf("Iterations %i : Error %2.4f : Error delta %2.4f : Centroid "
             "movement %2.4f\n",
             n + 1, e, ediff, cdiff);
      if (n && cdiff < std::numeric_limits<double>::epsilon() &&
          ediff < std::numeric_limits<double>::epsilon()) {
        break;
      }
      old_e = e;
    }
    std::map<int, int> data_to_cluster;
    for (unsigned int j = 0; j < centers; j++) {
      for (int k = 0; k < eigenvectors.rows(); k++) {
        if (post(k, j) == 1) {
          data_to_cluster[k] = j + 1;
        }
      }
    }
    for (int k = 0; k < eigenvectors.rows(); k++) {
      assignments.push_back(data_to_cluster[k]);
    }
  }
};
}
#endif
