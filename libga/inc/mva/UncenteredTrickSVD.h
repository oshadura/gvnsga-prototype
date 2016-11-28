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
    //RunUncenteredTrickSVDDiff();
    //RunUncenteredTrickSVDWithReductionOfComponents();
    //RunUncenteredTrickSVD_95Eigen();
    //RunUncenteredTrickSVDCutApproach();
    RunUncenteredTrickSVD2();
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
            Xtrick.conservativeResize(row + 1, Xtrick.cols());
          }
          if (X.cols() < col + 1) {
            X.conservativeResize(X.rows(), col + 1);
            Xtrick.conservativeResize(X.rows(), col + 1);
          }
          X(row, col) = std::atof(token.c_str());
          Xtrick(row, col) = std::atof(token.c_str());
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
    //Xtrick = X;
    X = X.rightCols(X.cols()); // was 2
    std::string sep = "\n----------------------------------------\n";
    //MatrixXd Y = Xtrick;
    /*
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
    std::cout << "Its left singular vectors are the columns of the thin U matrix:"
         << std::endl << svd.matrixU() << std::endl;
    std::cout << "Its right singular vectors are the columns of the thin V matrix:"
         << std::endl << svd.matrixV() << std::endl;
    */
  }

  template <typename F>
  void UnloadPopulation(Population<F> &newpop, MatrixXd &data) {
    typename F::Input ind;
    MatrixXd population;
    population.conservativeResize(X.rows(), X.cols());
    std::cout << data << "\n";
    //std::cout << Xtrick.leftCols(2) << "\n"; // was 2
    //population << Xtrick.leftCols(2), data; // was 2
    population << data;
    std::cout << "Finally..\n" << population << std::endl;
    std::vector<individual_t<F>> poplist;
    std::string sep = "\n----------------------------------------\n";
    for (int i = 0; i < population.rows(); ++i) {
      for (int j = 0; j < population.cols(); ++j) {
        // std::cout << "Gene to be added in a population[" << i << "," << j
        //          << "] is " << population(i, j) << std::endl;
        ind.push_back(population(i, j));
      }
      // std::cout << "New gene added." << std::endl;
      TGenes<F> newind = ind;
      poplist.push_back(std::make_shared<geantvmoop::TGenes<F>>(newind));
      ind.clear();
    }
    newpop = Population<F>(poplist);
  }

  void RunUncenteredSVD() {
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    auto s = svd.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svd.matrixV() << std::endl;
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
    int j = 0;
    MatrixXd Xsimple;
    Xsimple.conservativeResize(Xsimple.rows(), X.rows());
    Xsimple.conservativeResize(Xsimple.cols(), X.cols());
    while (j != X.rows()) {
      Xsimple +=
          /*std::sqrt(X.rows())*/ s(j) * lvec.col(j) * rvec.transpose().row(j);
      j++;
      std::cout << "Print matrix: \n" << Xsimple << std::endl;
    }
    // Transformed = X * eigenvectors;
    // Checkout if we are right
    // MatrixXd NewDataMatrix, NewDataMatrixTransposed;
    // NewDataMatrix = eigenvectors * Transformed.transpose();
    // NewDataMatrixTransposed = NewDataMatrix.transpose();
    // std::cout << "Check::Transformed back data matrix:\n"
    //          << NewDataMatrixTransposed << std::endl;
    X = Xsimple;
    std::cout << "Print final matrix\n" << X << std::endl;
  }

  void RunUncenteredTrickSVDWithReductionOfComponents() {
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    auto s = svd.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svd.matrixV() << std::endl;
    double totalvar = 0;
    int i = 0;
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
    for (unsigned int i = 0; i < rvec.cols(); i++) {
      if (normalise) {
        double norm = rvec.col(i).norm();
        rvec.col(i) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValues.push_back(std::make_pair(s(i), rvec.col(i)));
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
    // Transformed = X * eigenvectors;
    // Varince based selection (< 95 %)
    while (totalvar <= 0.85) {
      s(i) = fEigenValues[i].first;
      c += s(i);
      // cumulative(i) = c;
      rvec.col(i) = fEigenValues[i].second;
      totalvar = totalvar + (s(i) / s.sum());
      ++i;
    }
    // auto rveccut = rvec;
    // rveccut.conservativeResize(rveccut.rows(), i);
    // MatrixXd scut = s.asDiagonal();
    // scut.conservativeResize(scut.rows(), i);
    // scut.conservativeResize(scut.cols(), i);
    // Print();
    // std::cout << "Right vectors matrix cutted: \n" << rveccut << std::endl;
    // std::cout << "Singular values: \n" << scut << std::endl;
    MatrixXd eig = s.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "REVERSE SVD: " << std::endl;
    int j = 0;
    MatrixXd Xnew, Y;
    Xnew.conservativeResize(Xnew.rows(), X.rows());
    Xnew.conservativeResize(Xnew.cols(), X.cols());
    Y.conservativeResize(Y.rows(), X.rows());
    Y.conservativeResize(Y.cols(), X.cols());
    while (j != i) {
      Xnew =
          /*std::sqrt(X.rows())*/ s(j) * lvec.col(j) * rvec.transpose().row(j);
      Y += Xnew;
      j++;
      std::cout << "Iteration:\n " << j << std::endl;
      std::cout << "Matrix:\n " << Xnew << std::endl;
    }
    // eigenvectors.conservativeResize(eigenvectors.rows(), i);
    // Transformed.conservativeResize(Transformed.rows(), i);
    // TransformedCentered.conservativeResize(Transformed.rows(), i);
    // std::cout << "Reduced eigenvectors:\n" << eigenvectors << std::endl;
    // std::cout << "Reduced tranformed matrix \n" << Transformed << std::endl;
    // std::cout << "Total number of components to be used in transformed
    // matrix: "
    //          << i << std::endl;
    ////Transformed matrix
    // MatrixXd NewDataMatrix, NewDataMatrixTransposed;
    // NewDataMatrix = eigenvectors * Transformed.transpose();
    // NewDataMatrixTransposed = NewDataMatrix.transpose();
    // std::cout << "Transformed data matrix:\n"
    //          << NewDataMatrixTransposed << std::endl;
    // X = NewDataMatrixTransposed .array().abs();
    // std::cout << "Done." << std::endl;
    X = Y.array().abs();
    std::cout << "New matrix:\n" << X << std::endl;
  }

  void RunUncenteredTrickSVD_95Eigen() {
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    auto s = svd.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svd.matrixV() << std::endl;
    double totalvar, totalvarY, totalvarYnew = 0;
    int i = 0;
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;
    for (unsigned int i = 0; i < rvec.cols(); i++) {
      if (normalise) {
        double norm = rvec.col(i).norm();
        rvec.col(i) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValues.push_back(std::make_pair(s(i), rvec.col(i)));
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
    // variance based selection (< 95 %)
    while (totalvar <= 0.80) {
      s(i) = fEigenValues[i].first;
      c += s(i);
      // cumulative(i) = c;
      rvec.col(i) = fEigenValues[i].second;
      totalvar = totalvar + (s(i) / s.sum());
      ++i;
    }
    MatrixXd eig = s.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "REVERSE SVD: " << std::endl;
    int j = 0;
    MatrixXd Y, Ynew, Vecnew, rvecreduc, rvecreducSmall, red, redSmall, ZeroMatrixSmall, ZeroMatrixBig;
    // New matrix
    Y.conservativeResize(X.rows(), X.cols());
    Ynew.conservativeResize(X.rows(), X.cols());
    // Cutted matrix
    Vecnew.conservativeResize(rvec.rows(), rvec.cols() - i);
    Vecnew = rvec.rightCols(rvec.cols() - i); 
    // Reduction matrixes
    rvecreduc.conservativeResize(rvec.rows(), rvec.cols());
    rvecreducSmall.conservativeResize(rvec.rows(), rvec.cols());
    red = rvec;
    red.conservativeResize(red.rows(), i);
    // Reduction matrixes
    redSmall = rvec;
    redSmall = rvec.rightCols(rvec.cols() - i);
    // Zero matrixies 
    ZeroMatrixSmall.conservativeResize(rvec.rows(), rvec.cols() - i);
    ZeroMatrixSmall.fill(0);
    ZeroMatrixBig.conservativeResize(rvec.rows(), i);
    ZeroMatrixBig.fill(0);
    rvecreduc << red, ZeroMatrixSmall;
    rvecreducSmall << ZeroMatrixBig, Vecnew;
    std::cout << "Printing main matrix..." << std::endl;
    Y = X * rvecreduc * rvecreduc.transpose();
    // X = Y.array().abs();
    ////////////////////////////////////////////////////

    std::cout << "REVERSE SVD: " << std::endl;
    JacobiSVD<MatrixXd> svdY(Y, ComputeThinU | ComputeThinV);
    auto snewY = svdY.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svdY.singularValues() << std::endl;
    auto lvecnewY = svdY.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdY.matrixU() << std::endl;
    auto rvecnewY = svdY.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdY.matrixV() << std::endl;
    std::vector<std::pair<double, VectorXd>> fEigenValuesnewY;
    double cnewY = 0.0;
    int iy = 0;
    for (int iy = 0; iy < rvecnewY.cols(); iy++) {
      if (normalise) {
        double norm = rvecnewY.col(iy).norm();
        rvecnewY.col(iy) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValuesnewY.push_back(std::make_pair(snewY(iy), rvecnewY.col(iy)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(
        fEigenValuesnewY.begin(), fEigenValuesnewY.end(),
        [](std::pair<double, VectorXd> &a, std::pair<double, VectorXd> &b) {
          if (a.first > b.first)
            return true;
          if (a.first == b.first)
            return a.first > b.first;
          return false;
        });
    for (int iy = 0; iy < rvecnewY.cols(); iy++) {
      snewY(iy) = fEigenValuesnewY[iy].first;
      rvecnewY.col(iy) = fEigenValuesnewY[iy].second;
    }
    // variance based selection (< 95 %)
    std::cout << iy << std::endl;
    std::cout << snewY <<std::endl;
    while (totalvarY <= 0.80) {
      snewY(iy) = fEigenValuesnewY[iy].first;
      cnewY += snewY(iy);
      // cumulative(iy) = c;
      rvecnewY.col(iy) = fEigenValuesnewY[iy].second;
      totalvarY = totalvarY + (snewY(iy) / snewY.sum());
      ++iy;
    }
    MatrixXd eignewY = snewY.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    int jY = 0;
    MatrixXd XnewItY, YItY;
    XnewItY.conservativeResize(XnewItY.rows(), X.rows());
    XnewItY.conservativeResize(XnewItY.cols(), X.cols());
    YItY.conservativeResize(YItY.rows(), X.rows());
    YItY.conservativeResize(YItY.cols(), X.cols());
    while (jY != iy) {
      XnewItY =
          /*std::sqrt(X.rows())*/ snewY(jY) * lvecnewY.col(jY) *
          rvecnewY.transpose().row(jY);
      YItY += XnewItY;
      jY++;
      std::cout << "Iteration:\n " << jY << std::endl;
      std::cout << "Matrix:\n " << XnewItY << std::endl;
    }
    ////////////////////////////////////////////////////
    std::cout << "Printing side matrix..." << std::endl;
    Ynew = X * rvecreducSmall*rvecreducSmall.transpose();
    ////////////////////////////////////////////////////
    std::cout << "REVERSE SVD Ynew: " << std::endl;
    JacobiSVD<MatrixXd> svdYnew(Ynew, ComputeThinU | ComputeThinV);
    auto snewYnew = svdYnew.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svdYnew.singularValues() << std::endl;
    auto lvecnewYnew = svdYnew.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdYnew.matrixU() << std::endl;
    auto rvecnewYnew = svdYnew.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdYnew.matrixV() << std::endl;
    std::vector<std::pair<double, VectorXd>> fEigenValuesnewYnew;
    double cnewYnew = 0.0;
    int it = 0;
    for (int it = 0; it < rvecnewYnew.cols(); it++) {
      if (normalise) {
        double normYnew = rvecnewYnew.col(it).norm();
        rvecnewYnew.col(it) /= normYnew;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValuesnewYnew.push_back(std::make_pair(snewYnew(i), rvecnewYnew.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(
        fEigenValuesnewYnew.begin(), fEigenValuesnewYnew.end(),
        [](std::pair<double, VectorXd> &a, std::pair<double, VectorXd> &b) {
          if (a.first > b.first)
            return true;
          if (a.first == b.first)
            return a.first > b.first;
          return false;
        });
    for (unsigned int it = 0; it < rvecnewYnew.cols(); it++) {
      snewYnew(it) = fEigenValuesnewYnew[it].first;
      rvecnewYnew.col(it) = fEigenValuesnewYnew[it].second;
    }
    // variance based selection (< 95 %)
    while (totalvarYnew <= 0.50) {
      snewYnew(it) = fEigenValuesnewYnew[it].first;
      cnewYnew += snewYnew(it);
      // cumulative(i) = c;
      rvecnewYnew.col(it) = fEigenValuesnewYnew[it].second;
      totalvarYnew = totalvarYnew + (snewYnew(it) / snewYnew.sum());
      ++it;
    }
    MatrixXd eignewYnew = snewYnew.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    int jynew = 0;
    MatrixXd XnewItYnew, YItYnew;
    XnewItYnew.conservativeResize(XnewItYnew.rows(), X.rows());
    XnewItYnew.conservativeResize(XnewItYnew.cols(), X.cols());
    YItYnew.conservativeResize(YItYnew.rows(), X.rows());
    YItYnew.conservativeResize(YItYnew.cols(), X.cols());
    while (jynew != it) {
      XnewItYnew =
          /*std::sqrt(X.rows())*/ snewYnew(jynew) * lvecnewYnew.col(jynew) *
          rvecnewYnew.transpose().row(jynew);
      YItYnew += XnewItYnew;
      jynew++;
      std::cout << "Iteration:\n " << jynew << std::endl;
      std::cout << "Matrix:\n " << XnewItYnew << std::endl;
    }
    // Colecting matrix together
    X = Y + Ynew;
    std::cout << "New matrix:\n" << X << std::endl;
  }
  
 void RunUncenteredTrickSVDDiff() {
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    auto s = svd.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svd.matrixV() << std::endl;
    double totalvar, totalvarY, totalvarYnew = 0;
    int i = 0;
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;
    for (unsigned int i = 0; i < rvec.cols(); i++) {
      if (normalise) {
        double norm = rvec.col(i).norm();
        rvec.col(i) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValues.push_back(std::make_pair(s(i), rvec.col(i)));
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
    // variance based selection (< 95 %)
    while (totalvar <= 0.80) {
      s(i) = fEigenValues[i].first;
      c += s(i);
      // cumulative(i) = c;
      rvec.col(i) = fEigenValues[i].second;
      totalvar = totalvar + (s(i) / s.sum());
      ++i;
    }
    MatrixXd eig = s.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "REVERSE SVD: " << std::endl;
    int j = 0;
    MatrixXd Y, Ynew, Vecnew, rvecreduc, rvecreducSmall, red, redSmall, ZeroMatrixSmall, ZeroMatrixBig;
    // New matrix
    Y.conservativeResize(X.rows(), X.cols());
    Ynew.conservativeResize(X.rows(), X.cols());
    // Cutted matrix
    Vecnew.conservativeResize(rvec.rows(), rvec.cols() - i);
    Vecnew = rvec.rightCols(rvec.cols() - i); 
    // Reduction matrixes
    rvecreduc.conservativeResize(rvec.rows(), rvec.cols());
    rvecreducSmall.conservativeResize(rvec.rows(), rvec.cols());
    red = rvec;
    red.conservativeResize(red.rows(), i);
    // Reduction matrixes
    redSmall = rvec;
    redSmall = rvec.rightCols(rvec.cols() - i);
    // Zero matrixies 
    ZeroMatrixSmall.conservativeResize(rvec.rows(), rvec.cols() - i);
    ZeroMatrixSmall.fill(0);
    ZeroMatrixBig.conservativeResize(rvec.rows(), i);
    ZeroMatrixBig.fill(0);
    rvecreduc << red, ZeroMatrixSmall;
    rvecreducSmall << ZeroMatrixBig, Vecnew;
    std::cout << "Printing main matrix..." << std::endl;
    Y = X * rvecreduc * rvecreduc.transpose();
    JacobiSVD<MatrixXd> svdY(rvecreducSmall, ComputeThinU | ComputeThinV);
    auto snewY = svdY.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svdY.singularValues() << std::endl;
    auto lvecnewY = svdY.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdY.matrixU() << std::endl;
    auto rvecnewY = svdY.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdY.matrixV() << std::endl;
    std::vector<std::pair<double, VectorXd>> fEigenValuesnewY;
    double cnewY = 0.0;
    int iy = 0;
    for (int iy = 0; iy < rvecnewY.cols(); iy++) {
      if (normalise) {
        double norm = rvecnewY.col(iy).norm();
        rvecnewY.col(iy) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValuesnewY.push_back(std::make_pair(snewY(iy), rvecnewY.col(iy)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(
        fEigenValuesnewY.begin(), fEigenValuesnewY.end(),
        [](std::pair<double, VectorXd> &a, std::pair<double, VectorXd> &b) {
          if (a.first > b.first)
            return true;
          if (a.first == b.first)
            return a.first > b.first;
          return false;
        });
    for (int iy = 0; iy < rvecnewY.cols(); iy++) {
      snewY(iy) = fEigenValuesnewY[iy].first;
      rvecnewY.col(iy) = fEigenValuesnewY[iy].second;
    }
    // variance based selection (< 95 %)
    std::cout << iy << std::endl;
    std::cout << snewY <<std::endl;
    while (totalvarY <= 0.80) {
      snewY(iy) = fEigenValuesnewY[iy].first;
      cnewY += snewY(iy);
      // cumulative(iy) = c;
      rvecnewY.col(iy) = fEigenValuesnewY[iy].second;
      totalvarY = totalvarY + (snewY(iy) / snewY.sum());
      ++iy;
    }
    MatrixXd eignewY = snewY.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    int jY = 0;
    MatrixXd XnewItY, YItY;
    XnewItY.conservativeResize(XnewItY.rows(), X.rows());
    XnewItY.conservativeResize(XnewItY.cols(), X.cols());
    YItY.conservativeResize(YItY.rows(), X.rows());
    YItY.conservativeResize(YItY.cols(), X.cols());
    while (jY != iy) {
      XnewItY =
          /*std::sqrt(X.rows())*/ snewY(jY) * lvecnewY.col(jY) *
          rvecnewY.transpose().row(jY);
      YItY += XnewItY;
      jY++;
      std::cout << "Iteration:\n " << jY << std::endl;
      std::cout << "Matrix:\n " << XnewItY << std::endl;
    }
    ////////////////////////////////////////////////////
    std::cout << "Printing side matrix..." << std::endl;
    Ynew = X * rvecreducSmall*rvecreducSmall.transpose();
    ////////////////////////////////////////////////////
    
    // Colecting matrix together
    X = Y + Ynew;
    std::cout << "New matrix:\n" << X << std::endl;
  }

  void RunUncenteredTrickSVDWithReductionOfComponentsIteration() {
    JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
    auto s = svd.singularValues();
    std::cout << "Its singular values matrix X are:" << std::endl
              << svd.singularValues() << std::endl;
    auto lvec = svd.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix X:"
        << std::endl
        << svd.matrixU() << std::endl;
    auto rvec = svd.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix X:"
        << std::endl
        << svd.matrixV() << std::endl;
    double totalvar = 0;
    int i = 0;
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
    for (unsigned int i = 0; i < rvec.cols(); i++) {
      if (normalise) {
        double norm = rvec.col(i).norm();
        rvec.col(i) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValues.push_back(std::make_pair(s(i), rvec.col(i)));
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
    // variance based selection (< 95 %)
    /*
    while (totalvar <= 0.95) {
      s(i) = fEigenValues[i].first;
      c += s(i);
      // cumulative(i) = c;
      rvec.col(i) = fEigenValues[i].second;
      totalvar = totalvar + (s(i) / s.sum());
      ++i;
    }
    */
    // auto rveccut = rvec;
    // rveccut.conservativeResize(rveccut.rows(), i);
    // MatrixXd scut = s.asDiagonal();
    // scut.conservativeResize(scut.rows(), i);
    // scut.conservativeResize(scut.cols(), i);
    // Print();
    // std::cout << "Right vectors matrix cutted: \n" << rveccut << std::endl;
    // std::cout << "Singular values: \n" << scut << std::endl;
    MatrixXd eig = s.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    int f = 1;
    MatrixXd Xnew, X1new, Y;
    Xnew.conservativeResize(X.rows(), X.rows());
    // Xnew.conservativeResize(Xnew.cols(), X.cols());
    X1new.conservativeResize(X.rows(), X.rows());
    // X1new.conservativeResize(X1new.cols(), X.cols());
    Y.conservativeResize(X.rows(), X.rows());
    // Y.conservativeResize(Y.cols(), X.cols());
    while (f != X.cols()) {
      Xnew =
          /*std::sqrt(X.rows())*/ s(f) * lvec.col(f) * rvec.transpose().row(f);
      Y += Xnew;
      f++;
      std::cout << "Iteration:\n " << f << std::endl;
      std::cout << "Matrix:\n " << Xnew << std::endl;
    }
    X1new = s(0) * lvec.col(0) * rvec.transpose().row(0);
    std::cout << "REVERSE SVD: " << std::endl;
    JacobiSVD<MatrixXd> svdnew(Y, ComputeThinU | ComputeThinV);
    auto snew = svdnew.singularValues();
    std::cout << "Its singular values are:" << std::endl
              << svdnew.singularValues() << std::endl;
    auto lvecnew = svdnew.matrixU();
    std::cout
        << "Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdnew.matrixU() << std::endl;
    auto rvecnew = svdnew.matrixV();
    std::cout
        << "Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdnew.matrixV() << std::endl;
    std::vector<std::pair<double, VectorXd>> fEigenValuesnew;
    double cnew = 0.0;
    for (unsigned int i = 0; i < rvecnew.cols(); i++) {
      if (normalise) {
        double norm = rvecnew.col(i).norm();
        rvecnew.col(i) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValuesnew.push_back(std::make_pair(snew(i), rvecnew.col(i)));
    }
    // Sorting Eigen pairs [eigenvalue, eigenvector]
    std::sort(
        fEigenValuesnew.begin(), fEigenValuesnew.end(),
        [](std::pair<double, VectorXd> &a, std::pair<double, VectorXd> &b) {
          if (a.first > b.first)
            return true;
          if (a.first == b.first)
            return a.first > b.first;
          return false;
        });
    for (unsigned int i = 0; i < rvecnew.cols(); i++) {
      snew(i) = fEigenValuesnew[i].first;
      rvec.col(i) = fEigenValues[i].second;
    }
    // variance based selection (< 95 %)
    while (totalvar <= 0.95) {
      snew(i) = fEigenValuesnew[i].first;
      cnew += snew(i);
      // cumulative(i) = c;
      rvecnew.col(i) = fEigenValuesnew[i].second;
      totalvar = totalvar + (snew(i) / snew.sum());
      ++i;
    }
    MatrixXd eignew = snew.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    int j = 0;
    MatrixXd XnewIt, YIt;
    XnewIt.conservativeResize(XnewIt.rows(), X.rows());
    XnewIt.conservativeResize(XnewIt.cols(), X.cols());
    YIt.conservativeResize(YIt.rows(), X.rows());
    YIt.conservativeResize(YIt.cols(), X.cols());
    while (j != i) {
      XnewIt =
          /*std::sqrt(X.rows())*/ snew(j) * lvecnew.col(j) *
          rvecnew.transpose().row(j);
      YIt += XnewIt;
      j++;
      std::cout << "Iteration:\n " << j << std::endl;
      std::cout << "Matrix:\n " << XnewIt << std::endl;
    }
    // eigenvectors.conservativeResize(eigenvectors.rows(), i);
    // Transformed.conservativeResize(Transformed.rows(), i);
    // TransformedCentered.conservativeResize(Transformed.rows(), i);
    // std::cout << "Reduced eigenvectors:\n" << eigenvectors << std::endl;
    // std::cout << "Reduced tranformed matrix \n" << Transformed << std::endl;
    // std::cout << "Total number of components to be used in transformed
    // matrix: "
    //          << i << std::endl;
    ////Transformed matrix
    // MatrixXd NewDataMatrix, NewDataMatrixTransposed;
    // NewDataMatrix = eigenvectors * Transformed.transpose();
    // NewDataMatrixTransposed = NewDataMatrix.transpose();
    // std::cout << "Transformed data matrix:\n"
    //          << NewDataMatrixTransposed << std::endl;
    // X = NewDataMatrixTransposed .array().abs();
    // std::cout << "Done." << std::endl;
    X = YIt.array().abs();
    X = X + X1new;
    std::cout << "New matrix after all transformations:\n" << X << std::endl;
  }

    void RunUncenteredTrickSVDCutApproach() {
    double totalvar = 0;
    int i = 0;
    JacobiSVD<MatrixXd> svdX(X, ComputeThinU | ComputeThinV);
    auto s = svdX.singularValues();
    std::cout << "X Its singular values are:" << std::endl
              << svdX.singularValues() << std::endl;
    auto lvec = svdX.matrixU();
    std::cout
        << "X Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdX.matrixU() << std::endl;
    auto rvec = svdX.matrixV();
    std::cout
        << "X Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdX.matrixV() << std::endl;
    ///////////////////////////////////////////////////////////////////////////////
    //Checking again...
    ///////////////////////////////////////////////////////////////////////////////
    MatrixXd diag = s.asDiagonal();
    MatrixXd recheck = lvec*diag*rvec.transpose();
    MatrixXd  diff = recheck - X;
    std::cout << "diff:\n" << diff.array().abs().sum() << "\n";
    std::cout << "Print this thing: X = USV* :" << recheck << std::endl;
    /////////////////////////////////////////////////////////////////////////////// 
    // SVD transformation (hardcodded for two eigenvalues)
    //////////////////////////////////////////////////////////////////////////////
    MatrixXd TransformedU = X*rvec;
    MatrixXd Xreduced = X - s(1)*lvec.col(1)*rvec.transpose().row(1) - s(2)*lvec.col(2)*rvec.transpose().row(2);
    //////////////////////////////////////////////////////////////////////////////
    JacobiSVD<MatrixXd> svdredX(Xreduced, ComputeThinU | ComputeThinV);
    auto sred = svdredX.singularValues();
    std::cout << "Xreduced.. Its singular values are:" << std::endl
              << svdredX.singularValues() << std::endl;
    auto lvecred = svdredX.matrixU();
    std::cout
        << "Xreduced.. Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdredX.matrixU() << std::endl;
    auto rvecred = svdredX.matrixV();
    std::cout
        << "Xreduced.. Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdredX.matrixV() << std::endl;
    std::vector<std::pair<double, VectorXd>> fEigenValues;
    double c = 0.0;
    for (unsigned int i = 0; i < rvecred.cols(); i++) {
      if (normalise) {
        double norm = rvecred.col(i).norm();
        rvecred.col(i) /= norm;
        std::cout << "Vectors are normalised." << std::endl;
      }
      fEigenValues.push_back(std::make_pair(sred(i), rvecred.col(i)));
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
    for (unsigned int i = 0; i < rvecred.cols(); i++) {
      sred(i) = fEigenValues[i].first;
      rvecred.col(i) = fEigenValues[i].second;
    }
    // Printing current state information before  PC cutoff
    std::cout << "Printing original information after PCA" << std::endl;
    // Transformed = X * eigenvectors;
    // Varince based selection (< 95 %)
    while (totalvar <= 0.85) {
      sred(i) = fEigenValues[i].first;
      c += sred(i);
      // cumulative(i) = c;
      rvecred.col(i) = fEigenValues[i].second;
      totalvar = totalvar + (sred(i) / sred.sum());
      ++i;
    }
    // auto rveccut = rvec;
    // rveccut.conservativeResize(rveccut.rows(), i);
    // MatrixXd scut = s.asDiagonal();
    // scut.conservativeResize(scut.rows(), i);
    // scut.conservativeResize(scut.cols(), i);
    // Print();
    // std::cout << "Right vectors matrix cutted: \n" << rveccut << std::endl;
    // std::cout << "Singular values: \n" << scut << std::endl;
    MatrixXd eig = sred.asDiagonal();
    std::cout << "---------------------------\n" << std::endl;
    std::cout << "REVERSE SVD: " << std::endl;
    int j = 0;
    MatrixXd Xnew, Y;
    Xnew.conservativeResize(Xnew.rows(), X.rows());
    Xnew.conservativeResize(Xnew.cols(), X.cols());
    Y.conservativeResize(Y.rows(), X.rows());
    Y.conservativeResize(Y.cols(), X.cols());
    while (j != i) {
      Xnew =
          /*std::sqrt(X.rows())*/sred(j)*lvecred.col(j)*rvecred.transpose().row(j);
      Y += Xnew;
      j++;
      std::cout << "Iteration:\n " << j << std::endl;
      std::cout << "Matrix:\n " << Xnew << std::endl;
    }
    Y = Y + /*std::sqrt(X.rows())*/s(1)*lvec.col(1)*rvec.transpose().row(1) + /*std::sqrt(X.rows())*/s(2)*lvec.col(2)*rvec.transpose().row(2);
    X = Y.array().abs();
    std::cout << "New matrix:\n" << X << std::endl;
  }

void RunUncenteredTrickSVD2() {
    double totalvar = 0;
    int i = 0;
    JacobiSVD<MatrixXd> svdX(X, ComputeThinU | ComputeThinV);
    auto s = svdX.singularValues();
    std::cout << "X Its singular values are:" << std::endl
              << svdX.singularValues() << std::endl;
    auto lvec = svdX.matrixU();
    std::cout
        << "X Its left singular vectors are the columns of the thin U matrix:"
        << std::endl
        << svdX.matrixU() << std::endl;
    auto rvec = svdX.matrixV();
    std::cout
        << "X Its right singular vectors are the columns of the thin V matrix:"
        << std::endl
        << svdX.matrixV() << std::endl;
    ///////////////////////////////////////////////////////////////////////////////
    //Checking again...
    ///////////////////////////////////////////////////////////////////////////////
    MatrixXd diag = s.asDiagonal();
    MatrixXd recheck = lvec*diag*rvec.transpose();
    MatrixXd  diff = recheck - X;
    std::cout << "diff:\n" << diff.array().abs().sum() << "\n";
    std::cout << "Print this thing: X = USV* :" << recheck << std::endl;
    /////////////////////////////////////////////////////////////////////////////// 
    // SVD transformation (hardcodded for two eigenvalues)
    //////////////////////////////////////////////////////////////////////////////
    MatrixXd TransformedU = X*rvec;
    MatrixXd Y =  s(1)*lvec.col(1)*rvec.transpose().row(1) - s(2)*lvec.col(2)*rvec.transpose().row(2);
    X = Y.array().abs();
    std::cout << "New matrix:\n" << X << std::endl;
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
