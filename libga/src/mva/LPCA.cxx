#include "mva/LPCA.h"

namespace geantvmoop {

void LPCA::LoadData(const char *data, char sep) {
  // Read data
  unsigned int row = 0;
  // std::ifstream file(data);
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

void LPCA::RunLPCA() {
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
    eigen_pairs.push_back(std::make_pair(eigenvalues(i), eigenvectors.col(i)));
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

void LPCA::RunLPCAWithReductionOfComponents() {
  double totalvar = 0;
  int i = 0;
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
    eigen_pairs.push_back(std::make_pair(eigenvalues(i), eigenvectors.col(i)));
  }
  // Sorting Eigen pairs [eigenvalue, eigenvector]
  std::sort(eigen_pairs.begin(), eigen_pairs.end(),
            [](const std::pair<double, VectorXd> a,
               const std::pair<double, VectorXd> b)
                ->bool { return (a.first > b.first); });
  // Printing current state information before  PC cutoff
  std::cout << "Printing original information after PCA" << std::cout;
  Transformed = X * eigenvectors;
  //TransformedCentered = Xcentered * eigenvectors;
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
  eigenvectors.conservativeResize(eigenvectors.rows(), i);
  Transformed.conservativeResize(Transformed.rows(), i);
  //TransformedCentered.conservativeResize(Transformed.rows(), i);
  std::cout << "Reduced eigenvectors:\n" << eigenvectors << std::endl;
  std::cout << "Total number of components to be used in Transformed matrix: "
            << i << std::endl;
  // Transformed matrix
  MatrixXd NewDataMatrix, NewDataMatrixCentered;
  NewDataMatrix = eigenvectors * Transformed.transpose();
  //NewDataMatrixCentered = eigenvectors * TransformedCentered.transpose();
  std::cout << "New Transformed data:\n" << NewDataMatrix << std::endl;
  //std::cout << "New Transformed (centered?) data:\n" << NewDataMatrixCentered
  //          << std::endl;
}

void LPCA::Print() {
  std::cout << "Input data:\n" << X << std::endl;
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
  std::cout << "Transformed centred data:\n" << TransformedCentered
            << std::endl;
}

void LPCA::WriteTransformed(std::string file) {
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

void LPCA::WriteEigenvectors(std::string file) {
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
}
