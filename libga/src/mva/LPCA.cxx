#include "mva/LPCA.h"

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
    std::cout << "Failed to read file " << data << std::endl;
  }
}

void LPCA::UploadPopulation(Population<double> &pop) {
  for (int i = 0; i < pop.GetPopulationSize(); ++i) {
    for (int j = 0; j < pop.GetGenes(i).size(); ++j){
          auto ind = pop.GetGenes(i);
          auto gene = ind.GetGene(j);
          X.row(i) = VectorXd::Map(&gene, sizeof(gene));
    }
  }
}

void LPCA::LoadUpdatedPopulation(Population<double> &pop) {
}


void LPCA::RunLPCA() {
  Xcentered = X.rowwise() - X.colwise().mean();
  C = (Xcentered.adjoint() * Xcentered) / double(X.rows());
  EigenSolver<MatrixXd> edecomp(C);
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
    eigen_pairs.push_back(std::make_pair(eigenvalues(i), eigenvectors.col(i)));
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
  transformed = Xcentered * eigenvectors;
}

void LPCA::Print() {
  std::cout << "Input data: " << X << std::endl;
  std::cout << "Centered data: " << Xcentered << std::endl;
  std::cout << "Covariance matrix: " << C << std::endl;
  std::cout << "Eigenvalues: " << eigenvalues << std::endl;
  std::cout << "Eigenvectors: " << eigenvectors << std::endl;
  std::cout << "Sorted eigenvalues: " << std::endl;
  for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
    if (eigenvalues(i) > 0) {
      std::cout << "PC " << i + 1 << ": Eigenvalue: " << eigenvalues(i);
      printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",
             eigenvalues(i) / eigenvalues.sum(),
             cumulative(i) / eigenvalues.sum());
    }
  }
  std::cout << std::endl;
  std::cout << "Sorted eigenvectors:" << eigenvectors << std::endl;
  std::cout << "Transformed data:" << X *eigenvectors << std::endl;
  // std::cout << "Transformed centred data:"<< transformed << std::endl;
}

void LPCA::WriteTransformed(std::string file) {
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