#include "mva/KPCA.h"

namespace geantvmoop{

void KPCA::LoadData(const char *data, char sep) {
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

template <typename F> void KPCA::UploadPopulation(Population<F> &pop) {
  for (int i = 0; i < pop.GetPopulationSize(); ++i) {
    for (int j = 0; j < pop.GetGenes(i).size(); ++j){
          auto ind = pop.GetGenes(i);
          auto gene = ind.GetGene(j);
          X.row(i) = VectorXd::Map(&gene, sizeof(gene));
    }
  }
}

template <typename F> void KPCA::LoadUpdatedPopulation(Population<F> &pop) {
}


double KPCA::Kernel(const VectorXd &a, const VectorXd &b) {
  switch (kernel_type) {
  case 2:
    return (std::pow(a.dot(b) + constant, order));
  default:
    return (std::exp(-gamma * ((a - b).squaredNorm())));
  }
}

void KPCA::RunKpca() {
  // Fill kernel matrix
  K.resize(X.rows(), X.rows());
  for (unsigned int i = 0; i < X.rows(); i++) {
    for (unsigned int j = i; j < X.rows(); j++) {
      K(i, j) = K(j, i) = Kernel(X.row(i), X.row(j));
      // printf("k(%i,%i) = %f\n",i,j,K(i,j));
    }
  }
  // std::cout << std::endl << K << std::endl;
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
    eigen_pairs.push_back(std::make_pair(eigenvalues(i), eigenvectors.col(i)));
  }
  // http://stackoverflow.com/questions/5122804/sorting-with-lambda
  std::sort(eigen_pairs.begin(), eigen_pairs.end(),
            [](const std::pair<double, VectorXd> a, const std::pair<double, VectorXd> b)
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

  /*
  std::cout
 << "Input data:" << std::endl
 << X << std::endl
 << std::endl
;
  std::cout
 << "Centered data:"<< std::endl
 << Xcentered << std::endl
 << std::endl
;
  std::cout
 << "Centered kernel matrix:" << std::endl
 << Kcentered << std::endl
 << std::endl
;
  std::cout
 << "Eigenvalues:" << std::endl
 << eigenvalues << std::endl
 << std::endl
;
  std::cout
 << "Eigenvectors:" << std::endl
 << eigenvectors << std::endl
 << std::endl
;
  */
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

void KPCA::Print() {
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

void KPCA::WriteTransformed(std::string file) {
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

void KPCA::WriteEigenvectors(std::string file) {
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