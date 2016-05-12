#include "algorithms/MOEADWeights.h"

#include <algorithm>

namespace geantvmoop{

// static
std::vector<Weights> Weights::GetWeights(int n) {
  auto rndWeight = []() {
    // Checkout correct values!
    // auto r = Random::getInstance();
    // std::vector<double> v = { r->rndDouble(), r->rndDouble(), r->rndDouble()
    // };
    std::vector<double> v;
    v.resize(3);
    //
    double sum = 0;
    for (auto it = v.begin(); it != v.end(); ++it)
      sum += *it;
    for (unsigned int i = 0; i < v.size(); ++i)
      v[i] /= sum;

    return Weights{v[0], v[1], v[2]};
  };

  std::vector<Weights> weights;
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      if (i + j <= n) {
        int k = n - i - j;
        std::vector<double> v(3);
        v[0] = i / (double)n;
        v[1] = j / (double)n;
        v[2] = k / (double)n;
        Weights w{v[0], v[1], v[2]};
        // std::cout << v[0]<< ", " << v[1]<< ", " << v[2] << std::endl;
        weights.push_back(w);
      }
    }
  }

  // for (int i = 0; i < n; ++i) weights.push_back(rndWeight());
  return weights;
}

double Weights::GetDistance(const Weights &w) {
  double d = 0;
  if (this->size() != w.size())
    throw std::runtime_error("Weights has not the same length!");
  for (unsigned int i = 0; i < w.size(); ++i)
    d += ((*this)[i] - w[i]) * ((*this)[i] - w[i]);
  return std::sqrt(d);
}

std::vector<double> Weights::GetDistanceAll(const std::vector<Weights> &w) {
  std::vector<double> d;
  for (unsigned int i = 0; i < w.size(); ++i) {
    d.push_back(GetDistance(w[i]));
  }
  return d;
}

std::vector<int>
Weights::GetNearestNeighborByIndex(const std::vector<Weights> &w,
                                   unsigned int numOfNearest) {
  std::vector<double> d = GetDistanceAll(w);
  auto index = Sort::GetIndex(w.size());
  std::sort(index.begin(), index.end(),
            [&d](const int &lhs, const int &rhs) { return d[lhs] < d[rhs]; });
  while (index.size() > numOfNearest)
    index.pop_back();
  return index;
}

// static
double Weights::GetWeightedSum(const Weights &w,
                               const std::vector<double> &output) {
  double sum = 0;
  for (unsigned int i = 0; i < w.size(); ++i) {
    sum += w[i] * output[i];
  }
  return sum;
}

}