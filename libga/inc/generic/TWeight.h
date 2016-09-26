//===--- TWeight.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TWeight.h
 * @brief Implementation of weights for LibGA
 * prototype
 */
//
#pragma once

#ifndef __TWEIGHTS__
#define __TWEIGHTS__

#include "tools/GASort.h"
#include "tools/Random.h"
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace geantvmoop {

class Weights : public std::vector<double> {

public:
  Weights(std::initializer_list<double> list) : std::vector<double>(list) {}

  ~Weights() {}

  static std::vector<Weights> GetWeights(int n) {
    auto RandomWeight = []() {
      auto r = Random::GetInstance();
      std::vector<double> v = {r.RandomDouble(), r.RandomDouble(),
                               r.RandomDouble()};
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
          weights.push_back(w);
        }
      }
    }
    return weights;
  }

  double GetDistance(const Weights &w) {
    double d = 0;
    if (this->size() != w.size())
      throw std::runtime_error("Weights has not the same length!");
    for (unsigned int i = 0; i < w.size(); ++i)
      d += ((*this)[i] - w[i]) * ((*this)[i] - w[i]);
    return std::sqrt(d);
  }

  std::vector<double> GetDistanceAll(const std::vector<Weights> &w) {
    std::vector<double> distance;
    for (unsigned int i = 0; i < w.size(); ++i) {
      distance.push_back(GetDistance(w[i]));
    }
    return distance;
  }

  std::vector<int> GetNearestNeighbor(const std::vector<Weights> &w,
                                      unsigned int NNearest) {
    std::vector<double> distance = GetDistanceAll(w);
    auto index = SortUtil::GetIndex(w.size());
    std::sort(index.begin(), index.end(),
              [&distance](const int &lhs, const int &rhs) {
                return distance[lhs] < distance[rhs];
              });
    while (index.size() > NNearest)
      index.pop_back();
    return index;
  }

  static double GetWSum(const Weights &w, const std::vector<double> &output) {
    double sum = 0;
    for (unsigned int i = 0; i < w.size(); ++i) {
      sum += w[i] * output[i];
    }
    return sum;
  }
};
}

#endif
