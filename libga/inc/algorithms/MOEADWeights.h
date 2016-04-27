#ifndef MOEADWEIGHTS_H
#define MOEADWEIGHTS_H

#include <vector>
#include <initializer_list>
#include <stdexcept>
#include <cmath>
#include <unordered_map>

#include "tools/Sort.h"
#include "tools/Generator.h"

class Weights: public std::vector<double>{

public:
  Weights(std::initializer_list<double> list) : std::vector<double>(list) {}
  static std::vector<Weights> GetWeights(int n);
  double GetDistance(const Weights &w);
  std::vector<double> GetDistanceAll(const std::vector<Weights> &w);
  std::vector<int> GetNearestNeighborByIndex(const std::vector<Weights> &w,
                                             unsigned int numOfNearest);
  static double GetWeightedSum(const Weights &w,
                               const std::vector<double> &output);
};

#endif
