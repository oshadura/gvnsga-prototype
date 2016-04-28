#include <vector>
#include <ostream>
#include <string>
#include <utility>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <functional>

#include "tools/Generator.h"
#include "generic/TGenes.h"

void GeneratorDouble();

void GeneratorDouble(int fMin, int fMax) {
	/*
  std::vector<double> ind;
  std::random_device rnd_device;
  auto seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 mersenne_engine(seed);
  std::uniform_real_distribution<double> dist(fMin, fMax);
  auto gen = std::bind(dist, std::ref(mersenne_engine));
  std::generate_n(std::begin(ind) + i, 1, gen());
  */
}

void GeneratorDoubleVector(Genes<double> individual, int fMin, int fMax);

void GeneratorInt();

void GeneratorInt(int fMin, int fMax);

void GeneratorIntVector(Genes<int> individual,int fMin, int fMax);

void GeneratorVector();

void GeneratorSIMD(int fMin, int fMax);

void GeneratorSIMDVector(Genes<int> individual, int fMin, int fMax);