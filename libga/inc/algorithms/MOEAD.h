#ifndef __MOED__
#define __MOED__

#include "generic/Algorithm.h"
#include "MOEADWeights.h"
#include "gaoperators/SBXCrossover.h"
#include "gaoperators/PolynomialMutation.h"

/*

template <typename Trait> class MOEAD : public Algorithm<MOEAD<Trait>, Trait> {
private:
  std::vector<Weights> weights;
  std::vector<std::vector<int> > nearest;
  std::vector<double> refPoint;

  Population<Trait> pop;
  std::vector<double> fitness;

  ParetoFront<Trait> f;

  double propMutation = 0.2;
  RandomSelection selection;

  int counter = 0;

public:
  int populationSize = 1000;
  int T = 5;

  MOEAD(Trait problem) : Algorithm<MOEAD<Trait>, Trait>(problem) {}

  void init_() {

    weights = Weights::getWeights(populationSize);
    populationSize = weights.size();
    std::cout << populationSize << std::endl;

    if (T >= populationSize)
      throw std::runtime_error("Please set T lower than population size!");

    for (auto w : weights)
      nearest.push_back(w.getNearestNeighborByIndex(weights, T));

    pop = Population<Trait>{ populationSize };
    refPoint = getReferencePoint(pop);
    fitness =
        std::vector<double>(populationSize, std::numeric_limits<double>::max());
    for (int i = 0; i < pop.size(); ++i)
      fitness[i] = getFitness(weights[i], pop[i]->getOutput());
  }

  template <typename T> int randomIndex(std::vector<T> v) {
    return Random::getInstance()->rndInt(0, v.size());
  }

  void next_() {

    counter = 0;

    for (int i = 0; i < pop.size(); ++i) {

      auto a = pop[nearest[i][randomIndex(nearest[i])]];
      auto b = pop[nearest[i][randomIndex(nearest[i])]];
      IndividualPtr<Trait> off = SBXCrossover::crossover(a, b);
      if (Random::getInstance()->rndDouble() < propMutation)
        off = PolynomialMutation::mutate(off);

      updateReferencePoint(refPoint, off);
      auto out = off->getOutput();

      for (int ithNearest : nearest[i]) {
        double value = getFitness(weights[ithNearest], out);
        if (value < fitness[ithNearest]) {
          pop[ithNearest] = off;
          fitness[ithNearest] = value;
          ++counter;
        }
      }
      f.add(off);
    }

    // for (int k = 0; k < fitness.size(); ++k) std::cout << fitness[k] <<
    // std::endl;
    // Algorithm<MOEAD<Trait>, Trait>::waitForKey();
  }

  void info_(std::ostream &os) { os << counter << std::endl; }

  ParetoFront<Trait> front_() { return f; }

  static void updateReferencePoint(std::vector<double> &ref,
                                   const IndividualPtr<Trait> &ind) {
    int numOfObjectives = Trait::getNumOfObjectives();
    auto v = ind->getOutput();
    for (int i = 0; i < numOfObjectives; ++i) {
      ref[i] = std::min(ref[i], v[i]);
    }
  }

  static std::vector<double> getReferencePoint(const Population<Trait> &pop) {
    int numOfObjectives = Trait::getNumOfObjectives();
    std::vector<double> ref(numOfObjectives);
    for (int i = 0; i < numOfObjectives; ++i) {
      auto v = pop.getObjective(i);
      ref[i] = *(std::min_element(v.begin(), v.end()));
    }
    return ref;
  }

  double getFitness(const Weights &w, const std::vector<double> &output) {
    return getTchebichew(w, output);
  }

  double getTchebichew(const Weights &w, const std::vector<double> &output) {
    double maxDistance = 0;
    for (int i = 0; i < Trait::getNumOfObjectives(); ++i) {
      maxDistance = std::max(maxDistance, w[i] * (output[i] - refPoint[i]));
    }
    return maxDistance;
  }
};

*/
#endif
