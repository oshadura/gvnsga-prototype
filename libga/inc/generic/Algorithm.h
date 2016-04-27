#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <iostream>
#include "PF.h"

template <typename DerivedClass, typename Trait> class Algorithm {

private:
  Trait problem;

public:

  Algorithm(Trait problem) : problem(problem) {}

  /*
  virtual ParetoFront<Trait> solve(std::ostream *sInfo = nullptr,
                                   std::ostream *sOut = nullptr) {
    initialize();
    for (int i = 0; i < maxGeneration; ++i) {
      if (sOut != nullptr)
        front().json(jsonFront);
      next();
      if (sInfo != nullptr) {
        *sInfo << "generation: " << i << " -> ";
        info(*sInfo);
      }
    }
    if (sOut != nullptr)
      *sOut << jsonFront;
    return front();
  }
  */

  Trait getProblem() const { return problem; }

  void setProblem(Trait problem) { this->problem = problem; }

  static void waitForKey() {
    do {
      std::cout << '\n' << "Press a key to continue...";
    } while (std::cin.get() != '\n');
  }

  void next() { return static_cast<DerivedClass *>(this)->next_(); }

  void initialize() { return static_cast<DerivedClass *>(this)->init_(); }

  void info(std::ostream &os) {
    return static_cast<DerivedClass *>(this)->info_(os);
  }

  PF<Trait> front() {
    return static_cast<DerivedClass *>(this)->front_();
  }
};

#endif