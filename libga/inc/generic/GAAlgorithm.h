//===--- GAAlgorithm.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GAAlgorithm.h
 * @brief Implementation of generic interface for algorithm class for LibGA
 * prototype
 */
//

#ifndef __GAALGORITHM__
#define __GAALGORITHM__

#include <iostream>
#include "PF.h"

namespace geantvmoop {

template <typename Derived, typename F> class GAAlgorithm {

private:
  F problem;

public:
  int fMaxGeneration = 100;

  GAAlgorithm(F problem) : problem(problem) {}

  ~GAAlgorithm() {}

  virtual PF<F> SolvePF() {
    Initialize();
    for (int i = 0; i < fMaxGeneration; ++i) {
      Evolution();
    }
    return GetParetoFront();
  }

  F GetProblem() const { return problem; }

  void SetProblem(F problem) { this->problem = problem; }

  void Evolution() { return static_cast<Derived *>(this)->EvolutionImpl(); }

  void Initialize() { return static_cast<Derived *>(this)->InitializeImpl(); }

  void Print(std::ostream &os) {
    return static_cast<Derived *>(this)->PrintImpl(os);
  }

  PF<F> GetParetoFront() {
    return static_cast<Derived *>(this)->GetParetoFrontImpl();
  }
};
}

#endif
