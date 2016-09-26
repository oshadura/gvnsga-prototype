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
#pragma once

#ifndef __GAALGORITHM__
#define __GAALGORITHM__

#include "PF.h"
#include <iostream>

namespace geantvmoop {

template <typename Derived, typename F> class GAAlgorithm {

private:
  F problem;

public:
  int fMaxGeneration = 100;

  GAAlgorithm(F problem) : problem(problem) {}

  ~GAAlgorithm() = default;

  virtual PF<F> SolvePF() {
    Initialize();
    for (int i = 0; i < fMaxGeneration; ++i) {
      Evolution();
      // Print(std::cout);
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
