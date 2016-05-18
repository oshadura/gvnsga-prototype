//===--- GANSGA3.h - LibGA ----------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Nsga2.h
 * @brief Implementation of NSGA-III algorithms for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//


#ifndef __NSGA3__
#define __NSGA3__

#include "generic/GAAlgorithm.h"
#include "generic/PF.h"
#include "generic/ReferencePoint.h"
#include "addstructures/GANDRank.h"
#include "addstructures/GACD.h"
#include "addstructures/GAComparator.h"
#include "tools/Random.h"
#include <iostream>

#include "gaoperators/GATournamentSelection.h"
#include "gaoperators/GASBXCrossover.h"
#include "gaoperators/GAPolMutation.h"
#include "addstructures/GAComparator.h"

namespace geantvmoop {

template <typename F> class NSGA3 : public GAAlgorithm<NSGA3<F>, F> {

public:
  NSGA3(F problem) : GAAlgorithm<NSGA3<F>, F>(problem) {}

  void Initialize() {}

  void Evolution() {}

  void Print(std::ostream &os) {}

  PF<F> GetParetoFront() { return fFront; }

private:
  Population<F> pop;

  PF<F> fFront;

};
}

#endif
