//===--- TGenes.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TGenes.h
 * @brief Implementation of TGenes for LibGA
 * prototype
 */

#pragma once

#ifndef __TGENES__
#define __TGENES__

#include <functional>
#include <memory>
#include <vector>

#include <hwloc.h>
#include <unistd.h>

#include "hpc/GAEvaluate.h"
#include "hpc/GASequentualEvaluator.h"
#include "hpc/GASimpleEvaluator.h"

#include "instrumentation/CPUManager.h"

namespace geantvmoop {

template <typename F> class TGenes;
template <typename F> using individual_t = std::shared_ptr<TGenes<F> >;

template <typename F> class TGenes {

private:
  typename F::Input input;
  typename F::Output output;

public:
  TGenes() {}

  TGenes(const typename F::Input &i, bool fEvaluated = true) : input(i) {
    if (fEvaluated) {
      CPUManager cpumgr;
      cpumgr.InitCPU();
      hwloc_topology_t topology;
      int nbcores, ccores;
      unsigned int microseconds;
      microseconds = 3000;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      // printf("%d cores\n", nbcores);
      hwloc_topology_destroy(topology);
      ccores =
          nbcores - cpumgr.GetCurrentValue() / 100 * nbcores; // just a test
      if (ccores < 0.3) {
        usleep(microseconds);
      } else {
        Evaluate();
        //}
        // GASimpleEvaluator::GAEvaluate();
      }
    }
  }

  ~TGenes() {}

  bool IsDominating(const TGenes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (output[i] > other.output[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsDominated(const TGenes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (output[i] < other.output[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsEqual(const TGenes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (other.output[i] != output[i])
        return false;
    }
    return true;
  }

  void Evaluate() {
    // output = GASequentualEvaluator::Evaluate();
#ifdef FORK
      output = F::Evaluate(input);
#else
      output = F::Evaluate(input);
#endif
  }

  const typename F::Input &GetInput() const { return input; }

  const typename F::Output &GetOutput() const { return output; }

  friend std::ostream &operator<<(std::ostream &os, const TGenes<F> &ind) {
    auto indvector = ind.GetInput();
    for (int i = 0; i < indvector.size(); ++i)
      os << indvector[i] << " ";
    return os;
  }
};
}

#endif
