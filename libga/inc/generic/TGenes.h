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
//
 #pragma once

#ifndef __TGENES__
#define __TGENES__

#include <vector>
#include <memory>
#include <functional>

namespace geantvmoop {

template <typename F> class TGenes;
template <typename F> using individual_t = std::shared_ptr<TGenes<F> >;
// TBD: Add one more templated typename for evaluation procedure (fork-,
// MPI-based..)
template <typename F> class TGenes {

private:
  typename F::Input input;
  typename F::Output output;

public:
  TGenes() {}

  TGenes(const typename F::Input &i, bool fEvaluated = true) : input(i) {
    if (fEvaluated)
      // Templated class providing different evaluations..
      Evaluate();
  };

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

  void Evaluate() { output = F::Evaluate(input); }

  const typename F::Input &GetInput() const { return input; }

  const typename F::Output &GetOutput() const { return output; }

  friend std::ostream &operator<<(std::ostream &os,
                                  const TGenes<F> &ind) {
    auto indvector = ind.GetInput();
    for (int i = 0; i < indvector.size(); ++i)
      os << indvector[i] << " ";
    return os;
  }
};
}

#endif
