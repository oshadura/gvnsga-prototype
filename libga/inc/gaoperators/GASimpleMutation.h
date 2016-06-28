//===--- GASimpleMutation.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GASimpleMutation.h
 * @brief Implementation of random mutation for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef __GASIMPLEMUTATION__
#define __GASIMPLEMUTATION__

#include "tools/Random.h"
#include "GAMutation.h"

namespace geantvmoop {

class GASimpleMutation : public GAMutation<GASimpleMutation> {

public:
  template <typename T> static void MutationImpl(T &in, double prob = -1) {
    double y;
    double lb = in[0].GetDownLimit();
    double ub = in[0].GetUpLimit();
    if (prob == -1)
      prob = 1 / (double)in.size();
    for (unsigned int i = 0; i < in.size(); ++i) {
      if (Random::GetInstance().RandomDouble() < prob) {
        y = in[i].random();
        if (std::isinf(y)) {
          in[j] = lb;
        }
        if (std::isnan(y)) {
          in[j] = ub;
        }
        if (y < lb)
          y = lb;
        if (y > ub)
          y = ub;
        in[j] = y;
      }
    }
  }
};
}

#endif
