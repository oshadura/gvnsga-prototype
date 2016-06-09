//===--- PolMutation.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PolMutation.h
 * @brief Implementation of  polynomial mutation for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef __GAPOLMUTATION__
#define __GAPOLMUTATION__

#include "tools/Random.h"
#include "GAMutation.h"

namespace geantvmoop {

class GAPolMutation : public GAMutation<GAPolMutation> {

public:
  template <typename T> static void MutationImpl(T &in, double prob = -1) {
    if (prob == -1)
      prob = 1 / (double)in.size();
    for (unsigned int i = 0; i < in.size(); ++i) {
      if (Random::GetInstance().RandomDouble() < prob)
        in[i] = in[i].random();
    }
  }
};
}

#endif
