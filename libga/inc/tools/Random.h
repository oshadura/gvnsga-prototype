#pragma once

#ifndef MOO_RANDOM_H
#define MOO_RANDOM_H

#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <cstdlib>
#include <ctime>

namespace geantvmoop {

class Random {
public:
  static Random &GetInstance() {
    static Random Instance;
    return Instance;
  }

  Random(Random const &);
  void operator=(Random const &);

  double RandomDouble() { return ((double)std::rand() / (RAND_MAX)); }

  int RandomInt(int min, int max) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> unif(min, max);
    double a_random_int = unif(mt);
    return a_random_int;
  }

  bool RandomBool() { return rand() % 2 == 1; }

  double RandomDouble(int min, int max) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> unif(min, max);
    std::default_random_engine re;
    double a_random_double = unif(mt);
    return a_random_double;
  }

private:
  Random() { srand(time(NULL)); }
};
}

#endif
