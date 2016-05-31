#ifndef MOO_RANDOM_H
#define MOO_RANDOM_H

#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <cstdlib>

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
    return min + (std::rand() % (int)(max - min + 1));
  }

  bool RandomBool() { return rand() % 2 == 1; }

  double RandomDouble(int min, int max) {
    return min + (double)std::rand() / RAND_MAX * (max - min);
  }

private:
  Random() { srand(time(NULL)); }
};
}

#endif
