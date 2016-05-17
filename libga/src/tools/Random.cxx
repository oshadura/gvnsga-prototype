#include "tools/Random.h"
#include <cstdlib>

geantvmoop::Random *geantvmoop::Random::_singletonInst = new geantvmoop::Random();

double geantvmoop::Random::rndDouble() { return ((double)std::rand() / (RAND_MAX)); }

int geantvmoop::Random::rndInt(int min, int max) {
  return min + (std::rand() % (int)(max - min + 1));
}

bool geantvmoop::Random::rndBool() { return rand() % 2 == 1; }

double geantvmoop::Random::rndDouble(int min, int max) {
  return min + (double)std::rand() / RAND_MAX * (max - min);
}
