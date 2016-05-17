#include "Random.h"
#include <cstdlib>

moo::Random *moo::Random::_singletonInst = new moo::Random();

double moo::Random::rndDouble() { return ((double)std::rand() / (RAND_MAX)); }

int moo::Random::rndInt(int min, int max) {
  return min + (std::rand() % (int)(max - min + 1));
}

bool moo::Random::rndBool() { return rand() % 2 == 1; }

double moo::Random::rndDouble(int min, int max) {
  return min + (double)std::rand() / RAND_MAX * (max - min);
}
