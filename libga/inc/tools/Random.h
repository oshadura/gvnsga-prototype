#ifndef MOO_RANDOM_H
#define MOO_RANDOM_H

#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */

namespace geantvmoop {

class Random {
public:
  static Random *getInstance() { return _singletonInst; };

  double rndDouble();
  double rndDouble(int min, int max);
  int rndInt(int min, int max);
  bool rndBool();

private:
  static Random *_singletonInst;
  Random() { srand(time(NULL)); }
  Random(Random const &);
  void operator=(Random const &);
};
}

#endif // MOO_RANDOM_H
