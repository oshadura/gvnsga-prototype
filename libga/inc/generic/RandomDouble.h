#ifndef __RANDOMDOUBLE__
#define __RANDOMDOUBLE__

#include "generic/ConstrainedGaValue.h"

namespace geantvmoop{

class RandomDouble : public ConstrainedGaValue<double> {

public:
  RandomDouble(double v) : ConstrainedGaValue(v, 0, 1) {};
  RandomDouble(double v, double down, double up)
      : ConstrainedGaValue(v, down, up) {};
  RandomDouble(double down, double up) : ConstrainedGaValue(0, down, up) {};
  RandomDouble RandomSetup(){
  	auto generator = Generator::GetInstance();
  	double value = generator.RNGDouble(fDownValue, fUpValue);
  	return RandomDouble(value, fDownValue, fUpValue);
  }
};

}

#endif