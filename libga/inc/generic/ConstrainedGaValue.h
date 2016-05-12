#ifndef __CONSTRAINEDGAVALUE__
#define __CONSTRAINEDGAVALUE__

#include "generic/GaValue.h"

namespace geantvmoop{

template <typename T> class ConstrainedGaValue : public GaValue<T> {

protected:
  T fDownValue;
  T fUpValue;

public:
  ConstrainedGaValue() {}
  ConstrainedGaValue(T down, T up) : fDownValue(down), fUpValue(up) {}
  ConstrainedGaValue(T v, T down, T up)
      : GaValue<T>(v), fDownValue(down), fUpValue(up) {}
  void SetDownBound(T down) { fDownValue = down; }
  void SetUpBound(T up) { fUpValue = up; }
  T GetDownBound() const { return fDownValue; }
  T GetUpBound() const { return fUpValue; }
};

}

#endif