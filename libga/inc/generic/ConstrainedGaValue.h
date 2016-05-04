#ifndef __CONSTRAINEDGAVALUE__
#define __CONSTRAINEDGAVALUE__


#include "generic/GaValue.h"

template <typename T> class ConstrainedGaValue : public GaValue<T> {

protected:
  T fDownValue;
  T fUpValue;

public:
  ConstrainedGaValue(T v) {}
  ConstrainedGaValue(T down, T up) : fDownValue(down), fUpValue(up) {}
  ConstrainedGaValue(T v, T down, T up)
      : GaValue<T>(v), fDownValue(down), fUpValue(up) {}

private:
  void SetDownValue(T down) { fDownValue = down; }
  void SetUpValue(T up) { fUpValue = up; }
  T GetDownValue() const { return fDownValue; }
  T GetUpValue() const { return fUpValue; }
};

#endif