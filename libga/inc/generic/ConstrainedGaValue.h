#include "GaValue.h"

template <typename T> class ConstrainedGAValue : public GaValue<T> {

protected:
  T fDownValue;
  T fUpValue;

public:
  ConstrainedGAValue() {}
  ConstrainedGAValue(T down, T up) : fDownValue(down), fUpValue(up) {}
  ConstrainedGAValue(T v, T down, T up)
      : GaValue<T>(v), fDownValue(down), fUpValue(up) {}

private:
  void SetDownValue(T down) { fDownValue = down; }
  void SetUpValue(T up) { fUpValue = up; }
  T GetDownValue() const { return fDownValue; }
  T GetUpValue() const { return fUpValue; }
};