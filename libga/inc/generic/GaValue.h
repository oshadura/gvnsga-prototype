#ifndef __GAVALUE__
#define __GAVALUE__

template <typename T> class GaValue {

protected:
  T fValue;

public:
  GaValue() {}
  GaValue(T fValue) : fValue(fValue) {}

  T GetValue() const { return fValue; }

  virtual void SetValue(const T &value) { GaValue::fValue = value; }
};

template <typename T>
bool operator<(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.GetValue() < rhs.GetValue();
}

template <typename T>
bool operator>(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.GetValue() > rhs.GetValue();
}

template <typename T>
bool operator<=(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.GetValue() <= rhs.GetValue();
}

template <typename T>
bool operator>=(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.GetValue() >= rhs.GetValue();
}

template <typename T>
bool operator==(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.GetValue() == rhs.GetValue();
}

template <typename T>
bool operator!=(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.GetValue() != rhs.GetValue();
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const GaValue<T> &rhs) {
  return os << rhs.GetValue();
}

#endif
