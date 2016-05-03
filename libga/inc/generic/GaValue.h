template <typename T> class GaValue {

protected:
  T fValue;

public:
  GaValue() {}
  GaValue(T value) : fValue(value) {}

  T GetValue() const { return fValue; }

  virtual void SetValue(const T &value) { GaValue::value = value; }
};

template <typename T>
bool operator<(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.getValue() < rhs.getValue();
}

template <typename T>
bool operator>(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.getValue() > rhs.getValue();
}

template <typename T>
bool operator<=(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.getValue() <= rhs.getValue();
}

template <typename T>
bool operator>=(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.getValue() >= rhs.getValue();
}

template <typename T>
bool operator==(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.getValue() == rhs.getValue();
}

template <typename T>
bool operator!=(const GaValue<T> &lhs, const GaValue<T> &rhs) {
  return lhs.getValue() != rhs.getValue();
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const GaValue<T> &rhs) {
  return os << rhs.getValue();
}
