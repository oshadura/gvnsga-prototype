template <typename T> class GaVector : public std::vector<T> {

#include <vector>

public:
  GaVector() : std::vector<T>() {}
  GaVector(int n, const T &v) : std::vector<T>(n, v) {}
  //Constructor based on Limits and Generators ?
};