#ifndef __ALGORITHM__
#define __ALGORITHM__

#include <iostream>

template <typename Derived, typename F> class Algorithm {

private:
  F fFunction;

public:
  Algorithm(F fFunction) : fFunction(fFunction) {}

  F GetProblem() const { return fFunction; }

  void SetProblem(F fFunction) { this->fFunction = fFunction; }

  void Evolution() { return static_cast<Derived *>(this)->Evolution(); }

  void Initialize() { return static_cast<Derived *>(this)->Initialize(); }

  void Print(std::ostream &os) {
    return static_cast<Derived *>(this)->Print(os);
  }
};

#endif