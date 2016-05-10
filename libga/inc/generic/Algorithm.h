#ifndef __ALGORITHM__
#define __ALGORITHM__

#include <iostream>

#include "generic/PF.h"

template <typename Derived, typename F> class Algorithm {

private:
  F fFunction;

public:
  int fMaxGeneration = 100;

  Algorithm(F fFunction) : fFunction(fFunction) {}

  virtual PF<F> SolvePF(std::ostream *info = nullptr){
  	Initialize();
  	for (int i = 0; i < fMaxGeneration; ++i)
  	{
  		Evolution();
  		if (info != nullptr)
  		{
  			*info << "Generation" << i;
  			Print(*info);
  		}
  	}
  	return GetPF();
  }

  F GetProblem() const { return fFunction; }

  void SetProblem(F fFunction) { this->fFunction = fFunction; }

  void Evolution() { return static_cast<Derived *>(this)->Evolution(); }

  void Initialize() { return static_cast<Derived *>(this)->Initialize(); }

  void Print(std::ostream &os) {
    return static_cast<Derived *>(this)->Print(os);
  }

  PF<F> GetPF(){
  	return static_cast<Derived*>(this)->GetPF();
  }
};

#endif