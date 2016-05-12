#ifndef __FUNCTIONS__
#define __FUNCTIONS__

namespace geantvmoop {

template <typename DerivedClass> class Functions {

public:
  int fGen = 0;
  static int GetNObjectives() { return DerivedClass::GetOutput().size(); }
};
}

#endif
