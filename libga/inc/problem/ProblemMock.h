#ifndef __PROBLEMMOCK__
#define __PROBLEMMOCK__

#include "generic/Functions.h"
#include "tools/Generator.h"

#include <initializer_list>
#include <vector>

struct Mock : public std::vector<double> {
  Mock() {}
  Mock(std::initializer_list<double> list) : std::vector<double>(list) {}
  Mock(const std::vector<double> &v) {
    for (unsigned int i = 0; i < v.size(); ++i) {
      push_back(v[i]);
    }
  }

  Mock RandomSetup() const {
    std::vector<double> fIndividual;
    auto fGene = Generator::GetInstance();
    for (unsigned int i = 0; i < size(); ++i) {
      fIndividual.push_back(fGene.RNGDouble());
    }
    return Mock(fIndividual);
  }
};

class Problem : public Functions<Problem> {

public:
  typedef Mock Input;
  typedef std::vector<double> Output;
  static Output Evaluate(const Input &input) { return input; }
  static Input GetInput() { return Mock(std::vector<double>(2)); }
  static Output GetOutput() { return std::vector<double>(2); }
};

#endif 
