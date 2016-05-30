//===--- CSVManager.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file CSVManager.h
 * @brief Implementation of CSV class for LibGA
 * prototype
 */
//

#ifndef __CSVMANAGER__
#define __CSVMANAGER__

#include <iostream>
#include <fstream>

#include "csv.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"

namespace geantvmoop {

class CSVManager {

public:
  static CSVManager &GetInstance() {
    static CSVManager Instance;
    return Instance;
  }

  CSVManager(CSVManager const &) = delete;
  void operator=(CSVManager const &) = delete;

private:
  CSVManager() {}
  ~CSVManager() {}

public:
  template <typename F>
  void CSVOutput(std::string file, const Population<F> &pop) {
    std::ofstream populationcvs;
    populationcvs.open(file.c_str(), std::fstream::app);
    // Suppose to be variadic...
    populationcvs << "Gene1, Gene2, Gene3\n";
    // CVS format for R ananlysis..
    for (int i = 0; i < pop.size(); ++i) {
      auto individual = pop.GetTGenes(i);
      for (int i = 0; i < individual.size(); ++i) {
        auto parameter = individual[i];
        populationcvs << parameter << ",";
      }
      populationcvs << "\n";
    }
    populationcvs << "-=================-\n";
  }

  // Better to do it in variadic way..
  template <typename F>
  void LoadCSV(std::string file, const Population<F> &pop) {
    double Gene1, Gene2, Gene3;
    io::CSVReader<3> in(file);
    in.read_header(io::ignore_extra_column, "Gene1", "Gene2", "Gene3");
    while (in.read_row(Gene1, Gene2, Gene3)) {
      GAVector<GADouble> individual;
      individual.push_back(Gene1);
      individual.push_back(Gene2);
      individual.push_back(Gene3);
      auto fInd = std::make_shared<TGenes<F>>(individual);
      // TBF fixed for MAC OS X
      pop.push_back(fInd);
    }
  }
};
}

#endif
