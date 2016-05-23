#ifndef __CSVMANAGER__
#define __CSVMANAGER__

#include <iostream>
#include <fstream>

#include "csv.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/GAVector.h"

namespace geantvmoop {

template <typename F> class CSVManager {
  friend class Population<F>;

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
  void CSVOutput(std::string file, const Population<F> &pop) {
    std::ofstream populationcvs;
    populationcvs.open(file.c_str(), std::fstream::app);
    // populationcvs << "Gene1, Gene2, Gene3\n";
    /*
    for (int j = 0; j < pop.size(); ++j) {
      // CVS format for R
      populationcvs << j << ",";
      auto individual = pop[j];
      for (auto i : individual) {
        populationcvs << individual << ",";
      }
    }
    */
    populationcvs
        << "-==========================================================-\n";
  }

  // Better to do it in variadic way..
  void LoadCSV(std::string file, const Population<F> &pop) {
    double Gene1, Gene2, Gene3;
    io::CSVReader<3> in(file);
    in.read_header(io::ignore_extra_column, "Gene1", "Gene2", "Gene3");
    while (in.read_row(Gene1, Gene2, Gene3)) {
      GAVector<GADouble> individual;
      // GAVector<GADouble> fInd = static_cast<GAVector<GADouble>>(individual);
      individual.push_back(Gene1);
      /*
      std::cout << "Reading Gene1: " << fInd.push_back(Gene1) << std::endl;
      std::cout << "Reading Gene2: " << individual.push_back(Gene2)
                << std::endl;
      std::cout << "Reading Gene3: " << individual.push_back(Gene3)
                << std::endl;
      */
      auto fInd = std::make_shared<geantvmoop::TGenes<F>>(individual);
      pop.push_back(fInd);
    }
  }
};
}

#endif
