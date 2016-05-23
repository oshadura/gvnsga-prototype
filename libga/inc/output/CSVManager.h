#ifndef __CSVMANAGER__
#define __CSVMANAGER__

#include <iostream>
#include <fstream>

#include "csv.h"
#include "generic/Population.h"
#include "generic/TGenes.h"

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
  ~CSVManager();

public:
  void CSVOutput(std::string file, const Population<F> &pop) {
    std::ofstream populationcvs;
    populationcvs.open(file.c_str(), std::fstream::app);
    // populationcvs << "Gene1, Gene2, Gene3\n";
    for (int j = 0; j < pop.size(); ++j) {
      // CVS format for R
      populationcvs << j << ",";
      for (int i = 0; i < 3; ++i) {
        populationcvs << pop[j][i] << ",";
      }
    }
    populationcvs
        << "-==========================================================-\n";
  }

  // Better to do it in variadic way..
  void LoadCSV(std::string file, const Population<F> &pop) {
    io::CSVReader<3> in(file);
    in.read_header(io::ignore_extra_column, "Gene1", "Gene2", "Gene3");
    double Gene1, Gene2, Gene3;
    while (in.read_row(Gene1, Gene2, Gene3)) {
      TGenes<F> individual;
      std::cout << "Reading " << pop.push_back() << std::endl;
    }
  }
};
}

#endif
