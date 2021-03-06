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
#pragma once

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
  void CSVOutput(std::string file, const Population<F> &population, std::unordered_map<individual_t<F>, int> fIndRank, std::unordered_map<individual_t<F>, double> fIndCrowDist) {
    std::ofstream populationcvs;
    populationcvs.open(file.c_str(), std::fstream::app);
    // Suppose to be variadic...
    populationcvs << "Gene1, Gene2, Gene3, Gene4, Gene5, Gene6, Gene7, Gene8"
                     "Fitness1, Fitness2, Fitness3, Rank, CrowdingDistance\n";
    // CVS format for R ananlysis..
    for (std::size_t i = 0; i < population.size(); ++i) {
      for (std::size_t j = 0; j < population.GetTGenes(i).size(); ++j) {
        populationcvs << population.GetGeneValue(i, j) << ",";
      }
      for (std::size_t k = 0; k < population.GetTFitness(i).size(); ++k) {
        populationcvs << population.GetObjectiveValue(i, k) << ",";
      }
      auto ind = population[i];
      populationcvs << fIndRank[ind] << ",";
      populationcvs << fIndCrowDist[ind];
      populationcvs << "\n";
    }
    populationcvs << "-=================-\n";
  }

  // Better to do it in variadic way..
  /**
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
      // TBD fixed for MAC OS X
      pop.push_back(fInd);
    }
  }
  **/
};
}

#endif
