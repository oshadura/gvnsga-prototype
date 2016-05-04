#ifndef __CSVMANAGER__
#define __CSVMANAGER__


#include "csv.h"
#include "generic/Population.h"
#include "generic/TGenes.h"

template <typename F> class CSVManager {
public:

  static CSVManager& GetInstance(){
  	static CSVManager Instance;
  	return Instance;
  }

  CSVManager(CSVManager const&) = delete;
  void operator=(CSVManager const&) = delete;

private:

  CSVManager(){}
  ~CSVManager();

  // Old functions to be moved in a cxx
	/* 
  std::ofstream &Population<F>::CreateCVS(std::string file) {
    std::ofstream populationcvs;
    populationcvs.open(file.c_str(), std::fstream::app);
  }

  void CSVOutput(std::ofstream &populationcvs, const Population<F> &pop) {
    populationcvs
        << "Fitness1, Fitness2, Fitness3, Gene1, Gene2, Gene3, Gene4, "
           "Gene5, Gene6, Gene7, ConstrainViolation,Crowdingdistance, "
           "Rank\n";
    for (int j = 0; j < const_cast<Population<F> &>(pop).GetPopulationSize();
         ++j) {
      // CVS format for R
      populationcvs << j << ",";
      if (setupPop.fNObjectives > 0) {
        for (int i = 0; i < setupPop.fNObjectives; ++i) {
          populationcvs << const_cast<Population<F> &>(pop)
                               .GetGenes(j)
                               .GetFitness(i) << ",";
        }
      }
      if (setupPop.fNCons > 0) {
        for (int i = 0; i < setupPop.fNCons; ++i) {
          populationcvs << const_cast<Population<F> &>(pop)
                               .GetGenes(j)
                               .GetConstrain(i) << ",";
        }
      }
      if (setupPop.fNParam > 0) {
        for (int i = 0; i < setupPop.fNParam; ++i) {
          populationcvs << const_cast<Population<F> &>(pop).GetGenes(j).GetGene(
                               i) << ",";
        }
      }
      populationcvs << const_cast<Population<F> &>(pop)
                           .GetGenes(j)
                           .GetConsViol() << ",";
      populationcvs << const_cast<Population<F> &>(pop)
                           .GetGenes(j)
                           .GetCrowdingDistance() << ",";
      populationcvs << const_cast<Population<F> &>(pop).GetGenes(j).GetRank()
                    << ","
                    << "\n";
    }
    populationcvs
        << "-=========================================================="
           "===================-\n";
  }

  // Better to do it in variadic way..
  void LoadCSV(const Population<F> &pop) {
    io::CSVReader<13> in("PopulationDTLZ1.cvs");
    in.read_header(io::ignore_extra_column, "Fitness1", "Fitness2", "Fitness3",
                   "Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6",
                   "Gene7", "ConstrainViolation", "Crowdingdistance", "Rank");
    Double_t Fitness1, Fitness2, Fitness3, Gene1, Gene2, Gene3, Gene4, Gene5,
        Gene6, Gene7, ConstrainViolation, Crowdingdistance, Rank;
    while (in.read_row(Fitness1, Fitness2, Fitness3, Gene1, Gene2, Gene3, Gene4,
                       Gene5, Gene6, Gene7, ConstrainViolation,
                       Crowdingdistance, Rank)) {
      std::cout << "Reading " << Fitness1 << std::endl;
    }
  }

*/
};

#endif