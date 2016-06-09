#pragma once

#ifndef __MPIEVALUETE__
#define __MPIEVALUETE__

#include <boost/mpi.hpp>
#include "generic/Population.h"
#include "generic/TGenes.h"

namespace geantvmoop{

class MPIEvaluate {
protected:
  boost::mpi::environment &mpi_env;
  boost::mpi::communicator &world;
  Functions &func;
  int fNProcesses;
  int fNClients;
  // std::pair<std::vector<double>, std::vector<int> > decision_vars;
  // std::pair<std::vector<double>, std::vector<double> > objs_and_constraints;
  boost::mpi::content dv_c;
  boost::mpi::content oc_c;
  int max_tag;

public:
  MPIEvaluate(boost::mpi::environment &_mpi_env,
              boost::mpi::communicator &_world,
              ProblemDefinitions &_problem_defs)
      : mpi_env(_mpi_env), world(_world), problem_defs(_problem_defs),
        number_processes(world.size()), number_clients(number_processes - 1),
        max_tag(mpi_env.max_tag()) {}
};

}

#endif 
