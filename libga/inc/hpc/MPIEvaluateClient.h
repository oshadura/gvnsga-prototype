#ifndef MPIEvaluateClient_h
#define MPIEvaluateClient_h

#include <boost/mpi.hpp>
#include "generic/Population.h"
#include "generic/TGenes.h"

class ParallelEvaluateClient : private ParallelEvaluator {
  EvaluatorBase &eval;

public:
  ParallelEvaluateClient(boost::mpi::environment &_mpi_env,
                            boost::mpi::communicator &_world,
                            ProblemDefinitions &_problem_defs,
                            EvaluatorBase &_eval)
      : ParallelEvaluatorBase(_mpi_env, _world, _problem_defs), eval(_eval) {
    // Send skeleton of decision variable to make sending dvs to clients/slaves
    // more efficient
    boost::mpi::broadcast(world, boost::mpi::skeleton(decision_vars), 0);
    boost::mpi::broadcast(world, boost::mpi::skeleton(objs_and_constraints), 0);

    dv_c = boost::mpi::get_content(decision_vars);
    oc_c = boost::mpi::get_content(objs_and_constraints);
  }

  void operator()() {

    bool do_continue = true;
    while (do_continue) {
      //            std::cout << "waiting to receive" << std::endl;
      boost::mpi::status s = world.recv(0, boost::mpi::any_tag, dv_c);
      //            std::cout << " received " << decision_vars.first[0] << " "
      //            << decision_vars.first[1] << " for individual " << s.tag()
      //            << std::endl;
      if (s.tag() == max_tag) {
        do_continue = false;
      } else {
        // calc objective
        objs_and_constraints = eval(decision_vars.first, decision_vars.second);
        world.send(0, s.tag(), oc_c);
      }
    }
  }
};