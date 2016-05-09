#ifndef __MPIEVALUATESERVER__
#define __MPIEVALUATESERVER__

#include <boost/mpi.hpp>
#include "generic/Population.h"
#include "generic/TGenes.h"

class ParallelEvaluateServer : private ParallelEvaluator,
                                  public EvaluatePopulation{

public:
  ParallelEvaluateServer(boost::mpi::environment &_mpi_env,
                            boost::mpi::communicator &_world,
                            ProblemDefinitions &_problem_defs)
      : ParallelEvaluator(_mpi_env, _world, _problem_defs) {
    // Send skeleton of decision variable to make sending dvs to clients/slaves
    // more efficient
    // Send skeleton of decision variable to make sending dvs to clients/slaves
    // more efficient
    decision_vars = std::pair<std::vector<double>, std::vector<int> >(
        std::piecewise_construct,
        std::forward_as_tuple(
            std::vector<double>(problem_defs.real_lowerbounds.size(), 0.0)),
        std::forward_as_tuple(
            std::vector<int>(problem_defs.int_lowerbounds.size(), 0)));

    objs_and_constraints = std::pair<std::vector<double>, std::vector<double> >(
        std::piecewise_construct,
        std::forward_as_tuple(
            std::vector<double>(problem_defs.minimise_or_maximise.size(), 0.0)),
        std::forward_as_tuple(
            std::vector<double>(problem_defs.number_constraints, 0.0)));

    boost::mpi::broadcast(world, boost::mpi::skeleton(decision_vars), 0);
    boost::mpi::broadcast(world, boost::mpi::skeleton(objs_and_constraints), 0);

    dv_c = boost::mpi::get_content(decision_vars);
    oc_c = boost::mpi::get_content(objs_and_constraints);
  }

  ~ParallelEvaluateServer() {
    // Send signal to slaves to indicate shutdown.
    for (int i = 0; i < number_clients; ++i) {
      int client_id = i + 1;
      world.send(client_id, max_tag, dv_c);
    }
  }

  void operator()(PopulationSPtr population) {
    // Sanity check - that we can represent each individual by an mpi tag.
    if (population->populationSize() > (max_tag - 1)) {
      //            std::cout << "problem: max tag too small, population too
      //            large for mpi\n";
    }

    int individual = 0;
    std::vector<boost::mpi::request> reqs_out(number_clients);
    for (; individual < number_clients; ++individual) {
      decision_vars.first = (*population)[individual].getRealDVVector();
      decision_vars.second = (*population)[individual].getIntDVVector();
      int client_id = individual + 1;
      //            std::cout << "sending to " << client_id << " individual " <<
      //            individual << " with " << decision_vars.first[0] << " " <<
      //            decision_vars.first[1] << std::endl;
      world.send(client_id, individual, dv_c);
    }
    //        mpi::wait_all(reqs_out.begin(), reqs_out.end());

    while (individual < population->populationSize()) {
      boost::mpi::status s =
          world.recv(boost::mpi::any_source, boost::mpi::any_tag, oc_c);
      //            std::cout << "received from " << s.source() << " individual
      //            " << individual << " with " << objs_and_constraints.first[0]
      //            << " " << objs_and_constraints.first[1] << std::endl;
      (*population)[s.tag()].setObjectives(objs_and_constraints.first);
      (*population)[s.tag()].setConstraints(objs_and_constraints.second);

      decision_vars.first = (*population)[individual].getRealDVVector();
      decision_vars.second = (*population)[individual].getIntDVVector();
      world.send(s.source(), individual, dv_c);

      ++individual;
    }

    for (int i = 0; i < number_clients; ++i) {
      boost::mpi::status s =
          world.recv(boost::mpi::any_source, boost::mpi::any_tag, oc_c);
      (*population)[s.tag()].setObjectives(objs_and_constraints.first);
      (*population)[s.tag()].setConstraints(objs_and_constraints.second);
    }
  }
};