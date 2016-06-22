#pragma once

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

namespace geantvmoop {

class GAMPIEvaluator : public GAEvaluate<GAMPIEvaluator> {

public:
  // void Evaluate() { output = F::Evaluate(input); }
  template <typename T> static void Evaluate() {
    boost::mpi::environment env;
    boost::mpi::communicator world;
    Int eval; // Summary
    if (world.rank() == 0) {
      // Create evaluator server
      ParallelEvaluatePopServer eval_server(env, world,
                                            eval.getProblemDefinitions());
      // Send population to evaluator server..
      // Population<T> pop;
      eval_server(pop);
      // Boost for loop..
      BOOST_FOREACH(Individual & ind, *pop) {
        std::cout << ind.getRealDV(0) << "\t" << ind.getRealDV(1) << "\t"
                  << ind.getRealDV(2) << "\t" << ind.getObjective(0) << "\t"
                  << ind.getObjective(1) << std::endl;
      }

    } else {
      // Create evaluator client
      ParallelEvaluatePopClient eval_client(env, world,
                                            eval.getProblemDefinitions(), eval);
      eval_client();
    }
  }
};
}