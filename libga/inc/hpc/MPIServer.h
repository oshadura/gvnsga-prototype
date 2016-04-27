#ifndef MPIServer_h
#define MPIServer_h

#include <boost/mpi.hpp>
#include "generic/TGenes.h"

class MPIServer {
public:
  MPIServer(boost::mpi::environment &me, boost::mpi::communicator &w)
      : mpienv(me), world(w) {
    boost::mpi::broadcast(world, boost::mpi::skeleton(Genes<Double_t> & ind));
  }
  ~MPIServer();
  void Launch();

private:
  boost::mpi::environment &mpienv;
  boost::mpi::communicator &world;
};

#endif
