/*

#include "MPIServer.h"
#include "Population.h"
void MPIServer::Launch()(Population<Double_t> pop){
    std::vector<mpi::request> reqs_out(pop->GetPopulationSize());
    std::vector<mpi::request> reqs_in(pop->GetPopulationSize());
    for (int i = 0; i < pop->GetPopulationSize(); ++i){
        reqs_out[i] = boost::mpi::isend();
    }
    for (int i = 0; i < pop->GetPopulationSize(); ++i){
        //std::vector<double> objectives;
        //std::vector<double> constraints;
        reqs_in[i] = boost::mpi::irecv();
        //ind.setObjectives(objectives);
        //ind.setConstraints(constraints);
        }
    }
*/
