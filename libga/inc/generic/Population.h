//===--- Population.h - LibGA -------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Population.h
 * @brief Implementation of population class for LibGA
 * prototype
 */
//
#pragma once

#ifndef __POPULATION__
#define __POPULATION__

#include "GADouble.h"
#include "GAVector.h"
#include "TGenes.h"
#include "instrumentation/CPUManager.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <list>
#include <sstream>
#include <stack>
#include <unordered_map>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>

#include <boost/archive/archive_exception.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/interprocess/anonymous_shared_memory.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/anonymous_shared_memory.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <boost/container/static_vector.hpp>

#include <iostream>
#include <cstring>

#include <boost/filesystem.hpp>

#include <boost/asio.hpp>
#include <boost/thread.hpp>

#define READ 0
#define WRITE 1
namespace geantvmoop {

template <typename F, std::size_t SizePop>
class Population
    : public boost::container::static_vector<individual_t<F>, SizePop> {
  // template <typename F> class Population : public std::vector<individual_t<F>
  // > {

public:
  Population(std::initializer_list<individual_t<F> > list)
      : boost::container::static_vector<individual_t<F>, SizePop>(list) {}

  Population() : boost::container::static_vector<individual_t<F>, SizePop>() {}

  Population(const boost::container::static_vector<individual_t<F>, SizePop> &
                 individuals)
      : boost::container::static_vector<individual_t<F>, SizePop>(individuals) {
  }

private:
  individual_t<F> ind;
  boost::container::static_vector<individual_t<F>, SizePop> population;

  /*
  // Cereal library (serialization)
  friend class cereal::access;

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &boost::serialization::base_object<TGenes<F>>(*this);
    ar &ind;
  }
  */

public:
  void InitSerilizationPopulation(
      boost::container::static_vector<individual_t<F>, SizePop> &population) {
    pid_t fArrayDead[SizePop];
    boost::filesystem::path endpoint_name =
        boost::filesystem::unique_path("/tmp/ga-%%%%-%%%%-%%%%-%%%%");
    boost::asio::io_service io_service_;
    boost::asio::local::stream_protocol::endpoint ep(endpoint_name.native());
    boost::asio::local::stream_protocol::acceptor acceptor(io_service_, ep);
    std::cout << "Endpoint was created..." << std::endl;
    for (int i = 0; i < SizePop; ++i) {
      CPUManager cpumgr;
      hwloc_topology_t topology;
      double nbcores, ccores;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      hwloc_topology_destroy(topology);
      ccores = nbcores -
               cpumgr.GetCurrentValueAllCPU() / 100 * nbcores; // just a test
      std::cout << "Number of total free cores " << ccores << std::endl;
      if (ccores < 1) {
        std::cout << "Sleeping, because free CPU ratio  " << ccores
                  << " is low.." << std::endl;
        sleep(1);

      } else {
        // std::stringstream ss;
        io_service_.notify_fork(boost::asio::io_service::fork_prepare);
        pid_t pid = fork();
        fArrayDead[i] = pid;
        if (pid == 0) {
          std::cout << "Generating individual in a child #" << i << std::endl;
          typename F::Input gene = F::GetInput().random();
          auto individual = std::make_shared<TGenes<F> >(gene);
          boost::asio::local::stream_protocol::iostream stream(ep);
          boost::archive::binary_oarchive oa(stream);
          oa &i &BOOST_SERIALIZATION_NVP(individual);
          ;
          auto indvector = (*individual).GetInput();
          for (int i = 0; i < indvector.size(); ++i)
            std::cout << indvector[i] << " ";
          std::cout << std::endl;
          exit(EXIT_SUCCESS);
        } else if (pid < 0) {
          std::cout << "Error on fork" << std::endl;
        } else {
          std::cout << "..." << std::endl;
        }
      }
    }
    for (int i = 0; i < SizePop; ++i) {
      // int forki = 0;
      individual_t<F> transfer;
      std::cout << "=Creating a listener #" << i << " =" << std::endl;
      {
        boost::asio::local::stream_protocol::stream_protocol::iostream stream;
        boost::system::error_code ec;
        acceptor.accept(*stream.rdbuf(), ec);
        boost::archive::binary_iarchive ia(stream);
        ia &i &BOOST_SERIALIZATION_NVP(transfer);
      }
      std::cout << "Pushing a new gene in a parent..." << std::endl;
      this->push_back(transfer);
      auto indvector = (*transfer).GetInput();
      for (int i = 0; i < indvector.size(); ++i)
        std::cout << indvector[i] << " ";
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.."
                << std::endl;
      waitpid(fArrayDead[i], NULL, 0);
      std::cout << "PID: " << fArrayDead[i] << " has shut down.." << std::endl;
      io_service_.notify_fork(boost::asio::io_service::fork_parent);
    }
    boost::filesystem::remove(endpoint_name);
    std::fill(fArrayDead, fArrayDead + SizePop, 0);
  }

  void InitPipePopulation(
      boost::container::static_vector<individual_t<F>, SizePop> &population) {
    int pipega[SizePop][2];
    pid_t cpid;
    ssize_t result;
    pipe(pipega[SizePop]);
    double readervalue;
    pid_t fArrayDead[SizePop];
    for (int i = 0; i < SizePop; ++i) {
      CPUManager cpumgr;
      hwloc_topology_t topology;
      double nbcores, ccores;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      hwloc_topology_destroy(topology);
      ccores = nbcores -
               cpumgr.GetCurrentValueAllCPU() / 100 * nbcores; // just a test
      std::cout << "Number of total free cores " << ccores << std::endl;
      if (ccores < 1) {
        std::cout << "Sleeping, because free CPU ratio  " << ccores
                  << " is low.." << std::endl;
        sleep(10);
      } else {
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        pid_t pid = fork();
        fArrayDead[i] = pid;
        ////////////////////////////=CHILD=//////////////////////
        if (pid == 0) {
          typename F::Input gene = F::GetInput().random();
          auto individual = std::make_shared<TGenes<F> >(gene);
          auto indvector = (*individual).GetInput();
          auto fit = (*individual).GetOutput();
          // for (int i = 0; i < indvector.size(); ++i)
          //  std::cout << indvector[i] << " ";
          // std::cout << std::endl;
          close(pipega[SizePop][READ]);
          for (auto it : indvector) {
            // lockf(pipega[SizePop][WRITE], F_LOCK, 0);
            write(pipega[SizePop][WRITE], &it, sizeof(double));
            // lockf(pipega[SizePop][WRITE], F_ULOCK, 0);
            std::cout << "Individual part to be send: " << it << std::endl;
          }
          for (auto iter : fit) {
            // lockf(pipega[SizePop][WRITE], F_LOCK, 0);
            write(pipega[SizePop][WRITE], &iter, sizeof(double));
            // lockf(pipega[SizePop][WRITE], F_ULOCK, 0);
            std::cout << "Fitness part to be send: " << iter << std::endl;
          }
          close(pipega[WRITE]); // close the read-end of the pipe
          wait(NULL);
          exit(EXIT_SUCCESS);
        } else if (pid < 0) {
          ////////////////////////////=NOPE=//////////////////////
          exit(EXIT_FAILURE);
          std::cout << "Error on fork" << std::endl;
        } else {
          ////////////////////////////=PARENT=////////////////////
          /*
          typename F::Input tmpinput;
          typename F::Output tmpoutput;
          fArrayDead[i] = cpid;
          close(pipega[WRITE]);
          std::cout << "We are starting to read on master job.." << std::endl;
          memset(&tmpoutput, 0, sizeof(tmpoutput));
          memset(&tmpinput, 0, sizeof(tmpinput));
          while (read(pipega[READ], &readervalue, sizeof(double)) > 0) {
            std::cout << "I am in a reader loop for " <<SizePop<< std::endl;
            for (int i = 0; i < tmpinput.size(); ++i) {
              std::cout << "Individual part to be received: " << readervalue
                        << std::endl;
              tmpinput.push_back(readervalue);
            }
            for (int i = 0; i < tmpoutput.size(); ++i) {
              std::cout << "Fitness part to be received: " << readervalue
                        << std::endl;
              tmpoutput.push_back(readervalue);
            }
          }
          auto forkedindividual =
              std::make_shared<TGenes<F> >(tmpinput, tmpoutput);
          auto indv = (*forkedindividual).GetInput();
          std::cout << "--------------------------------------" << std::endl;

          for (int i = 0; i < indv.size(); ++i)
            std::cout << indv[i] << " ";
          std::cout << "--------------------------------------" << std::endl;
          this->push_back(forkedindividual);
          // std::cout << "We are stoping to read.." << std::endl;
          close(pipega[READ]);
        */
        }
      }

      for (int i = 0; i < SizePop; ++i) {
        typename F::Input tmpinput;
        typename F::Output tmpoutput;
        // fArrayDead[i] = cpid;
        close(pipega[SizePop][WRITE]);
        std::cout << "We are starting to read on master job.." << std::endl;
        memset(&tmpoutput, 0, sizeof(tmpoutput));
        memset(&tmpinput, 0, sizeof(tmpinput));
        while (read(pipega[SizePop][READ], &readervalue, sizeof(double)) > 0) {
          std::cout << "I am in a reader loop for " << SizePop << std::endl;
          for (int i = 0; i < tmpinput.size(); ++i) {
            std::cout << "Individual part to be received: " << readervalue
                      << std::endl;
            tmpinput.push_back(readervalue);
          }
          for (int i = 0; i < tmpoutput.size(); ++i) {
            std::cout << "Fitness part to be received: " << readervalue
                      << std::endl;
            tmpoutput.push_back(readervalue);
          }
        }
        auto forkedindividual =
            std::make_shared<TGenes<F> >(tmpinput, tmpoutput);
        auto indv = (*forkedindividual).GetInput();
        std::cout << "--------------------------------------" << std::endl;

        for (int i = 0; i < indv.size(); ++i)
          std::cout << indv[i] << " ";
        std::cout << "--------------------------------------" << std::endl;
        this->push_back(forkedindividual);
        // std::cout << "We are stoping to read.." << std::endl;
        close(pipega[SizePop][READ]);
        std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.."
                  << std::endl;
        waitpid(fArrayDead[i], NULL, 0);
        std::cout << "PID: " << fArrayDead[i] << " has shut down.."
                  << std::endl;
      }
      std::fill(fArrayDead, fArrayDead + SizePop, 0);
    }
  }

  void InitSharedMemPopulation(
      boost::container::static_vector<individual_t<F>, SizePop> &population) {
    int pipega[SizePop][2];
    pid_t cpid;
    ssize_t result;
    pipe(pipega[SizePop]);
    double readervalue;
    pid_t fArrayDead[SizePop];
    using namespace boost::interprocess;
    for (int i = 0; i < SizePop; ++i) {
      CPUManager cpumgr;
      hwloc_topology_t topology;
      double nbcores, ccores;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      hwloc_topology_destroy(topology);
      ccores = nbcores -
               cpumgr.GetCurrentValueAllCPU() / 100 * nbcores; // just a test
      std::cout << "Number of total free cores " << ccores << std::endl;
      if (ccores < 1) {
        std::cout << "Sleeping, because free CPU ratio  " << ccores
                  << " is low.." << std::endl;
        sleep(10);
      } else {
        //////////////////////////////////////////////////////////
        pid_t pid = fork();
        fArrayDead[i] = pid;
        if (pid == 0) {
          typename F::Input gene = F::GetInput().random();
          auto individual = std::make_shared<TGenes<F> >(gene);
          auto indv = (*individual).GetInput();
          std::cout << "--------------------------------------" << std::endl;
          for (int i = 0; i < indv.size(); ++i)
            std::cout << indv[i] << " ";
          std::cout << "--------------------------------------" << std::endl;
          population.push_back(individual);
          wait(NULL);
          exit(EXIT_SUCCESS);
        } else if (pid < 0) {
          exit(EXIT_FAILURE);
          std::cout << "Error on fork" << std::endl;
        } else {
          sleep(50);
        }
      }
      for (int i = 0; i < SizePop; ++i) {
        std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.."
                  << std::endl;
        waitpid(fArrayDead[i], NULL, 0);
        std::cout << "PID: " << fArrayDead[i] << " has shut down.."
                  << std::endl;
      }
      std::fill(fArrayDead, fArrayDead + SizePop, 0);
    }
    // std::cout << population;
  }

  void InitSimplePopulation(
      boost::container::static_vector<individual_t<F>, SizePop> &population) {
    for (int i = 0; i < SizePop; ++i) {
      CPUManager cpumgr;
      cpumgr.InitCPU();
      hwloc_topology_t topology;
      double nbcores, ccores;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      hwloc_topology_destroy(topology);
      ccores =
          nbcores - cpumgr.GetCurrentValueCPU() / 100 * nbcores; // just a test
      std::cout << " Number of total free cores " << ccores << std::endl;
      if (ccores < 0.3) {
        std::cout << "Sleeping, because free CPU ratio  " << ccores
                  << " is low.." << sleep(50);
      } else {
        typename F::Input gene = F::GetInput().random();
        auto individual = std::make_shared<TGenes<F> >(gene);
        population.push_back(individual);
      }
      // std::cout << &&population;
    }
  }

  ~Population() {}

  /*
  #if defined __clang__
   void push_back(individual_t<F> ind) const { (*this).push_back(ind); }
  #endif
  */

  template <class Archive> void save(Archive &ar) const {
    // ar(ind);
    ar << ind;
  }

  template <class Archive> void load(Archive &ar) {
    // ar(ind);
    ar >> ind;
  }

  const typename F::Input &GetTGenes(int i) const {
    return (*this)[i]->GetInput();
  }

  const typename F::Output &GetTFitness(int i) const {
    return (*this)[i]->GetOutput();
  }

  typename F::Input GetGene(int gene) const {
    typename F::Input v;
    for (unsigned int j = 0; j < this->size(); ++j)
      v.push_back(GetGeneValue(j, gene));
    return v;
  }

  double GetGeneValue(int index, int gene) const {
    auto value = (*this)[index]->GetInput()[gene];
    return value.GetGAValue();
  }

  typename F::Output GetObjective(int objective) const {
    typename F::Output v;
    for (unsigned int j = 0; j < this->size(); ++j)
      v.push_back(GetObjectiveValue(j, objective));
    return v;
  }

  double GetObjectiveValue(int index, int objective) const {
    return (*this)[index]->GetOutput()[objective];
  }

  bool IsNonDominated(const individual_t<F> &ind) const {
    for (auto entry : *this) {
      if (ind->IsDominatedBy(*entry))
        return true;
    }
    return false;
  }

  const individual_t<F> PopBack() {
    individual_t<F> a = this->back();
    this->pop_back();
    return a;
  }

  void Remove(const Population<F, SizePop> &pop) {
    for (auto entry : pop) {
      this->erase(std::remove(this->begin(), this->end(), entry), this->end());
    }
  }

  template <typename T, std::size_t SizeGene>
  void SortVector(const boost::container::static_vector<T, SizeGene> &v, bool isDescending = false) {
    std::unordered_map<individual_t<F>, T> m;
    for (int i = 0; i < this->size(); ++i)
      m[(*this)[i]] = v[i];
    SortMap(m, isDescending);
  }

  template <typename T, std::size_t SizeGene>
  std::vector<int> SortIndex(const boost::container::static_vector<T, SizeGene> &v,
                             bool isDescending = false) const {
    std::vector<int> index = GetIndex();
    if (isDescending) {
      std::sort(
          index.begin(), index.end(),
          [&v](const int &lhs, const int &rhs) { return v[lhs] > v[rhs]; });
    } else
      std::sort(
          index.begin(), index.end(),
          [&v](const int &lhs, const int &rhs) { return v[lhs] < v[rhs]; });
    return index;
  }

  template <typename T>
  void SortMap(std::unordered_map<individual_t<F>, T> &m,
               bool isDescending = false) {
    if (isDescending) {
      std::sort(this->begin(), this->end(),
                [&m](const individual_t<F> &lhs,
                     const individual_t<F> &rhs) { return m[lhs] > m[rhs]; });
    } else
      std::sort(this->begin(), this->end(),
                [&m](const individual_t<F> &lhs,
                     const individual_t<F> &rhs) { return m[lhs] < m[rhs]; });
  }

  void SortObj(int objective, bool isDescending = false) {
    SortVec(GetObjective(objective), isDescending);
  }

  std::vector<int> GetIndex() const {
    std::vector<int> index(this->size());
    for (unsigned int k = 0; k < this->size(); ++k)
      index[k] = k;
    return index;
  }

  friend std::ostream &operator<<(std::ostream &s,
                                  const Population<F, SizePop> &pop) {
    std::cout << "---------------------------\n";
    std::cout << "Size of population: " << pop.size() << std::endl;
    std::cout << "---------------------------\n" << std::endl;
    for (int i = 0; i < pop.size(); ++i) {
      std::cout << "Individual " << i << std::endl;
      for (int j = 0; j < pop.GetTGenes(i).size(); ++j) {
        std::cout << pop.GetGeneValue(i, j) << "|";
      }
      std::cout << "\nFitness function value: " << std::endl;
      for (int k = 0; k < pop.GetTFitness(i).size(); ++k) {
        std::cout << pop.GetObjectiveValue(i, k) << "|";
      }
    }
    std::cout << "---------------------------\n" << std::endl;
    return s;
  }
};
}

#endif
