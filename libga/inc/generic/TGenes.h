//===--- TGenes.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TGenes.h
 * @brief Implementation of TGenes for LibGA
 * prototype
 */

#pragma once

#ifndef __TGENES__
#define __TGENES__

#include <functional>
#include <memory>
#include <vector>

#include <hwloc.h>
#include <unistd.h>

#include "hpc/GAEvaluate.h"
#include "hpc/GASequentualEvaluator.h"
#include "hpc/GASimpleEvaluator.h"

#include "instrumentation/CPUManager.h"

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>

#define READ 0
#define WRITE 1

namespace geantvmoop {

template <typename F> class TGenes;

template <typename F> using individual_t = std::shared_ptr<TGenes<F> >;

template <typename F> class TGenes {

private:
  typename F::Input input;
  typename F::Output output;
  // Temporary object..
  typename F::Output tmpoutput;

public:
  TGenes() {}

  TGenes(const typename F::Input &i, bool fEvaluated = true) : input(i) {
    if (fEvaluated) {
      CPUManager cpumgr;
      cpumgr.InitCPU();
      hwloc_topology_t topology;
      double nbcores, ccores;
      unsigned int microseconds;
      microseconds = 3000;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      // printf(" Number of free cores %d cores\n", nbcores);
      hwloc_topology_destroy(topology);
      ccores =
          nbcores - cpumgr.GetCurrentValueCPU() / 100 * nbcores; // just a test
      std::cout << " Number of total free cores " << ccores << std::endl;
      if (ccores < 0.3) {
        sleep(50);
      } else {
        Evaluate();
      }
#ifdef EVOLUTION
      GASimpleEvaluator::GAEvaluate();
#endif
    }
  }

  ~TGenes() {}

private:

  friend class cereal::access;

  //template <class Archive> void serialize(Archive &ar, TGenes<F> tg) { ar(tg.GetInput(), tg.GetOutput()); }

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &output;
    ar &input;
  }

public:
  bool IsDominating(const TGenes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (output[i] > other.output[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsDominated(const TGenes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (output[i] < other.output[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsEqual(const TGenes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (other.output[i] != output[i])
        return false;
    }
    return true;
  }

  void Evaluate() {
#ifdef EVOLUTION
    output = GASequentualEvaluator::Evaluate();
#endif
#ifdef ENABLE_GEANTV
    size_t sizeofOutput = sizeof(output) + sizeof(double) * output.capacity();
    const int fNumberChildren = 1;
    int pipega[fNumberChildren + 1];
    pid_t fArray[fNumberChildren];
    pid_t cpid;
    ssize_t result;
    pipe(pipega);
    double fitness;
    // Forking a child process - should be in loop too
    cpid = fork();
    // Loop if we have more children
    for (int i = 0; i < fNumberChildren; ++i) {
      //  cpid = fork();
      if (cpid > 0) {
        std::cout << "Starting father process: " << std::endl;
        fArray[i] = cpid;
        close(pipega[WRITE]);
        std::cout << "We are starting to read.." << std::endl;
        memset(&tmpoutput, 0, sizeof(tmpoutput));
        // for (int i = 0; i < output.size(); ++i) {
        // read(pipega[READ], &fitness, sizeof(double));
        while (read(pipega[READ], &fitness, sizeof(double)) > 0) {
          tmpoutput.push_back(fitness);
          std::cout << "Parent read next value: " << fitness << std::endl;
        }
        output = tmpoutput;
        std::cout << "We are stoping to read.." << std::endl;
        close(pipega[READ]);
        for (int i = 0; i < fNumberChildren; ++i) {
          std::cout << "Waiting for PID: " << fArray[i] << " to finish.."
                    << std::endl;
          waitpid(fArray[i], NULL, 0);
          std::cout << "PID: " << fArray[i] << " has shut down.." << std::endl;
        }
      } else if (cpid < 0) {
        std::cerr << "Fork for evaluation was failed." << std::endl;
        exit(EXIT_FAILURE);
      } else {
        std::cout << "Starting child.." << std::endl;
        output = F::Evaluate(input);
        close(pipega[READ]);
        for (auto it : output) {
          write(pipega[WRITE], &it, sizeof(double));
          std::cout << "Vector part to be send: " << it << std::endl;
        }
        close(pipega[WRITE]); // close the read-end of the pipe
        wait(NULL);
        exit(EXIT_SUCCESS);
      }
      std::cout << "We are back to master job::" << std::endl;
    }
    // Cleaning array from previos pids info
    std::fill(fArray, fArray + fNumberChildren, 0);
#else
    output = F::Evaluate(input);
#endif
  }

  const typename F::Input &GetInput() const { return input; }

  const typename F::Output &GetOutput() const { return output; }

  friend std::ostream &operator<<(std::ostream &os, const TGenes<F> &ind) {
    auto indvector = ind.GetInput();
    for (int i = 0; i < indvector.size(); ++i)
      os << indvector[i] << " ";
    return os;
  }
};
}

#endif
