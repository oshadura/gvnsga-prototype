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

#define READ 0
#define WRITE 1

namespace geantvmoop {

template <typename F> class TGenes;
template <typename F> using individual_t = std::shared_ptr<TGenes<F> >;

template <typename F> class TGenes {

private:
  typename F::Input input;
  typename F::Output output;
  typename F::Output tmpoutput;

public:
  TGenes() {}

  TGenes(const typename F::Input &i, bool fEvaluated = true) : input(i) {
    if (fEvaluated) {
      /*
      CPUManager cpumgr;
      cpumgr.InitCPU();
      hwloc_topology_t topology;
      int nbcores, ccores;
      unsigned int microseconds;
      microseconds = 3000;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      // printf(" Number of free cores %d cores\n", nbcores);
      hwloc_topology_destroy(topology);
      ccores =
          nbcores - cpumgr.GetCurrentValue() / 100 * nbcores; // just a test
      printf(" Number of free cores %d cores\n", ccores);
      if (ccores < 0.3) {
        usleep(microseconds);
      } else {
        */
        Evaluate();
//}
#ifdef EVOLUTION
        GASimpleEvaluator::GAEvaluate();
#endif
      //}
    }
  }

  ~TGenes() {}

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
    std::cout << "Size of expected buffer for fitness container is :"
              << sizeofOutput << std::endl;
    const int fNumberChildren = 1;
    int pipega[fNumberChildren + 1];
    pid_t fArrayDead[fNumberChildren];
    pid_t cpid;
    ssize_t result;
    pipe(pipega);
    /*
    cpid = fork();
    for (int i = 0; i < fNumberChildren; ++i) {
      //  cpid = fork();
      if (cpid > 0) {
        std::cout << "Starting father.." << std::endl;
        fArrayDead[i] = cpid;
        close(pipega[WRITE]);
        std::cout << "=======New fitness just created:========" << std::endl;
        for (auto i : tmpoutput)
          std::cout << i << ' ' << std::endl;
        std::cout << "===============" << std::endl;
        std::cout << "We are starting to read.." << std::endl;
        while (read(pipega[READ], &output, sizeofOutput * 2) > 0) {
          std::cout << "=======Parent reads:========" << std::endl;
          for (auto i : output)
            std::cout << i << ' ' << std::endl;
          std::cout << "===============" << std::endl;
        }
        std::cout << "===============" << std::endl;
        std::cout << "We are stoping to read.." << std::endl;
        close(pipega[READ]);
        for (int i = 0; i < fNumberChildren; ++i) {
          std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.."
                    << std::endl;
          waitpid(fArrayDead[i], NULL, 0);
          std::cout << "PID: " << fArrayDead[i] << " has shut down.."
                    << std::endl;
        }
      } else if (cpid < 0) {
        std::cerr << "Fork for evaluation was failed." << std::endl;
        exit(EXIT_FAILURE);
      } else {
        std::cout << "Starting child.." << std::endl;
        output = F::Evaluate(input);
        close(pipega[READ]);
        write(pipega[WRITE], &output, sizeofOutput);
        std::cout << "=======Child writes:========" << std::endl;
        std::cout << "===============" << std::endl;
        close(pipega[WRITE]); // close the read-end of the pipe
        wait(NULL);
        exit(EXIT_SUCCESS);
      }
      std::cout << "We are back to master job::" << std::endl;
    }
    */
    double fitness; 
    // Forking a child process - should be in loop too
    cpid = fork();
    // Loop if we have more children
    for (int i = 0; i < fNumberChildren; ++i) {
      //  cpid = fork();
      if (cpid > 0) {
        std::cout << "Starting father process: " << std::endl;
        fArrayDead[i] = cpid;
        close(pipega[WRITE]);
        std::cout << "We are starting to read.." << std::endl;
        memset(&tmpoutput, 0, sizeof(tmpoutput));
        //for (int i = 0; i < output.size(); ++i) {
          //read(pipega[READ], &fitness, sizeof(double));
      while (read(pipega[READ], &fitness, sizeofOutput * 2) > 0) {
          tmpoutput.push_back(fitness);
          std::cout << "Parent read next value: " << fitness << std::endl;
        }
        output = tmpoutput;
        std::cout << "We are stoping to read.." << std::endl;
        close(pipega[READ]);
        for (int i = 0; i < fNumberChildren; ++i) {
          std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.."
                    << std::endl;
          waitpid(fArrayDead[i], NULL, 0);
          std::cout << "PID: " << fArrayDead[i] << " has shut down.."
                    << std::endl;
        }
      } else if (cpid < 0) {
        std::cerr << "Fork for evaluation was failed." << std::endl;
        exit(EXIT_FAILURE);
      } else {
        std::cout << "Starting child.." << std::endl;
        output = F::Evaluate(input);
        close(pipega[READ]);
        memset(&tmpoutput, 0, sizeof(tmpoutput));
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
    std::fill(fArrayDead, fArrayDead + fNumberChildren, 0);
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
