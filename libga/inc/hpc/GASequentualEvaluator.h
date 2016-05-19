#ifndef __SEQUENTUALEVALUATOR__
#define __SEQUENTUALEVALUATOR__

#include "GACrossover.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"
#include "tools/Random.h"
#include <cmath>

namespace geantvmoop {

class SequentualEvaluator : public GAEvaluate<SequentualEvaluator> {

public:
  // void Evaluate() { output = F::Evaluate(input); }
  template <typename T> static void Evaluate() {
	// Changes according new function ... -> recheck.
    size_t sizeofOutput = sizeof(fOutput) + sizeof(T) * fOutput.capacity();
    std::cout << "Size of expected buffer for fitness is :" << sizeofOutput
              << std::endl;
    const int fNumberChildren = 1;
    int pipeGA[fNumberChildren + 1];
    // Array of pids to be killed after
    pid_t fArrayDead[fNumberChildren];
    // Pid
    pid_t cpid;
    ssize_t result;
    pipe(pipeGA);
    ///////////////////////////////////////////////////////////////////////
    std::vector<T> tempFitness;
    tempFitness.resize(ind.setup->fNObjectives);
    // Forking a child process - should be in loop too
    cpid = fork();
    // Loop if we have more children
    for (int i = 0; i < fNumberChildren; ++i) {
      //  cpid = fork();
      if (cpid > 0) {
        std::cout << "Starting father.." << std::endl;
        fArrayDead[i] = cpid;
        close(pipeGA[WRITE]);
        std::cout << "=======New fitness just created:========" << std::endl;
        for (auto i : tempFitness)
          std::cout << i << ' ' << std::endl;
        std::cout << "===============" << std::endl;
        //////////////////////////////////////
        std::cout << "We are starting to read.." << std::endl;
        while (read(pipeGA[READ], &tempFitness, sizeofOutput * 2) > 0) {
          std::cout << "=======Parent reads:========" << std::endl;
          for (auto i : tempFitness)
            std::cout << i << ' ' << std::endl;
          std::cout << "===============" << std::endl;
          ind.SetFitness(tempFitness);
        }
        std::cout << "===============" << std::endl;
        std::cout << "We are stoping to read.." << std::endl;
        close(pipeGA[READ]);
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
        ind.SetFitness((setup.evfunc)(ind));
        close(pipeGA[READ]);
        tempFitness = ind.GetFitnessVector();
        write(pipeGA[WRITE], &tempFitness, sizeofOutput);
        std::cout << "=======Child writes:========" << std::endl;
        for (auto it = ind.GetFitnessVector().begin();
             it != ind.GetFitnessVector().end(); ++it) {
          std::cout << *it << std::endl;
        }
        std::cout << "===============" << std::endl;
        close(pipeGA[WRITE]); // close the read-end of the pipe
        wait(NULL);
        exit(EXIT_SUCCESS);
      }
      std::cout << "We are back to master job::" << std::endl;
    }
    if (setup.fNCons) {
      ConstViol = 0;
    }
    fEvaluated = true;
    // Cleaning array of previos pids
    fArrayDead[fNumberChildren] = 0;
    return fOutput;
  }
};

}