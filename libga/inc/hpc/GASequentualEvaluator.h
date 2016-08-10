#pragma once

#ifndef __SEQUENTUALEVALUATOR__
#define __SEQUENTUALEVALUATOR__

#include "GAEvaluate.h"
#include "generic/GADouble.h"
#include "generic/GAVector.h"
#include "tools/Random.h"

#include <cmath>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#define READ 0
#define WRITE 1

namespace geantvmoop {

class GASequentualEvaluator : public GAEvaluate<GASequentualEvaluator> {

public:
  template <typename F> static void EvaluateImpl() {
    typename F::Input input;
    typename F::Output output;
    size_t sizeofOutput = sizeof(output) + sizeof(F) * output.capacity();
    std::cout << "Size of expected buffer for fitness container is :"
              << sizeofOutput << std::endl;
    const int fNumberChildren = 1;
    int pipega[fNumberChildren + 1];
    pid_t fArrayDead[fNumberChildren];
    pid_t cpid;
    ssize_t result;
    pipe(pipega);
    cpid = fork();
    for (int i = 0; i < fNumberChildren; ++i) {
      //  cpid = fork();
      if (cpid > 0) {
        std::cout << "Starting father.." << std::endl;
        fArrayDead[i] = cpid;
        close(pipega[WRITE]);
        std::cout << "=======New fitness just created:========" << std::endl;
        for (auto i : output)
          std::cout << i << ' ' << std::endl;
        std::cout << "===============" << std::endl;
        //////////////////////////////////////
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
    // Cleaning array from previos pids info
    std::fill(fArrayDead, fArrayDead + fNumberChildren, 0);
  }
};
}

#endif
