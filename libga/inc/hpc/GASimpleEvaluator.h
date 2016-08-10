#pragma once

#ifndef __SIMPLEEVALUATOR__
#define __SIMPLEEVALUATOR__

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

namespace geantvmoop {

class GASimpleEvaluator : public GAEvaluate<GASimpleEvaluator> {

public:
  template <typename F> static void EvaluateImpl() {
    typename F::Input input;
    typename F::Output output;
    output = F::Evaluate(input);
    return output;
  }
};
}

#endif