#ifndef __SBXCROSSOVER__
#define __SBXCROSSOVER__

#include "Crossover.h"
#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/GaVector.h"
#include "generic/RandomDouble.h"

#include "random"

namespace geantvmoop {

class SBXCrossover : public Crossover<SBXCrossover> {

public:
  template <typename Individual>
  static Individual CrossoverGA(const Individual &ind1,
                                const Individual &ind2) {
    GaVector<RandomDouble> *offspring1 = static_cast<Individual &>(ind1);
    GaVector<RandomDouble> *offspring2 = static_cast<Individual &>(ind2);
    for (unsigned int i = 0; i < offspring1->size(); ++i) {
      SBXCrossover::CrossoverEvolution(offspring1[i], offspring2[i], 0.5);
    }
    return offspring1;
  }

  static void CrossoverEvolution(RandomDouble &a, RandomDouble &b,
                                 double fCrossDistributionIndex) {
    double fX0 = a.GetValue();
    double fX1 = b.GetValue();
    double fDeltaX = fabs(fX1 - fX0);
    double fLowBound = a.GetDownBound();
    double fUpBound = b.GetUpBound();
    double fBoundedValue1, fBoundedValue2, fnewa, fnewb, b1, b2;
    if (fX0 < fX1) {
      fBoundedValue1 = 1 + 2 * (fX0 - fLowBound) / fDeltaX;
      fBoundedValue2 = 1 + 2 * (fUpBound - fX1) / fDeltaX;
    } else {
      fBoundedValue1 = 1 + 2 * (fX1 - fLowBound) / fDeltaX;
      fBoundedValue2 = 1 + 2 * (fUpBound - fX0) / fDeltaX;
    }
    if (fBoundedValue1 < fBoundedValue2) {
      fBoundedValue2 = fBoundedValue1;
    } else {
      fBoundedValue1 = fBoundedValue2;
    }
    double fRandomProbUniform1 =
        1 - 1 / (2 * std::pow(fBoundedValue1, fCrossDistributionIndex + 1));
    double fRandomProbUniform2 =
        1 - 1 / (2 * std::pow(fBoundedValue2, fCrossDistributionIndex + 1));
    double fRandomProb = Generator::GetInstance().RNGDouble();
    if (fRandomProb == 1.0) {
      fRandomProb = std::nextafter(fRandomProb, -1);
    }
    double u1 = fRandomProb * fRandomProbUniform1;
    double u2 = fRandomProb * fRandomProbUniform2;
    if (u1 <= 0.5) {
      b1 = std::pow(2 * u1, 1 / (fCrossDistributionIndex + 1));
    } else {
      b1 = std::pow(0.5 / (1 - u1), 1 / (fCrossDistributionIndex + 1));
    }
    if (u2 <= 0.5) {
      b2 = std::pow(2 * u2, 1 / (fCrossDistributionIndex + 1));
    } else {
      b2 = std::pow(0.5 / (1 - u2), 1 / (fCrossDistributionIndex + 1));
    }
    if (fX0 < fX1) {
      fnewa = 0.5 * (fX0 + fX1 + b1 * (fX0 - fX1));
      fnewb = 0.5 * (fX0 + fX1 + b2 * (fX1 - fX0));
    } else {
      fnewa = 0.5 * (fX0 + fX1 + b2 * (fX0 - fX1));
      fnewb = 0.5 * (fX0 + fX1 + b1 * (fX1 - fX0));
    }
    if (fnewa < fLowBound) {
      a.SetValue(fLowBound);
    } else if (fnewa > fUpBound) {
      a.SetValue(fUpBound);
    } else {
      a.SetValue(fnewa);
    }
    if (fnewb < fLowBound) {
      b.SetValue(fLowBound);
    } else if (fnewb > fUpBound) {
      b.SetValue(fUpBound);
    } else {
      b.SetValue(fnewb);
    }
    if (Generator::GetInstance().RNGBool()) {
      double temp = a.GetValue();
      a.SetValue(b.GetValue());
      b.SetValue(temp);
    }
  }
};
}

#endif
