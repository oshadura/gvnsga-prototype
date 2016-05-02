#ifndef __POLYNOMIALMUTATION__
#define __POLYNOMIALMUTATION__

#include "Mutation.h"

class PolynomialMutation : public Mutation<PolynomialMutation> {

public:
  template <typename F, typename T>
  static void MutateGA(const Functions *setup) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rand(0, 1);
    Double_t fRnd, fDelta1, fDelta2, fMutPow, fDelta, fValue;
    Double_t y, LimitDown, LimitUp, xy;
    Int_t fNMut = 0;
    for (Int_t j = 0; j < setup->fNParam; ++j) {
      if (rand(gen) <= setup->fPMut) {
        y = fGenes[j];
        LimitDown = setup->GetIntervalLimit(j).first;
        LimitUp = setup->GetIntervalLimit(j).second;
        fDelta1 = (y - LimitDown) / (LimitUp - LimitDown);
        fDelta2 = (LimitUp - y) / (LimitUp - LimitDown);
        fRnd = rand(gen);
        fMutPow = 1.0 / (setup->fEtaMut + 1.0);
        if (fRnd <= 0.5) {
          xy = 1.0 - fDelta1;
          fValue = 2.0 * fRnd +
                   (1.0 - 2.0 * fRnd) * (pow(xy, (setup->fEtaMut + 1.0)));
          fDelta = pow(fValue, fMutPow) - 1.0;
        } else {
          xy = 1.0 - fDelta2;
          fValue = 2.0 * (1.0 - fRnd) +
                   2.0 * (fRnd - 0.5) * (pow(xy, (setup->fEtaMut + 1.0)));
          fDelta = 1.0 - (pow(fValue, fMutPow));
        }
        y = y + fDelta * (LimitUp - LimitDown);
        // std::cout << "New part of Gene<T> for Mutation() " << y << std::endl;
        if (y < LimitDown)
          y = LimitDown;
        if (y > LimitUp)
          y = LimitUp;
        SetGene(j, y);
        fNMut = fNMut + 1;
      }
    }
    // std::cout << "Print number of mutations in Gene: " << fNMut << std::endl;
    return fNMut;
  }
};

#endif
