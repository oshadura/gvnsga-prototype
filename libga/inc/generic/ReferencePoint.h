//===--- ReferencePoint.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ReferencePoint.h
 * @brief Implementation of RP class for LibGA
 * prototype
 */
//
#pragma once

#ifndef __REFERENCEPOINT__
#define __REFERENCEPOINT__

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "generic/Population.h"
#include "generic/TGenes.h"

#include "TObject.h"

namespace geantvmoop {

class ReferencePoint {
public:
  ReferencePoint();
  ReferencePoint(Int_t fSize);
  ~ReferencePoint();
  void ClearRP();
  void AddMemberRP();
  void AddPotentialMemberRP(std::size_t ind, Double_t distance);
  Int_t FindClosestMemberRP();
  Int_t RandomMemberRP();
  void RemovePotentialMemberRP(std::size_t ind);
  void GenerateRecursivelyRP(std::vector<ReferencePoint> *fRP,
                             ReferencePoint *p, size_t fNumberObjects,
                             size_t fLeft, size_t fTotal, size_t fElement);
  void GenerateRP(std::vector<ReferencePoint> *fRP, size_t fSize,
                  const std::vector<std::size_t> &p);
  void AssociateRP(std::vector<ReferencePoint> *fRP,
                   const Population<Double_t> &pop,
                   std::vector<std::vector<Int_t>> &fFront);
  Double_t PerpedicularDistance(const std::vector<double> &fDirection,
                                const std::vector<double> &fPoint);
  std::vector<Double_t> Position() const { return fPosition; }

private:
  std::vector<Double_t> fPosition;
  std::vector<std::pair<std::size_t, Double_t>> fPossibleSolutions;
  std::size_t fMemberSize;
};
}

#endif
