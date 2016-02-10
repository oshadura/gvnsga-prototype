#ifndef __GEANTVFITNESS__
#define __GEANTVFITNESS__

#include <typeinfo>
#include <TSystem.h>
#include <vector>

class GeantVFitness {
public:
  GeantVFitness(): fMemoryVector(0) {}
  virtual ~GeantVFitness() {}
  void LogMemoryFitness();
  void LogTimeFitness();
  void HistOutputFitness();

private:
  std::vector<ProcInfo_t> fMemoryVector;

  ClassDef(GeantVFitness, 1)
};

#endif
