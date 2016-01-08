#ifndef __GEANTVFITNESS__
#define __GEANTVFITNESS__

#include <typeinfo>
#include <TSystem.h>
#include <vector>

class GeantVFitness{
public:
	GeantVFitness(){}
	virtual ~GeantVFitness(){}
	void LogMemoryFitness();
	void LogTimeFitness();
	void HistOutputFitness();

private:
	std::vector<ProcInfo_t> fMemoryVector;

};	

#endif
