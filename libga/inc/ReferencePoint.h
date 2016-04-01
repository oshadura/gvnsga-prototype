#ifndef __REFERENCEPOINT__
#define __REFERENCEPOINT__

#include <vector>
#include <utility>
#include <limits>

#include "TObject.h"

class ReferencePoint{
public:
	ReferencePoint();
	~ReferencePoint();
	void ClearRP();
	void AddMemberRP();
	void AddPotentialMemberRP(std::size_t ind, Double_t distance);
	Int_t FindClosestMemberRP();
	Int_t RandomMemberRP();
	void RemovePotentialMemberRP(std::size_t ind);
private:
	std::vector<Double_t> fPosition;
	std::vector<std::pair<std::size_t,Double_t>> fPossibleSolutions;
	std::size_t fMemberSize;
};

#endif