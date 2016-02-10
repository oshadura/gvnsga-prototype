#ifndef __ALGORITHMNSGA3__
#define __ALGORITHMNSGA3__

class AlgorithmNSGA3
{
public:
	AlgorithmNSGA3(arguments);
	~AlgorithmNSGA3();

///////////////////////////////////////// MATH ////////////////////////////////////////////////////////
// ASF(): achievement scalarization function
double ASF(const std::vector<double> &objs, const std::vector<double> &weight);

// GuassianElimination(): used to calculate the hyperplane
void GuassianElimination(std::vector<double> *px, std::vector< std::vector<double> > A, const std::vector<double> &b);

// PerpendicularDistance(): calculate the perpendicular distance from a point to a line
double PerpendicularDistance(const std::vector<double> &direction, const std::vector<double> &point);

///////////////////////////////////// Reference Point ////////////////////////////////////////////////
/ ----------------------------------------------------------------------------------
//		CReferencePoint
//
// Reference points play very important roles in NSGA-III. Individuals in the population
// are associated with reference points, and the survivors in the environmental selection
// are determined based on the niche count of the reference points.
//
// Check Algorithms 1-4 in the orignal paper for the usage of reference points.
// ----------------------------------------------------------------------------------

const std::vector<double> & pos() const { return fPosition; }
std::vector<double> & pos() { return fPosition; }
std::size_t MemberSize() const { return fMmembersize; }
bool HasPotentialMember() const { return !fPotentialMembers.empty(); }
void Clear();
void AddMember();
void AddPotentialMember(std::size_t member_ind, double distance);
int FindClosestMember() const;
int RandomMember() const;
void RemovePotentialMember(std::size_t member_ind);

// ----------------------------------------------------------------------------------
// GenerateReferencePoints():
//
// Given the number of objectives (M) and the number of divisions (p), generate the set of 
// reference points. Check Section IV-B and equation (3) in the original paper.

void GenerateReferencePoints(std::vector<ReferencePoint> *rps, std::size_t M, const std::vector<std::size_t> &p);
// ----------------------------------------------------------------------------------
// Associate():
//
// Associate individuals in the population with reference points.
// Check Algorithm 3 in the original paper.
class CPopulation;
void Associate(std::vector<ReferencePoint> *prps, const Population &pop, const std::vector<std::vector<T>> v; &fronts);
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------
//	The environmental selection mechanism is the key innovation of 
//  the NSGA-III algorithm.
//
//  Check Algorithm I in the original paper of NSGA-III.
// ----------------------------------------------------------------------

void EnvironmentalSelection(Population *pnext, // population in the next generation
							Population *pcur,  // population in the current generation
							std::vector<ReferencePoint> rps, // the set of reference points
							std::size_t PopSize);

private:

	/* data */
};

#endif
