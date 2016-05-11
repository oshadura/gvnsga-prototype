#include "generic/ReferencePoint.h"

ReferencePoint::ReferencePoint() {}

ReferencePoint::ReferencePoint(Int_t fSize) {}

ReferencePoint::~ReferencePoint() {}

void ReferencePoint::ClearRP() {
  fMemberSize = 0;
  fPossibleSolutions.clear();
}

void ReferencePoint::AddMemberRP() { fMemberSize += 1; }

void ReferencePoint::AddPotentialMemberRP(std::size_t ind, Double_t distance) {
  fPossibleSolutions.push_back(std::make_pair(ind, distance));
}

Int_t ReferencePoint::FindClosestMemberRP() {
  Double_t fMinDist = std::numeric_limits<Double_t>::max();
  Int_t fMinInd = -1;
  for (int i = 0; i < fPossibleSolutions.size(); ++i) {
    if (fPossibleSolutions[i].second < fMinDist) {
      fMinDist = fPossibleSolutions[i].second;
      fMinInd = fPossibleSolutions[i].first;
    }
  }
  return fMinInd;
}

Int_t ReferencePoint::RandomMemberRP() {
  if (fPossibleSolutions.size() > 0) {
    return fPossibleSolutions[rand() % fPossibleSolutions.size()].first;
  } else {
    return -1;
  }
}

void ReferencePoint::RemovePotentialMemberRP(std::size_t ind) {
  for (int i = 0; i < fPossibleSolutions.size(); ++i) {
    if (fPossibleSolutions[i].first == ind) {
      fPossibleSolutions.erase(fPossibleSolutions.begin() + 1);
      return;
    }
  }
}

void ReferencePoint::GenerateRecursivelyRP(std::vector<ReferencePoint> *fRP,
                                           ReferencePoint *fP,
                                           size_t fNumberObjects, size_t fLeft,
                                           size_t fTotal, size_t fElement) {
  if (fElement == fNumberObjects - 1) {
    fP->Position()[fElement] = static_cast<Double_t>(fLeft) / fTotal;
    fRP->push_back(*fP);
  } else {
    for (int i = 0; i < fLeft; ++i) {
      fP->Position()[fElement] = static_cast<Double_t>(i) / fTotal;
      GenerateRecursivelyRP(fRP, fP, fNumberObjects, (fLeft - i), fTotal,
                            (fElement + 1));
    }
  }
}

void ReferencePoint::GenerateRP(std::vector<ReferencePoint> *fRP, size_t fSize,
                                const std::vector<std::size_t> &p) {
  ReferencePoint fP(fSize); // TBD - to be implemented
  GenerateRecursivelyRP(fRP, &fP, fSize, p[0], p[0], 0);
  if (p.size() > 1) {
    std::vector<ReferencePoint> fInsideReferencePoint;
    GenerateRecursivelyRP(&fInsideReferencePoint, &fP, fSize, p[1], p[1], 0);
    Double_t fCenter = 1 / fSize;
    for (int i = 0; i < fInsideReferencePoint.size(); ++i) {
      for (int j = 0; i < fInsideReferencePoint[j].Position().size(); ++j) {
        fInsideReferencePoint[i].Position()[j] =
            (fCenter + fInsideReferencePoint[i].Position()[j]) / 2;
      }
      fRP->push_back(fInsideReferencePoint[i]);
    }
  }
}

void ReferencePoint::AssociateRP(std::vector<ReferencePoint> *fRP,
                                 const Population<Double_t> &pop,
                                 std::vector<std::vector<Int_t>> &fFront) {
  std::vector<ReferencePoint> &fRPCopy = *fRP;
  for (int t = 0; t < fFront.size(); ++t) {
    for (int i = 0; i < fFront[t].size(); ++i) {
      size_t fMinReferencePoint = fRPCopy.size();
      Double_t fMinDist = std::numeric_limits<Double_t>::max();
      for (int r = 0; r < fRPCopy.size(); ++r) {
        // Temporary just compile solution
        std::vector<double> point;
        point.push_back(1);
        Double_t d = PerpedicularDistance(
            fRPCopy[r].Position(), point /*pop[fFront[t][i] ].conv_objs()*/);
        if (d < fMinDist) {
          fMinDist = d;
          fMinReferencePoint = r;
        }
      }
      if (t + 1 != fFront.size()) {
        fRPCopy[fMinReferencePoint].AddMemberRP();
      } else {
        fRPCopy[fMinReferencePoint].AddPotentialMemberRP(fFront[t][i],
                                                         fMinDist);
      }
    }
  }
}

// ---------------------------------------------------------------------
// PerpendicularDistance:
//
// Given a direction vector (w1, w2) and a point P(x1, y1),
// we want to find a point Q(x2, y2) on the line connecting (0, 0)-(w1, w2)
// such that (x1-x2, y1-y2) is perpendicular to (w1, w2).
//
// Since Q is on the line (0, 0)-(w1, w2), it should be (w1*k, w2*k).
// (x1-w1*k, y1-w2*k).(w1, w2) = 0. (inner product)
// => k(w1^2 + w2^2) = w1x1 + w2x2
// => k = (w1x1 + w2x2)/(w1^2 +w2^2).
//
// After obtaining k, we have Q = (w1*k, w2*k) and the distance between P and Q.
//
// Code example:
//    vector<double> dir{1, 3}, point{5.5, 1.5};
//    cout << PerpendicularDistance(dir, point) << endl;
// ---------------------------------------------------------------------
Double_t
ReferencePoint::PerpedicularDistance(const std::vector<double> &fDirection,
                                     const std::vector<double> &fPoint) {
  Double_t fNumerator = 0, fDenominator = 0;
  for (size_t i = 0; i < fDirection.size(); i += 1) {
    fNumerator += fDirection[i] * fPoint[i];
    fDenominator += std::pow(fDirection[i], 2);
  }
  Double_t k = fNumerator / fDenominator;
  Double_t d = 0;
  for (size_t i = 0; i < fDirection.size(); i += 1) {
    d += std::pow((k * fDirection[i] - fPoint[i]), 2);
  }
  return sqrt(d);
}
/*
// ----------------------------------------------------------------------
// TranslateObjectives():
//
// 1. Find the ideal point
// 2. Translate the objective values
// 3. Return the ideal point
//
// Check steps 1-3 in Algorithm 2 in the original paper of NSGAIII.
// ----------------------------------------------------------------------
std::vector<Double_t> TranslateObjectives(Population<Double_t> *pop, const
std::vector<std::vector<Int_t>> &fFront)
{
        Population &p = *pop;
        std::vector<double> fIdealPoint(pop[0].objs().size());

        const size_t NumObj = pop[0].objs().size();
        for (size_t f=0; f<NumObj; f+=1)
        {
                double minf = numeric_limits<double>::max();
                for (size_t i=0; i<fronts[0].size(); i+=1) // min values must
appear in the first front
                {
                        minf = std::min(minf, pop[ fronts[0][i] ].objs()[f]);
                }
                ideal_point[f] = minf;

                for (size_t t=0; t<fronts.size(); t+=1)
                {
                        for (size_t i=0; i<fronts[t].size(); i+=1)
                        {
                                size_t ind = fronts[t][i];
                                pop[ind].conv_objs().resize(NumObj);
                                pop[ind].conv_objs()[f] = pop[ind].objs()[f] -
minf;
                        }
                }
        }

        return ideal_point;

}

// ----------------------------------------------------------------------
// FindExtremePoints():
//
// Find the extreme points along each objective axis.
// The extreme point has the minimal ASF value.
// Return the indices of extreme individuals in the population.
//
// Check step 4 in Algorithm 2 and eq. (4) in the original paper.
// ----------------------------------------------------------------------
void FindExtremePoints(std::vector<size_t> *extreme_points, const CPopulation
&pop, const CNondominatedSort::TFronts &fronts)
{
        vector<size_t> &exp = *extreme_points;
        exp.clear();

        for (size_t f=0; f<pop[0].objs().size(); f+=1)
        {
                vector<double> w(pop[0].objs().size(), 0.000001);
                w[f] = 1.0;

                double min_ASF = numeric_limits<double>::max();
                size_t min_indv = fronts[0].size();

                for (size_t i=0; i<fronts[0].size(); i+=1)  // only consider the
individuals in the first front
                {
                        double asf = MathAux::ASF(pop[ fronts[0][i]
].conv_objs(), w); // nsga3cpp 1.11 (2015.04.26 thanks to Vivek Nair for his
correction by email)

                        if ( asf < min_ASF )
                        {
                                min_ASF = asf;
                                min_indv = fronts[0][i];
                        }
                }

                exp.push_back(min_indv);
        }

}// FindExtremePoints()


// ----------------------------------------------------------------------
// FindMaxObjectives(): (added in nsga3cpp v1.2)
//
// Find the maximal objective values among the current solutions.
// Intercepts are set by these values when we cannot construct the
// hyperplane.
//
// This method follows the implementation in the following paper:
//
// Yuan Yuan, Hua Xu, Bo Wang,
// "An Experimental Investigation of Variation Operators in
//  Reference-Point Based Many-Objective Optimization,"
// GECCO, pp. 775-782, 2015
//
// I think Jerry Lee for finding out the difference between the two
// implementations and recommending the modification.
//
// ----------------------------------------------------------------------
vector<double> FindMaxObjectives(const CPopulation &pop)
{
        const size_t NumObj = pop[0].objs().size();

        vector<double> max_point(NumObj, -numeric_limits<double>::max());
        for (size_t i=0; i<pop.size(); i+=1)
        {
                for (size_t f=0; f<NumObj; f+=1)
                {
                        max_point[f] = std::max(max_point[f], pop[i].objs()[f]);
                }
        }

        return max_point;
}

// ----------------------------------------------------------------------
// ConstructHyperplane():
//
// Given the extreme points, construct the hyperplane.
// Then, calculate the intercepts.
//
// Check step 6 in Algorithm 2 in the original paper.
// ----------------------------------------------------------------------
void ConstructHyperplane(vector<double> *pintercepts, const CPopulation &pop,
const vector<size_t> &extreme_points)
{
        // Check whether there are duplicate extreme points.
        // This might happen but the original paper does not mention how to deal
with it.
        bool duplicate = false;
        for (size_t i=0; !duplicate && i<extreme_points.size(); i+=1)
        {
                for (size_t j=i+1; !duplicate && j<extreme_points.size(); j+=1)
                {
                        duplicate = (extreme_points[i] == extreme_points[j]);
                }
        }

        vector<double> &intercepts = *pintercepts;
        intercepts.assign(pop[0].objs().size(), 0);

        bool negative_intercept = false;
        if (!duplicate)
        {
                // Find the equation of the hyperplane
                vector<double> b(pop[0].objs().size(), 1.0);
                vector< vector<double> > A;
                for (size_t p=0; p<extreme_points.size(); p+=1)
                {
                        A.push_back(pop[ extreme_points[p] ].conv_objs()); //
v1.11: objs() -> conv_objs()
                }
                vector<double> x;
                MathAux::GuassianElimination(&x, A, b);

                // Find intercepts
                for (size_t f=0; f<intercepts.size(); f+=1)
                {
                        intercepts[f] = 1.0/x[f];

                        if(x[f] < 0)
                        {
                                negative_intercept = true;
                                break;
                        }
                }
        }

        if (duplicate || negative_intercept) // v1.2: follow the method in Yuan
et al. (GECCO 2015)
        {
                vector<double> max_objs = FindMaxObjectives(pop);
                for (size_t f=0; f<intercepts.size(); f+=1)
                {
                        intercepts[f] = max_objs[f];
                }
        }
}

// ----------------------------------------------------------------------
// NormalizeObjectives():
//
// Normalize objective values with respect to the intercepts and the ideal
point.
// Check step  7 in Algorithm 2 and eq. (5) in the original paper.
// ----------------------------------------------------------------------
void NormalizeObjectives(CPopulation *ppop, const CNondominatedSort::TFronts
&fronts, const vector<double> &intercepts, const vector<double> &ideal_point)
{
        CPopulation &pop = *ppop;

        for (size_t t=0; t<fronts.size(); t+=1)
        {
                for (size_t i=0; i<fronts[t].size(); i+=1)
                {
                        size_t ind = fronts[t][i];
                        for (size_t f=0; f<pop[ ind ].conv_objs().size(); f+=1)
                        {
                                if ( fabs(intercepts[f])>10e-10 ) // avoid the
divide-by-zero error
                                        pop[ ind ].conv_objs()[f] = pop[ ind
].conv_objs()[f]/(intercepts[f]); // v1.11: fixed
                                else
                                        pop[ ind ].conv_objs()[f] = pop[ ind
].conv_objs()[f]/10e-10;
                        }
                }
        }

}// NormalizeObjectives()

// ----------------------------------------------------------------------
// FindNicheReferencePoint():
//
// Find the reference point with the minimal cluster size.
// Return one randomly if there is more than one point.
//
// Check steps 3-4 in Algorithm 4 in the original paper.
// ----------------------------------------------------------------------
size_t FindNicheReferencePoint(const vector<CReferencePoint> &rps)
{
        // find the minimal cluster size
        size_t min_size = numeric_limits<size_t>::max();
        for (size_t r=0; r<rps.size(); r+=1)
        {
                min_size = std::min(min_size, rps[r].MemberSize());
        }

        // find the reference points with the minimal cluster size Jmin
        vector<size_t> min_rps;
        for (size_t r=0; r<rps.size(); r+=1)
        {
                if (rps[r].MemberSize() == min_size)
                {
                        min_rps.push_back(r);
                }
        }

        // return a random reference point (j-bar)
        return min_rps[rand()%min_rps.size()];
}

// ----------------------------------------------------------------------
// SelectClusterMember():
//
// Select a potential member (an individual in the front Fl) and associate
// it with the reference point.
//
// Check the last two paragraphs in Section IV-E in the original paper.
// ----------------------------------------------------------------------
int SelectClusterMember(const CReferencePoint &rp)
{
        int chosen = -1;
        if (rp.HasPotentialMember())
        {
                if (rp.MemberSize() == 0) // currently has no member
                {
                        chosen =  rp.FindClosestMember();
                }
                else
                {
                        chosen =  rp.RandomMember();
                }
        }
        return chosen;
}

// ----------------------------------------------------------------------
// EnvironmentalSelection():
//
// Check Algorithms 1-4 in the original paper.
// ----------------------------------------------------------------------
void EnvironmentalSelection(CPopulation *pnext, CPopulation *pcur,
vector<CReferencePoint> rps, size_t PopSize)
{
        CPopulation &cur = *pcur, &next = *pnext;
        next.clear();

        // ---------- Step 4 in Algorithm 1: non-dominated sorting ----------
        CNondominatedSort::TFronts fronts = NondominatedSort(cur);

        // ---------- Steps 5-7 in Algorithm 1 ----------
        vector<size_t> considered; // St
        size_t last = 0, next_size = 0;
        while (next_size < PopSize)
        {
                next_size += fronts[last].size();
                last += 1;
        }
        fronts.erase(fronts.begin()+last, fronts.end()); // remove useless
individuals

        for (size_t t=0; t<fronts.size()-1; t+=1)
        {
                for (size_t i=0; i<fronts[t].size(); i+=1)
                {
                        next.push_back(cur[ fronts[t][i] ]);
                }
        }

        // ---------- Steps 9-10 in Algorithm 1 ----------
        if (next.size() == PopSize) return;


        // ---------- Step 14 / Algorithm 2 ----------
        vector<double> ideal_point = TranslateObjectives(&cur, fronts);

        vector<size_t> extreme_points;
        FindExtremePoints(&extreme_points, cur, fronts);

        vector<double> intercepts;
        ConstructHyperplane(&intercepts, cur, extreme_points);

        NormalizeObjectives(&cur, fronts, intercepts, ideal_point);

        // ---------- Step 15 / Algorithm 3, Step 16 ----------
        Associate(&rps, cur, fronts);

        // ---------- Step 17 / Algorithm 4 ----------
        while (next.size() < PopSize)
        {
                size_t min_rp = FindNicheReferencePoint(rps);

                int chosen = SelectClusterMember(rps[min_rp]);
                if (chosen < 0) // no potential member in Fl, disregard this
reference point
                {
                        rps.erase(rps.begin()+min_rp);
                }
                else
                {
                        rps[min_rp].AddMember();
                        rps[min_rp].RemovePotentialMember(chosen);
                        next.push_back(cur[chosen]);
                }
        }

}
// ----------------------------------------------------------------------
*/
