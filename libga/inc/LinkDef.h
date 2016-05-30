#ifdef __ROOTCLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;


#pragma link C++ namespace geantvmoop;

#pragma link C++ class GACD + ;
#pragma link C++ class GAComparator + ;
#pragma link C++ class GAGenericIndicator + ;
#pragma link C++ class GANDRank + ;

#pragma link C++ class CMAES + ;
#pragma link C++ class GAMOEAD + ;
#pragma link C++ class GANSGA2 + ;
#pragma link C++ class GANSGA3 + ;

#pragma link C++ class GACrossover + ;
#pragma link C++ class GAMutation + ;
#pragma link C++ class GAPolMutation + ;
#pragma link C++ class GASBXCrossover + ;
#pragma link C++ class GASelection + ;
#pragma link C++ class GASimpleSelection + ;
#pragma link C++ class GATournamentSelection + ;
#pragma link C++ class GANoiseReduction + ;
#pragma link C++ class PCAinvPCA + ;

#pragma link C++ class Functions + ;
#pragma link C++ class FunctiosProxy + ;
#pragma link C++ class GAAlgorithm + ;
#pragma link C++ class GAConstrainedValue + ;
#pragma link C++ class GADouble + ;
#pragma link C++ class GAValue + ;
#pragma link C++ class GAVector + ;
#pragma link C++ class PF + ;
#pragma link C++ class Population + ;
#pragma link C++ class ReferencePoint + ;
#pragma link C++ class TGenes + ;
#pragma link C++ class TWeight + ;

#pragma link C++ class GAEvaluate + ;
#pragma link C++ class GAMPIEvaluator + ;
#pragma link C++ class GASequentualEvaluator + ;
#pragma link C++ class MPIClient + ;
#pragma link C++ class MPIEvaluateClient + ;
#pragma link C++ class MPIEvaluateServer + ;
#pragma link C++ class MPIOrganization + ;
#pragma link C++ class MPIServer + ;

#pragma link C++ class KPCA + ;
#pragma link C++ class LPCA + ;
#pragma link C++ class PCA + ;

// csv.h
#pragma link C++ class CSVManager + ;

#pragma link C++ class HistogramManager<RunGeantV> + ;
#pragma link C++ class HistogramManager<DTLZ1> + ;
#pragma link C++ class HistogramManager<DTLZ2> + ;
#pragma link C++ class HistogramManager<DTLZ3> + ;
#pragma link C++ class HistogramManager<DTLZ4> + ;
#pragma link C++ class HistogramManager<DTLZ5> + ;
#pragma link C++ class HistogramManager<DTLZ6> + ;
#pragma link C++ class HistogramManager<DTLZ7> + ;


#pragma link C++ class DTLZ1;
#pragma link C++ class DTLZ2;
#pragma link C++ class DTLZ3;
#pragma link C++ class DTLZ4;
#pragma link C++ class DTLZ5;
#pragma link C++ class DTLZ6;
#pragma link C++ class DTLZ7;

#pragma link C++ class RunGeantV + ;

#pragma link C++ class GeantVFitness + ;

#pragma link C++ class GASort + ;
#pragma link C++ class Random + ;

#endif
