#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class Population < Double_t > +;
#pragma link C++ class Population < double > +;
#pragma link C++ class HistogramManager + ;
#pragma link C++ class Genes < Double_t > +;
#pragma link C++ class Genes < double > +;
#pragma link C++ class AlgorithmNSGA + ;
#pragma link C++ class Functions + ;
#pragma link C++ class ExceptionMessenger + ;
#pragma link C++ class GeantVFitness + ;
#pragma link C++ class PFMWatch + ;
#pragma link C++ class ExceptionMessenger + ;

#pragma link C++ class std::vector < Double_t > +;
#pragma link C++ class std::vector < double > +;
#pragma link C++ class std::vector < std::vector < Double_t >> +;
#pragma link C++ class std::vector < std::vector < double >> +;
#endif
