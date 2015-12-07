#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class Population<Double_t> + ;
#pragma link C++ class HistogramManager + ;
#pragma link C++ class Genes<Double_t> + ;
#pragma link C++ class AlgorithmNSGA + ;
#pragma link C++ class Functions + ;

#pragma link C++ class vector < Double_t > +;
#pragma link C++ class std::vector < std::vector < Double_t >> +;

#endif
