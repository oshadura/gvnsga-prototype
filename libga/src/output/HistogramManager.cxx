#include <vector>

#include "generic/TGenes.h"
#include "generic/Population.h"
#include "algorithms/GANSGA2.h"
#include "output/HistogramManager.h"

#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF2.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "Rtypes.h"

using namespace geantvmoop;

// ClassImpT(HistogramManager,F)
templateClassImp(HistogramManager);

// Moved to .h...
/*
template <class F>
bool HistogramManager<F>::HistoFill(Population<F> &pop, char *hfile) {
  double ind, fitness;
  TH1F *PopDist =
      new TH1F("PopDist", "Population distribution", pop.size(), 0., 1.);
  PopDist->GetXaxis()->SetTitle("Population distribution");
  PopDist->GetYaxis()->SetTitle("N of bins");

  TH1F *PopFitnessDist = new TH1F(
      "PopFitnessDist", "Population fitness distribution", pop.size(), 0., 1.);
  PopFitnessDist->GetXaxis()->SetTitle("Population fitness distribution");
  PopFitnessDist->GetYaxis()->SetTitle("N of bins");

  TH2F *XScatter = new TH2F("XScatter", "N events versus size vector", 50, 0.,
                            1., 50, 0., 1.);
  XScatter->GetXaxis()->SetTitle("N Events");
  XScatter->GetYaxis()->SetTitle("Size of vector");

  for (int i = 0; i < pop.size(); ++i) {
    // Bad! I need to have one more variable!
    for (int j = 0; j < pop[0].size(); ++j) {
      ind = pop.GetIJGenes(i, j);
      fitness = pop.GetIJFitness(i, j);
      PopDist->Fill(ind);
      PopFitnessDist->Fill(fitness);
    }
  }
  TCanvas *output = new TCanvas("output", "DTLZ1 genetic optimisation", 1);
  output->Divide(2, 2);
  output->cd(1);
  PopDist->Draw();
  output->cd(2);
  XScatter->Draw();
  output->cd(3);
  PopFitnessDist->Draw();
  return true;
}
*/

template <class F> void HistogramManager<F>::Reset() { /*Instance = 0;*/
}
