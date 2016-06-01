#ifndef __HISTOGRAMMANAGER__
#define __HISTOGRAMMANAGER__

//===--- HistogramManager.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file HistogramManager.h
 * @brief Implementation of CSV class for LibGA
 * prototype
 */
//

#include "generic/Population.h"
#include "TObject.h"
#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF2.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "Rtypes.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <string>

#include <stdio.h>
#include <math.h>

using namespace geantvmoop;

template <typename F> class HistogramManager : public TObject {

private:
  F problem;

public:
  static HistogramManager &GetInstance() {
    static HistogramManager Instance;
    return Instance;
  }

  HistogramManager(HistogramManager const &) = delete;

  void operator=(HistogramManager const &) = delete;

  std::vector<int> ScatterCombinationCalculator(int N, int K) {
    std::vector<int> ScattComb;
    ScattComb.reserve(N);
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0);      // N-K trailing 0's
    // print integers and permute bitmask
    do {
      for (int i = 1; i <= N; ++i) {
        if (bitmask[i]) {
          // std::cout << " " << i;
          ScattComb.push_back(i);
        }
      }
      // std::cout << std::endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return ScattComb;
  }

  bool HistoFill(Population<F> &pop, char *hfile, int generation) {
    char namepop[20], namefitn[20], namefolder[20], namescatter[20], x1str[10],
        x2str[10];
    std::vector<int> ScatterCombination;
    // TDirectory *folder;
    TObjArray HList(0);
    sprintf(namepop, "%s%d", "PopDist", generation);
    sprintf(namefitn, "%s%d", "PopFitnessDist", generation);
    sprintf(namefolder, "%s%d", "PopulationStatisticsGeneration", generation);
    TFile file(hfile, "update");
    file.mkdir(namefolder);
    file.cd(namefolder);
    /////////////////////////////////////////////////////////////////
    TH1F *PopDist =
        new TH1F(namepop, "Population distribution", pop.size(), 0., 1.);
    PopDist->GetXaxis()->SetTitle("TGenes / bins");
    PopDist->GetYaxis()->SetTitle("N");
    /////////////////////////////////////////////////////////////////
    TH1F *PopFitnessDist = new TH1F(namefitn, "Population fitness distribution",
                                    pop.size(), 0., 1.);
    PopFitnessDist->GetXaxis()->SetTitle("TGenes / bins");
    PopFitnessDist->GetYaxis()->SetTitle("N");
    /////////////////////////////////////////////////////////////////
    std::cout << "Lets check Scatter combination vector" << std::endl;
    ScatterCombination =
        ScatterCombinationCalculator(pop.GetTGenes(0).size() + 1, 2);
    for (auto i : ScatterCombination)
      std::cout << i << ' ';
    /////////////////////////////////////////////////////////////////

    for (int i = 0; i < pop.size(); ++i) {
      // X Scatter plots
      for (int it = 0; it < ScatterCombination.size(); it+=2) {
        //Taking correct X ID
        auto valueX1 = ScatterCombination.at(it);
        auto valueX2 = ScatterCombination.at(it+1);

        sprintf(namescatter, "%s%d%s%d", "X", valueX1, "vsX", valueX2);
        sprintf(x1str, "%s%d", "X", valueX1);
        sprintf(x2str, "%s%d", "X", valueX2);

        auto x1 = pop.GetGeneValue(i, it);
        auto x2 = pop.GetGeneValue(i, it + 1);
        TString histoname = TString::Format(namescatter);
        TH2F *myhist = ((TH2F *)(HList.FindObject(histoname)));
        if (!myhist) {
          myhist = new TH2F(TString::Format(namescatter),
                            "Scatter plot of different TGenes", pop.size(), 0.,
                            1., pop.size(), 0., 1.);
          HList.Add(myhist);
        }
        myhist->GetXaxis()->SetTitle(x1str);
        myhist->GetYaxis()->SetTitle(x2str);
        myhist->Fill(x1, x2);
      }
      // Distribution plots
      for (int j = 0; j < pop.GetTGenes(0).size(); ++j) {
        auto ind = pop.GetGeneValue(i, j);
        auto fitness = pop.GetObjectiveValue(i, j);
        PopDist->Fill(ind);
        PopFitnessDist->Fill(fitness);
      }
    }
    PopDist->Draw();
    PopDist->Write();
    HList.Write();
    // XScatter->Draw();
    // XScatter->Write();
    PopFitnessDist->Draw();
    PopFitnessDist->Write();

    return true;
  }

  void Reset();

private:
  HistogramManager() {}
  ~HistogramManager() {}

  ClassDef(HistogramManager, 1)
};

// ClassDefT2(HistogramManager,F)

#endif
