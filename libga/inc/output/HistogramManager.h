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
#include "TStyle.h"
#include "TGraph2D.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "Rtypes.h"
#include "TF3.h"
#include "TError.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"

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
      for (int i = 0; i < N; ++i) {
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
    std::cout << "Building output statistics for generation " << generation
              << std::endl;
    char namepop[20], namefitn[20], namefolder[20], namescatter[20], x1str[10],
        x2str[10], y1str[10], y2str[10], nameFitLand[20];
    std::vector<int> ScatterCombinationX, ScatterCombinationY;
    // TDirectory *folder;
    TObjArray HXList(0);
    TObjArray HYList(0);
    sprintf(namepop, "%s%d", "PopDist", generation);
    sprintf(namefitn, "%s%d", "PopFitnessDist", generation);
    sprintf(namefolder, "%s%d", "PopulationStatisticsGeneration", generation);
    sprintf(nameFitLand, "%s%d", "FitLand", generation);
    TFile file(hfile, "update");
    file.mkdir(namefolder);
    file.cd(namefolder);

    //////// Population distribution
    TH1F *PopDist =
        new TH1F(namepop, "Population distribution", pop.size(), 0., 1.);
    PopDist->GetXaxis()->SetTitle("TGenes / bins");
    PopDist->GetYaxis()->SetTitle("N");

    /////// Fitness landscape
    TF3 *FitLand = new TF3(nameFitLand, "[0] * x + [1] * y  + [3] * z - 0.5",
                           0, 100, 0, 100, 0, 100);
    FitLand->SetParameters(2, 2, 2);
    ROOT::Fit::Fitter fitter;
    // wrapped the TF1 in a IParamMultiFunction interface for teh Fitter class
    ROOT::Math::WrappedMultiTF1 wrapperfunction(*FitLand, 3);
    fitter.SetFunction(wrapperfunction);
    ROOT::Fit::BinData data(pop.size(), 3);
    double function[3];
    double predictor[pop.size()];
    double error = 0.1;
    ////// Preparation for scatter combination
    // std::cout << "Lets check Scatter combination vector" << std::endl;
    ScatterCombinationX =
        ScatterCombinationCalculator(pop.GetTGenes(0).size(), 2);
    ScatterCombinationY =
        ScatterCombinationCalculator(pop.GetTFitness(0).size(), 2);

    //////// Fitness distribution
    TH1F *PopFitnessDist = new TH1F(namefitn, "Population fitness distribution",
                                    pop.size(), 0., 1.);
    PopFitnessDist->GetXaxis()->SetTitle("TGenes / bins");
    PopFitnessDist->GetYaxis()->SetTitle("N");

    /////////////////////////////////////////////////////////////////
    /*
    for (auto i : ScatterCombination)
      std::cout << i << ' ';
    std::cout << std::endl;
    */
    /////////////////////////////////////////////////////////////////
    for (int i = 0; i < pop.size(); ++i) {
      // X Scatter plots
      for (int it = 0; it < ScatterCombinationX.size(); it += 2) {
        // Taking correct X ID
        auto valueX1 = ScatterCombinationX.at(it);
        auto valueX2 = ScatterCombinationX.at(it + 1);

        sprintf(namescatter, "%s%d%s%d", "X", valueX1, "vsX", valueX2);
        sprintf(x1str, "%s%d", "X", valueX1);
        sprintf(x2str, "%s%d", "X", valueX2);

        auto x1 = pop.GetGeneValue(i, it);
        auto x2 = pop.GetGeneValue(i, it + 1);
        TString histoname = TString::Format(namescatter);
        TH2F *myhistx = ((TH2F *)(HXList.FindObject(histoname)));
        if (!myhistx) {
          myhistx = new TH2F(TString::Format(namescatter),
                             "Scatter plot of different TGenes", pop.size(), 0.,
                             1., pop.size(), 0., 1.);
          HXList.Add(myhistx);
        }
        myhistx->GetXaxis()->SetTitle(x1str);
        myhistx->GetYaxis()->SetTitle(x2str);
        myhistx->Fill(x1, x2);
      }
      // Y Scatter plots
      for (int it = 0; it < ScatterCombinationY.size(); it += 2) {
        // Taking correct Y ID
        auto valueY1 = ScatterCombinationY.at(it);
        auto valueY2 = ScatterCombinationY.at(it + 1);

        sprintf(namescatter, "%s%d%s%d", "Y", valueY1, "vsY", valueY2);
        sprintf(y1str, "%s%d", "Y", valueY1);
        sprintf(y2str, "%s%d", "Y", valueY2);

        auto y1 = pop.GetObjectiveValue(i, it);
        auto y2 = pop.GetObjectiveValue(i, it + 1);
        TString histoname = TString::Format(namescatter);
        TH2F *myhisty = ((TH2F *)(HYList.FindObject(histoname)));
        if (!myhisty) {
          myhisty = new TH2F(TString::Format(namescatter),
                             "Scatter plot of different TGenes", pop.size(), 0.,
                             50., pop.size(), 0., 50.);
          HYList.Add(myhisty);
        }
        myhisty->GetXaxis()->SetTitle(y1str);
        myhisty->GetYaxis()->SetTitle(y2str);
        myhisty->Fill(y1, y2);
      }
      // Distribution plots
      for (int j = 0; j < pop.GetTGenes(0).size(); ++j) {
        auto ind = pop.GetGeneValue(i, j);
        auto fitness = pop.GetObjectiveValue(i, j);
        PopDist->Fill(ind);
        PopFitnessDist->Fill(fitness);
      }
      // for TGraph2D - ONLY DTLZ1 - 3 objectives
      auto x = pop.GetObjectiveValue(i, 1);
      auto y = pop.GetObjectiveValue(i, 2);
      auto z = pop.GetObjectiveValue(i, 3);
      function[0] = x;
      function[1] = y;
      function[2] = z;
      // predictor DTZ1
      predictor[i] = x + y + z - 0.5;
      // add the 3d-data coordinate, the predictor value  and its errors
      data.Add(function, predictor[i], error);
    }
    PopDist->Draw();
    PopDist->Write();
    ////////////////////
    HXList.Write();
    HYList.Write();
    ////////////////////
    PopFitnessDist->Draw();
    PopFitnessDist->Write();
    ////////////////////
    bool ret = fitter.Fit(data);
    // if (ret) {
    const ROOT::Fit::FitResult &res = fitter.Result();
    // print result (should be around 1)
    res.Print(std::cout);
    // copy all fit result info (values, chi2, etc..) in TF3
    FitLand->SetFitResult(res);
    // test fit p-value (chi2 probability)
    double prob = res.Prob();
    FitLand->Draw();
    FitLand->Write();
    if (prob < 1.E-2)
      Error("HistoFill", "Bad data fit - fit p-value is %f", prob);
    //} else {
    // Error("HistoFill", "3D fit failed");
    //}
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
