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
 * @brief Implementation of Histogramme class for LibGA
 * prototype
 */
//
#pragma once

#include "generic/Population.h"
#include "generic/Functions.h"
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
#include "TH3.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "Rtypes.h"
#include "TF3.h"
#include "TError.h"
#include "TProfile.h"
#include "TImage.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include "TRandom2.h"

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
    char namepop[20], namefitn[20], namefile[10], namefolder[20],
        namescatter[20], x1str[10], x2str[10], y1str[10], y2str[10],
        nameFitLand[20], name3dhist[20], name3dhistx[20], name3dhisty[20];
    std::vector<int> ScatterCombinationX, ScatterCombinationY;
    // TDirectory *folder;
    TObjArray HXList(0);
    TObjArray HYList(0);
    TRandom2 random;
    sprintf(namepop, "%s%d", "PopDist", generation);
    sprintf(namefitn, "%s%d", "PopFitnessDist", generation);
    sprintf(namefolder, "%s%d", "PopulationStatisticsGeneration", generation);
    sprintf(nameFitLand, "%s%d", "FitLand", generation);
    sprintf(name3dhist, "%s%d", "h3a", generation);
    sprintf(name3dhistx, "%s%d", "h3x", generation);
    sprintf(name3dhisty, "%s%d", "h3y", generation);
    // sprintf(namefile, "%s%d%s", "StatisticsX", generation, ".png");
    TFile file(hfile, "update");
    file.mkdir(namefolder);
    file.cd(namefolder);
    //////// Population distribution
    TH1F *PopDist =
        new TH1F(namepop, "Population distribution", pop.size(), 0., 1.);
    PopDist->GetXaxis()->SetTitle("TGenes / bins");
    PopDist->GetYaxis()->SetTitle("N");
    PopDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopDist->SetMarkerStyle(20);
    PopDist->SetMarkerColor(kBlue);
    PopDist->SetMarkerSize(.6); //
    /////// Population 3 DHistogram
    TH3F *h3a = new TH3F(name3dhist, "3D Population", 20, 0, 100, 20, 0, 100,
                         20, 0, 100);
    /////// Fitness landscape
    // TF3 *FitLand = new TF3(nameFitLand, "[0] * x + [1] * y  + [3] * z - 0.5",
    // 0,
    //                      100, 0, 100, 0, 100, 3);
    h3a->SetFillColor(kYellow); // Fill fill color to yellow
    h3a->SetFillColor(kYellow); // Fill fill color to yellow
    h3a->SetMarkerStyle(20);
    h3a->SetMarkerColor(kBlue);
    h3a->SetMarkerSize(.6); //
    ////////////////////////////////
    /////// Population 3 DHistogram
    TH3F *h3y =
        new TH3F(name3dhistx, "Y1/Y2/Y3", 20, 0, 10, 20, 0, 10, 20, 0, 10);
    h3y->SetFillColor(kYellow); // Fill fill color to yellow
    h3y->SetFillColor(kYellow); // Fill fill color to yellow
    h3y->SetMarkerStyle(20);
    h3y->SetMarkerColor(kBlue);
    h3y->SetMarkerSize(.6); //
    ////////////////////////////////
    /////// Population 3 DHistogram
    TH3F *h3x =
        new TH3F(name3dhisty, "X1/X2/X3", 20, 0, 1, 20, 0, 1, 20, 0, 1);
    h3x->SetFillColor(kYellow); // Fill fill color to yellow
    h3x->SetFillColor(kYellow); // Fill fill color to yellow
    h3x->SetMarkerStyle(20);
    h3x->SetMarkerColor(kBlue);
    h3x->SetMarkerSize(.6); //
    ////////////////////////////////
    TF3 *FitLand = new TF3(nameFitLand, F::TruePF, 0, 50, 0, 50, 0, 50, 3);
    FitLand->SetParameters(1, 1, 1);
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
    PopFitnessDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopFitnessDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopFitnessDist->SetMarkerStyle(20);
    PopFitnessDist->SetMarkerColor(kBlue);
    PopFitnessDist->SetMarkerSize(.6); //
    //////// Canvas distribution
    // TCanvas *Xcanvas = new TCanvas("Xcanvas");
    ////////////////////////////////////////////////////////////////
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
          myhistx->SetFillColor(kYellow); // Fill fill color to yellow
          myhistx->SetFillColor(kYellow); // Fill fill color to yellow
          myhistx->SetMarkerStyle(20);
          myhistx->SetMarkerColor(kBlue);
          myhistx->SetMarkerSize(.6); //
          HXList.Add(myhistx);
        }
        myhistx->GetXaxis()->SetTitle(x1str);
        myhistx->GetYaxis()->SetTitle(x2str);
        myhistx->Fill(x1, x2);
        /*
        Xcanvas->Divide(2, 2);
        Xcanvas->cd(1);
        myhistx->DrawClone("Contz");
        Xcanvas->cd(2);
        myhistx->DrawClone("surf3");
        Xcanvas->cd(3);
        myhistx->ProfileX()->DrawClone();
        Xcanvas->cd(4);
        myhistx->ProfileY()->DrawClone();
        Xcanvas->Draw();
        */
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
                             10., pop.size(), 0., 10.);
          myhisty->SetFillColor(kYellow); // Fill fill color to yellow
          myhisty->SetFillColor(kYellow); // Fill fill color to yellow
          myhisty->SetMarkerStyle(20);
          myhisty->SetMarkerColor(kBlue);
          myhisty->SetMarkerSize(.6); //
          HYList.Add(myhisty);
        }
        myhisty->GetXaxis()->SetTitle(y1str);
        myhisty->GetYaxis()->SetTitle(y2str);
        myhisty->Fill(y1, y2);
        // HYList.Draw("surf3");
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
      h3a->Fill(x, y, z);
      h3y->Fill(x, y, z);
      ///Just to check...!!!
      auto X1 = pop.GetGeneValue(i, 1);
      auto X2 = pop.GetGeneValue(i, 2);
      auto X3 = pop.GetGeneValue(i, 3);
      h3x->Fill(X1, X2, X3);
      // predictor DTZ1
      predictor[i] = x + y + z - 0.5 + random.Gaus(0, error);
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
    h3a->Draw("surf3");
    h3a->Fit(FitLand);
    h3a->Write();
    ///////////////////
    h3x->Draw();
    h3x->Write();
    ///////////////////
    h3y->Draw();
    h3y->Write();
    ///////////
    bool ret = fitter.Fit(data);
    if (ret) {
      const ROOT::Fit::FitResult &res = fitter.Result();
      // print result (should be around 1)
      res.Print(std::cout);
      // copy all fit result info (values, chi2, etc..) in TF3
      FitLand->SetFitResult(res);
      // test fit p-value (chi2 probability)
      double prob = res.Prob();
      // FitLand->Draw();
      // FitLand->Write();
      if (prob < 1.E-2)
        Error("HistoFill", "Bad data fit - fit p-value is %f", prob);
      else
        std::cout << "Good fit : p-value  = " << prob << std::endl;
    } else
      Error("HistoFill", "3D fit failed");
    return true;
  }
  void Reset();

private:
  HistogramManager() {}

  ~HistogramManager() {}

  // ClassDef(HistogramManager, 1)
};

// ClassDefT2(HistogramManager,F)

#endif
