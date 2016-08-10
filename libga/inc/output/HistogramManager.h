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

#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include "Rtypes.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TDirectory.h"
#include "TDirectory.h"
#include "TError.h"
#include "TF2.h"
#include "TF3.h"
#include "TFile.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TImage.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMultiGraph.h"
#include "TObject.h"
#include "TObject.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom2.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TVirtualFitter.h"
#include "generic/Functions.h"
#include "generic/Population.h"

#include "TGraphQQ.h"

#include <algorithm>
#include <iostream>
#include <iostream>
#include <string>
#include <vector>

#include <math.h>
#include <stdio.h>

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

  bool HistoFill(Population<F> &pop, char const *hfile, int generation) {
    std::cout << "Building output statistics for generation " << generation
              << std::endl;
    char namepop[20], namefitn[20], /* namefile[10],*/ namefolder[20],
        /* namefolderprevious[20],*/ namescatter[20], x1str[10], x2str[10],
        y1str[10], y2str[10], nameFitLand[20], name3dhist[20], name3dhistx[20],
        name3dhisty[20]/*, name3dhistyprevious[20]*/;
    std::vector<int> ScatterCombinationX, ScatterCombinationY;
    TObjArray HXList(0);
    TObjArray HYList(0);
    TRandom2 random;
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
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
        new TH1F(namepop, "Population distribution", pop.size(), 0, 1);
    PopDist->GetXaxis()->SetTitle("TGenes / bins");
    PopDist->GetYaxis()->SetTitle("N");
    PopDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopDist->SetMarkerStyle(20);
    PopDist->SetMarkerColor(kBlue);
    PopDist->SetMarkerSize(.6); //
                                /////// Population 3D Histogram
    TH3F *h3a =
        new TH3F(name3dhist, "3D Population", 20, -10, 30, 20, -10, 30, 20, -10, 30);
    // TH3F *h3a = new TH3F(name3dhist, "3D Population", 20, -30, 5, 20, -30, 5,
    //                     20, -30, 5);
    h3a->SetFillColor(kYellow); // Fill fill color to yellow
    h3a->SetFillColor(kYellow); // Fill fill color to yellow
    h3a->SetMarkerStyle(20);
    h3a->SetMarkerColor(kBlue);
    h3a->SetMarkerSize(.6); //
                            ////////////////////////////////
                            /////// Population 3 DHistogram
    TH3F *h3y = new TH3F(name3dhisty, "Y1/Y2/Y3", 20, -20, 30, 20, -10, 30, 20, -10, 30);
    // TH3F *h3y =
    //    new TH3F(name3dhisty, "Y1/Y2/Y3", 20, -30, 5, 20, -30, 5, 20, -30, 5);
    h3y->SetFillColor(kYellow); // Fill fill color to yellow
    h3y->SetFillColor(kYellow); // Fill fill color to yellow
    h3y->SetMarkerStyle(20);
    h3y->SetMarkerColor(kBlue);
    h3y->SetMarkerSize(.6); //
                            ////////////////////////////////
                            /////// Population 3 DHistogram
    TH3F *h3x = new TH3F(name3dhistx, "X1/X2/X3", 20, 0, 1, 20, 0, 1, 20, 0, 1);
    // TH3F *h3x =
    //    new TH3F(name3dhistx, "X1/X2/X3", 20, -5, 5, 20, -5, 5, 20, -5, 5);
    h3x->SetFillColor(kYellow); // Fill fill color to yellow
    h3x->SetFillColor(kYellow); // Fill fill color to yellow
    h3x->SetMarkerStyle(20);
    h3x->SetMarkerColor(kBlue);
    h3x->SetMarkerSize(.6); //
    ROOT::Fit::BinData data(pop.size(), 3);
    double function[3];
    double predictor[pop.size()];
    double parameterspredictor[3];
    // Predictor values for all DTLZ
    parameterspredictor[0] = 1;
    parameterspredictor[1] = 1;
    parameterspredictor[2] = 1;
    double error = 0.1;
    ////////////////////////////////
    // DTLZ b.
    TF3 *FitLand = new TF3(nameFitLand, F::TruePF, -10, 20, -10, 20, -10, 20, 3);
    // TF3 *FitLand = new TF3(nameFitLand, F::TruePF, -30, 5, -30, 5, -30, 5,
    // 3);

    FitLand->SetParameters(1, 1, 1);
    ROOT::Fit::Fitter fitter;
    // wrapped the TF1 in a IParamMultiFunction interface for teh Fitter class
    ROOT::Math::WrappedMultiTF1 wrapperfunction(*FitLand, 3);
    fitter.SetFunction(wrapperfunction);
    ////// Preparation for scatter combination
    ScatterCombinationX =
        ScatterCombinationCalculator(pop.GetTGenes(0).size(), 2);
    ScatterCombinationY =
        ScatterCombinationCalculator(pop.GetTFitness(0).size(), 2);

    //////// Fitness distribution
    // DTLZ b.
    TH1F *PopFitnessDist = new TH1F(namefitn, "Population fitness distribution ",pop.size(), -10., 20.);
    // TH1F *PopFitnessDist = new TH1F(namefitn, "Population fitness
    // distribution",
    //                                pop.size(), -30., 5.);
    PopFitnessDist->GetXaxis()->SetTitle("TGenes / bins");
    PopFitnessDist->GetYaxis()->SetTitle("N");
    PopFitnessDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopFitnessDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopFitnessDist->SetMarkerStyle(20);
    PopFitnessDist->SetMarkerColor(kBlue);
    PopFitnessDist->SetMarkerSize(.6); //
    /////////////////////////////////////////////////////////////////
    for (std::size_t i = 0; i < pop.size(); ++i) {
      double genearray[pop.GetTGenes(0).size()];
      // X Scatter plots
      for (std::size_t it = 0; it < ScatterCombinationX.size(); it += 2) {
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
                             "Scatter plot of different TGenes", pop.size(), 0,
                             1, pop.size(), 0, 1);
          // myhistx = new TH2F(TString::Format(namescatter),
          //                   "Scatter plot of different TGenes", pop.size(),
          // -5,
          //                   5, pop.size(), -5, 5);
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
        ///////////////////////////////
        // Previous generation
        /*
        sprintf(namefolderprevious, "%s%d", "PopulationStatisticsGeneration",
                generation - 1);
        TDirectory *dir = file.GetDirectory(namefolderprevious);
        if (dir == 0) {
          std::cout << "Fatal error: folder " << namefolderprevious
                    << " does not exist." << std::endl;
          return 0;
        }
        dir->cd();
        TH3F *myhistxold = (TH3F *)gROOT->FindObject(name3dhistyprevious);
        file.cd(namefolder);
        */
      }
      // Y Scatter plots
      for (std::size_t it = 0; it < ScatterCombinationY.size(); it += 2) {
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
          // DTLZ b.
          myhisty = new TH2F(TString::Format(namescatter),
                             "Scatter plot of different TGenes", pop.size(), 0,
                             100., pop.size(), -10, 20.);
          // Kursawe b.
          // myhisty = new TH2F(TString::Format(namescatter),
          //                   "Scatter plot of different TGenes", pop.size(),
          //                   -30, 5., pop.size(), -30, 5.);
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
      for (std::size_t j = 0; j < pop.GetTGenes(0).size(); ++j) {
        auto ind = pop.GetGeneValue(i, j);
        auto fitness = pop.GetObjectiveValue(i, j);
        genearray[j] = ind;
        PopDist->Fill(ind);
        PopFitnessDist->Fill(fitness);
      }
      // for TGraph2D - ONLY DTLZ - 3 objectives
      //////////////////////////////////////////
      // DATA + PREDICTOR for FIT
      auto x = pop.GetObjectiveValue(i, 0);
      auto y = pop.GetObjectiveValue(i, 1);
      auto z = pop.GetObjectiveValue(i, 2);
      function[0] = x;
      function[1] = y;
      function[2] = z;
      h3a->Fill(x, y, z);
      h3y->Fill(x, y, z);
      auto X1 = pop.GetGeneValue(i, 0);
      auto X2 = pop.GetGeneValue(i, 1);
      auto X3 = pop.GetGeneValue(i, 2);
      h3x->Fill(X1, X2, X3);
      predictor[i] =
          F::TruePF(genearray, parameterspredictor) + random.Gaus(0, error);
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
    /*
    TVirtualFitter *fit = TVirtualFitter::GetFitter();
    fit->PrintResults(2, 0.);
    TMatrixD *covMatrix =
        new TMatrixD(pop.GetTGenes(0).size(), pop.GetTGenes(0).size(),
    fit->GetCovarianceMatrix());
    covMatrix->Print();
    */
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

// ClassDefT2(HistogramManager, F)

#endif
