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
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TKey.h"
#include "TFile.h"
#include "THStack.h"
#include "TMatrixD.h"

#include "TGraphQQ.h"

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

  // get a TCanvas with a QQ plot with confidence band
  TCanvas *GetQQPlotCanvas(TH3F *h1, TH3F *h2) {

    // Quantiles: 20
    const int nq = 20;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    Double_t yq1[nq]; // array to contain the quantiles
    Double_t yq2[nq]; // array to contain the quantiles

    for (Int_t i = 0; i < nq; i++)
      xq[i] = Float_t(i + 1) / nq;
    h1->GetQuantiles(nq, yq1, xq);
    h2->GetQuantiles(nq, yq2, xq);

    Double_t xq_plus[nq];
    Double_t xq_minus[nq];
    Double_t yq2_plus[nq];
    Double_t yq2_minus[nq];

    /* KS_cv: KS critical value

               1.36
    KS_cv = -----------
             sqrt( N )

    Where 1.36 is for alpha = 0.05 (confidence level 1-5%=95%, about 2 sigma)

    For 1 sigma (alpha=0.32, CL=68%), the value in the nominator is 0.9561, it
    is gotten by
    GetCriticalValue(1, 1 - 0.68).

    NOTE:
    o For 1-sample KS test (data and theoretic), N should be n
    o For 2-sample KS test (2 data set), N should be sqrt(m*n/(m+n))! Here is
    the case
      m or n (size of samples) should be effective size for a histogram
    o Critival value here is valid for only for sample size >= 80 (some
    references say 35)
      which means, for example, for a unweighted histogram, it must have more
    than 80 (or 35)
      entries filled and then confidence band is reliable.

    */

    float esum1 = GetEffectiveSampleSize(h1);
    float esum2 = GetEffectiveSampleSize(h2);
    std::cout << esum1 << std::endl;

    // one sigma band
    float KS_cv = GetCriticalValue(1, 1 - 0.68) /
                  TMath::Sqrt((esum1 * esum2) / (esum1 + esum2));

    for (Int_t i = 0; i < nq; i++)
      xq_plus[i] = (Float_t)(xq[i] + KS_cv); // upper limit
    for (Int_t i = 0; i < nq; i++)
      xq_minus[i] = (Float_t)(xq[i] - KS_cv); // lower limit

    h2->GetQuantiles(nq, yq2_plus, xq_plus);
    h2->GetQuantiles(nq, yq2_minus, xq_minus);

    double yq2_err_plus[nq];
    double yq2_err_minus[nq];
    for (Int_t i = 0; i < nq; i++) {
      yq2_err_plus[i] = yq2_plus[i] - yq2[i];
      yq2_err_minus[i] = yq2[i] - yq2_minus[i];
    }

    TCanvas *c = new TCanvas("c", "QQ with CL", 600, 450);
    TGraph *gr = new TGraph(
        nq - 1, yq1, yq2); // forget the last point, so number of points: (nq-1)
    gr->SetLineColor(kRed + 2);
    gr->SetMarkerColor(kRed + 2);
    gr->SetMarkerStyle(20);
    gr->SetTitle(Form("QQ plot; %s Quantile; %s Quantile", h1->GetName(),
                      h2->GetName()));
    gr->Draw("ap");
    double x_min = gr->GetXaxis()->GetXmin();
    double x_max = gr->GetXaxis()->GetXmax();
    double y_min = gr->GetXaxis()->GetXmin();
    double y_max = gr->GetXaxis()->GetXmax();
    c->Clear();

    // some debug codes:
    //   printf("x_min: %f\n", (float)x_min);
    //   printf("x_max: %f\n", (float)x_max);
    //   printf("y_min: %f\n", (float)y_min);
    //   printf("y_max: %f\n", (float)y_max);

    // add confidence level band in gray
    TGraphAsymmErrors *ge = new TGraphAsymmErrors(nq - 1, yq1, yq2, 0, 0,
                                                  yq2_err_minus, yq2_err_plus);
    ge->SetFillColor(17);

    ///////////////////
    // put all together
    TMultiGraph *mg = new TMultiGraph("mg", "");
    mg->SetMinimum(y_min);
    mg->SetMaximum(y_max);
    mg->Add(gr, "ap");
    mg->Add(ge, "3");
    mg->Add(gr, "p");
    mg->Draw();

    // a straight line y=x to be a reference
    TF1 *f_dia = new TF1("f_dia", "x", h1->GetXaxis()->GetXmin(),
                         h1->GetXaxis()->GetXmax());
    f_dia->SetLineColor(9);
    f_dia->SetLineWidth(2);
    f_dia->SetLineStyle(2);
    f_dia->Draw("same");

    TLegend *leg = new TLegend(0.52, 0.15, 0.87, 0.35);
    leg->SetFillColor(0);
    leg->SetShadowColor(17);
    leg->SetBorderSize(3);
    leg->AddEntry(gr, "QQ points", "p");
    leg->AddEntry(ge, "68% CL band", "f");
    leg->AddEntry(f_dia, "Diagonal line", "l");
    leg->Draw();

    return c;
  }

  // calculate effective sample size for a histogram
  // same way as ROOT does.
  float GetEffectiveSampleSize(TH3F *h) {
    TAxis *axis = h->GetXaxis();
    Int_t last = axis->GetNbins();
    float esum = 0;
    float sum = 0, ew = 0, w = 0;
    for (int bin = 1; bin <= last; bin++) {
      sum += h->GetBinContent(bin);
      ew = h->GetBinError(bin);
      w += ew * ew;
    }
    esum = sum * sum / w;
    return esum;
  }

  //====================================================
  // the routine is used to calculate critical value given
  // n and p, confidential level = 1 - p.
  // Original Reference:
  // http://velveeta.che.wisc.edu/octave/lists/archive//octave-sources.2003/msg00031.html
  // I just checked it, but it is not available now...
  double GetCriticalValue(int n, double p) {
    double dn = 1;
    double delta = 0.5;
    double res;
    res = TMath::KolmogorovProb(dn * sqrt(n));
    while (res > 1.0001 * p || res < 0.9999 * p) {
      if (res > 1.0001 * p)
        dn = dn + delta;
      if (res < 0.9999 * p)
        dn = dn - delta;
      delta = delta / 2.;
      res = TMath::KolmogorovProb(dn * sqrt(n));
    }
    return dn;
  }
  //========================================================

  bool HistoFill(Population<F> &pop, char *hfile, int generation) {
    std::cout << "Building output statistics for generation " << generation
              << std::endl;
    char namepop[20], namefitn[20], namefile[10], namefolder[20],
        namefolderprevious[20], namescatter[20], x1str[10], x2str[10],
        y1str[10], y2str[10], nameFitLand[20], name3dhist[20], name3dhistx[20],
        name3dhisty[20], name3dhistyprevious[20];
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
        new TH1F(namepop, "Population distribution", pop.size(), 0., 1.);
    PopDist->GetXaxis()->SetTitle("TGenes / bins");
    PopDist->GetYaxis()->SetTitle("N");
    PopDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopDist->SetFillColor(kYellow); // Fill fill color to yellow
    PopDist->SetMarkerStyle(20);
    PopDist->SetMarkerColor(kBlue);
    PopDist->SetMarkerSize(.6); //
    /////// Population 3 DHistogram
    TH3F *h3a =
        new TH3F(name3dhist, "3D Population", 20, 0, 1, 20, 0, 1, 20, 0, 1);
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
    TH3F *h3y = new TH3F(name3dhistx, "Y1/Y2/Y3", 20, 0, 5, 20, 0, 5, 20, 0, 5);
    h3y->SetFillColor(kYellow); // Fill fill color to yellow
    h3y->SetFillColor(kYellow); // Fill fill color to yellow
    h3y->SetMarkerStyle(20);
    h3y->SetMarkerColor(kBlue);
    h3y->SetMarkerSize(.6); //
    ////////////////////////////////
    /////// Population 3 DHistogram
    TH3F *h3x = new TH3F(name3dhisty, "X1/X2/X3", 20, 0, 1, 20, 0, 1, 20, 0, 1);
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
        ///////////////////////////////
        // Previous generation
        sprintf(namefolderprevious, "%s%d", "PopulationStatisticsGeneration",
                generation - 1);
        sprintf(name3dhistyprevious, "%s%d", "h3y", generation - 1);
        // TString rootname = namefolderprevious + "/" + name3dhistyprevious;
        TDirectory *dir = file.GetDirectory(namefolderprevious);
        if (dir == 0) {
          std::cout << "Fatal error: folder " << namefolderprevious
                    << " does not exist." << std::endl;
          return 0;
        }
        dir->cd();
        TH3F *myhistxold = (TH3F *)gROOT->FindObject(name3dhistyprevious);
        THStack *hs =
            new THStack("hs", "2 distributions POpulatiion N-1 & Population N");
        hs->Add(h3y);
        hs->Add(myhistxold);
        hs->Draw("nostack");
        // draw legend
        TLegend *leg = new TLegend(1, 1, 1, 1);
        leg->SetFillColor(0);
        leg->AddEntry(h3y, h3y->GetTitle(), "pl");
        leg->AddEntry(myhistxold, myhistxold->GetTitle(), "pl");
        leg->Draw();

        TCanvas *can_qq = GetQQPlotCanvas(h3y, myhistxold);
        can_qq->Draw();
        file.cd(namefolder);
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
      /// Just to check...!!!
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
    TVirtualFitter *fit = TVirtualFitter::GetFitter();
    fit->PrintResults(2, 0.);
    TMatrixD *covMatrix =
        new TMatrixD(pop.GetTGenes(0).size(), pop.GetTGenes(0).size(), fit->GetCovarianceMatrix());
    covMatrix->Print();
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
