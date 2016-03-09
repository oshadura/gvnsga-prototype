#include <vector>

#include "TGenes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

ClassImp(HistogramManager)

    HistogramManager *HistogramManager::HistoInstance = 0;

HistogramManager *HistogramManager::Instance() {
  if (HistoInstance == 0)
    HistoInstance = new HistogramManager();
  return HistoInstance;
}

#ifdef ENABLE_GEANTV
HistogramManager::HistogramManager()
    : hAllev(0), hBuffev(0), hThread(0), hPriority(0), hSteps(0) {}
#else
HistogramManager::HistogramManager() : hx(0) {}
#endif

HistogramManager::~HistogramManager() { HistoInstance = 0; }

Bool_t
HistogramManager::CheckValue(ROOT::Internal::TTreeReaderValueBase *value) {
  if (value->GetSetupStatus() < 0) {
    std::cerr << "Error " << value->GetSetupStatus() << "setting up reader for "
              << value->GetBranchName() << '\n';
    return false;
  }
  return true;
}

Bool_t HistogramManager::HistoFill(Population<Double_t> &pop, char *hfile) {
  std::cout << "Filling histograms.." << std::endl;
  TFile *file = TFile::Open(hfile);
  TTreeReader reader("Population", file);
  TTreeReaderArray<std::vector<std::vector<Double_t>>> Population(reader,
                                                                  "Population");
  if (!CheckValue(&Population)) {
    return false;
  }
  while (reader.Next()) {
    if (reader.GetEntryStatus() == TTreeReader::kEntryValid) {
      std::cout << "Loaded entry " << reader.GetCurrentEntry() << '\n';
    } else {
      switch (reader.GetEntryStatus()) {
      case TTreeReader::kEntryValid:
        std::cerr << "Errorrrrrrrrrrrrrr!\n";
        break;
      case TTreeReader::kEntryNotLoaded:
        std::cerr << "Error: TTreeReader has not loaded any data yet!\n";
        break;
      case TTreeReader::kEntryNoTree:
        std::cerr << "Error: TTreeReader cannot find a tree names -> Genetic "
                     "Algorithm TTree!\n";
        break;
      case TTreeReader::kEntryNotFound:
        std::cerr << "Error: The entry number doe not exist\n";
        break;
      case TTreeReader::kEntryChainSetupError:
        std::cerr << "Error: TTreeReader cannot access a chain element, "
                     "e.g.file without the tree\n";
        break;
      case TTreeReader::kEntryChainFileError:
        std::cerr << "Error: TTreeReader cannot open a chain element, "
                     "e.g.missing file\n";
        break;
      case TTreeReader::kEntryDictionaryError:
        std::cerr
            << "Error: TTreeReader cannot find the dictionary for some data\n";
        break;
      }
      return false;
    }
    TCanvas *mon =
        (TCanvas *)gROOT->GetListOfCanvases()->FindObject("Monitoring");
    if (pop.GetNParam() == 1)
      mon->Divide(1, 1);
    else
      mon->Divide(2, pop.GetNParam());
    int ipad = 0;
// Creating monitoring canvas
#ifdef ENABLE_GEANTV
    hAllev = new TH1F("hAllev", "Totall number of events", 100, 0, 100);
    hAllev->SetLineColor(kMagenta);
    hAllev->SetStats(false);
    hAllev->Fill();
    mon->cd(++ipad);
    hAllev->Draw("SAME");
    // Check also SCAT
    //////////////////////////
    hBuffev = new TH1F("hBuffev", "Buffered events", 100, 0, 99);
    // hBuffev->SetFillColor(kRed);
    // hBuffev->SetFillStyle(3001);
    hBuffev->SetLineColor(0);
    hBuffev->SetStats(false);
    hBuffev->Fill();
    mon->cd(++ipad);
    hBuffev->Draw("SAME");
    /////////////////////////
    hThread = new TH1F("hThread", "Number of threads", 100, 0, 16);
    // hThread->SetFillColor(kRed);
    // hThread->SetFillStyle(3001);
    hThread->SetLineColor(0);
    hThread->SetStats(false);
    hThread->FIll();
    mon->cd(++ipad);
    hThread->Draw("SAME");
    //////////////////////////
    hPriority = new TH1F("hPriority", "Priority", 100, 0, 1);
    // hPriority->SetFillColor(kRed);
    // hPriority->SetFillStyle(3001);
    hPriority->SetLineColor(0);
    hPriority->SetStats(false);
    hPriority->Fill();
    mon->cd(++ipad);
    hPriority->Draw("SAME");
    //////////////////////////
    hSteps = new TH1F("hSteps", "Number steps for learning", 10000, 0, 10000);
    // hSteps->SetFillColor(kRed);
    // hSteps->SetFillStyle(3001);
    hSteps->SetLineColor(0);
    hSteps->SetStats(false);
    hSteps->Fill();
    mon->cd(++ipad);
    hSteps->Draw("SAME");
    //////////////////////////
    hVector = new TH1F("hVector", "Vector size", 64, 0, 64);
    // hVector->SetFillColor(kRed);
    // hVector->SetFillStyle(3001);
    hVector->SetLineColor(0);
    hVector->SetStats(false);
    hVector->Fill();
    mon->cd(++ipad);
    hVector->Draw("SAME");
    //////////////////////////
    hMaxVector = new TH1F("hMaxVector", "Max Vector size", 512, 0, 512);
    // hMaxVector->SetFillColor(kRed);
    // hMaxVector->SetFillStyle(3001);
    hMaxVector->SetLineColor(0);
    hMaxVector->SetStats(false);
    hMaxVector->Fill();
    mon->cd(++ipad);
    hMaxVector->Draw("SAME");
#else
    hx = new TH1F("hx", "DTLZx benchmark", 100, 0, 1);
    // hx->SetFillColor(kRed);
    // hx->SetFillStyle(3001);
    hx->SetLineColor(0);
    hx->SetStats(false);
    hx->Fill();
    mon->cd(++ipad);
    hx->Draw("SAME");
#endif
    mon->Update();
    double stamp = 0.;
    int i, j, bin;
    int nmaxtot;
    while (1) { // exit condition here
      i = int(stamp);
      ipad = 0;
      gSystem->Sleep(50); // millisec
      // Fill histograms
      for (auto it = pop.fPopulation.begin(); it != pop.fPopulation.end();
           ++it) {
        if (stamp > 100) {
#ifdef ENABLE_GEANTV
          if (hAllev) {
            hAllev->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hAllev->SetBinContent(j + 1, it->GetAllev(*it));
          }
          if (hBuffev) {
            hBuffev->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hBuffev->SetBinContent(j + 1, it->GetBuffev(*it));
          }
          if (hThread) {
            hThread->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hThread->SetBinContent(j + 1, it->GetThread(*it));
          }
          if (hPriority) {
            hPriority->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hAllev->SetBinContent(j + 1, it->GetPriority(*it));
          }
          if (hSteps) {
            hSteps->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hSteps->SetBinContent(j + 1, it->GetSteps(*it));
          }
          if (hVector) {
            hVector->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hAllev->SetBinContent(j + 1, it->GetVector(*it));
          }
          if (hMaxVector) {
            hMaxVector->GetXaxis()->Set(100, stamp - 100, stamp);
            for (j = 0; j < 100; j++)
              hAllev->SetBinContent(j + 1, it->GetMaxVector(*it));
          }
#else
          if (hx) {
            hx->SetBinContent(i + 1, it->GetGene(0));
          }
#endif
        } else {
#ifdef ENABLE_GEANTV
          if (hAllev) {
            hAllev->SetBinContent(i + 1, it->GetAllev(*it));
          }
          if (hBuffev) {
            hBuffev->SetBinContent(i + 1, it->GetBuffev(*it));
          }
          if (hThread) {
            hThread->SetBinContent(i + 1, it->GetThread(*it));
          }
          if (hPriority) {
            hPriority->SetBinContent(i + 1, it->GetPriority(*it));
          }
          if (hSteps) {
            hSteps->SetBinContent(i + 1, it->GetSteps(*it));
          }
          if (hVector) {
            hVector->SetBinContent(i + 1, it->GetVector(*it));
          }
          if (hMaxVector) {
            hMaxVector->SetBinContent(i + 1, it->GetMaxVector(*it));
          }
#else
          if (hx) {
            hx->SetBinContent(i + 1, it->GetGene(0));
          }
#endif
        }
#ifdef ENABLE_GEANTV
        if (hAllev) {
          mon->cd(++ipad);
          hAllev->Draw();
        }
        if (hBuffev) {
          mon->cd(++ipad);
          hBuffev->Draw();
        }
        if (hThread) {
          mon->cd(++ipad);
          hThread->Draw();
        }
        if (hPriority) {
          mon->cd(++ipad);
          hThread->Draw();
        }
        if (hSteps) {
          mon->cd(++ipad);
          hSteps->Draw();
        }
        if (hVector) {
          mon->cd(++ipad);
          hVector->Draw();
        }
        if (hMaxVector) {
          mon->cd(++ipad);
          hMaxVector->Draw();
        }
#else
        if (hx) {
          mon->cd(++ipad);
          hx->Draw();
        }
#endif
      }
      mon->Modified();
      mon->Update();
      stamp += 1;
    }
  }
  gSystem->Sleep(100); // millisec
  return true;
}

void HistogramManager::Reset() {
#ifdef ENABLE_GEANTV
  hAllev = 0;
  hBuffev = 0;
  hThread = 0;
  hPriority = 0;
  hSteps = 0;
  hVector = 0;
  hMaxVector = 0;
#else
  hx = 0;
#endif
  HistoInstance = 0;
}
