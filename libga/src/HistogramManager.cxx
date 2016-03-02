#include <vector>

#include "TGenes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"
#include "TDirectory.h"
#include "TObject.h"

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
HistogramManager::HistogramManager():hAllev(0),hBuffev(0),hThread(0), hPriority(0), hSteps(0) {
/*
  hAllev = new TH1F("hAllev", "Totall number of events", 100, 50, 50);
  hBuffev = new TH1F("hBuffev", "Buffered events", 100, 50, 50);
  hThread = new TH1F("hThread", "Number of threads", 100, -16, 16);
  hPriority = new TH1F("hPriority", "Priority", 100, -1, 1);
  hSteps = new TH1F("hSteps", "Number steps for learning", 1000, 0, 1000);
  hVector = new TH1F("hVector", "Vector size", 100, 0, 100);
*/
}
#else
HistogramManager::HistogramManager():hx(0) {
/*
  hx = new TH1F("hx", "DTLZx benchmark", 100, 50, 50);
*/
}
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

Bool_t HistogramManager::HistoFill(Population<Double_t> &pop, char *file) {
  TFile *fHisto = TFile::Open(file);
  TTreeReader reader("Population", fHisto);
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
    // Creating monitoring canvas
    for (auto it = pop.fPopulation.begin(); it != pop.fPopulation.end(); ++it) {
#ifdef ENABLE_GEANTV
    hAllev = new TH1F("hAllev", "Totall number of events", 100, 50, 50);
    hAllev->GetYaxis()->SetRangeUser(0, 1);
    hAllev->GetXaxis()->SetNdivisions(nthreads + 1, true);
    hAllev->GetYaxis()->SetNdivisions(10, true);
    hAllev->SetFillColor(kGreen);
    hAllev->SetLineColor(0);
    hAllev->SetStats(false);
    hAllev->Draw();

    hBuffev = new TH1F("hBuffev", "Buffered events", 100, 50, 50);
    hBuffev->SetFillColor(kRed);
    hBuffev->SetFillStyle(3001);
    hBuffev->SetLineColor(0);
    hBuffev->SetStats(false);
    //cmon->cd(++ipad);
    hBuffev->Draw();
#else

#endif

    for (auto it = pop.fPopulation.begin(); it != pop.fPopulation.end(); ++it) {
#ifdef ENABLE_GEANTV
      hAllev->Fill(it->GetAllev(*it));
      hBuffev->Fill(it->GetBuffev(*it));
      hThread->Fill(it->GetThread(*it));
      hPriority->Fill(it->GetPriority(*it));
      hSteps->Fill(it->GetSteps(*it));
      hVector->Fill(it->GetVector(*it));
      hMaxVector->Fill(it->GetMaxVector(*it));
#else
      hx->Fill(it->GetGene(0));
#endif
    }
  }
  fHisto->cd();
  fHisto->Write();
  fHisto->Close();
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

}
