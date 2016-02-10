#include <vector>

#include "TGenes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"
#include "TDirectory.h"

ClassImp(HistogramManager)

    HistogramManager::HistogramManager(TDirectory *dir) {
  TDirectory *saved = gDirectory;
  dir->cd();
  TH1F hNIndividual("hNIndividual", "Number of individuals", 100, 500, 500);
  TH1F hAllev("hAllev", "Totall number of events", 100, 50, 50);
  TH1F hBuffev("hBuffev", "Buffered events", 100, 50, 50);
  TH1F hThread("hThread", "Number of threads", 100, -16, 16);
  TH1F hPriority("hPriority", "Priority", 100, -1, 1);
  TH1F hSteps("hSteps", "Number steps for learning", 1000, 0, 1000);
  TH1F hVector("hVector", "Vector size", 100, 0, 100);
  saved->cd();
}

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

  TFile *ga = TFile::Open(file);
  TTreeReader reader("Population", ga);
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
    for (auto it = pop.fPopulation.begin(); it != pop.fPopulation.end(); ++it) {
      hAllev->Fill(it->GetAllev(*it));
      hBuffev->Fill(it->GetBuffev(*it));
      hThread->Fill(it->GetThread(*it));
      hPriority->Fill(it->GetPriority(*it));
      hSteps->Fill(it->GetSteps(*it));
      hVector->Fill(it->GetVector(*it));
      hMaxVector->Fill(it->GetMaxVector(*it));
    }
  }
}
