#include <vector>

#include "Genes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"
#include "TDirectory.h"

HistogramManager::HistogramManager(TDirectory *dir) {
  TDirectory *saved = gDirectory;
  dir->cd();
  fPopulationSize =
      new TH1F("hNIndividual", "Number of individuals", 100, 500, 500);
  fAllev = new TH1F("hAllev", "Totall number of events", 100, 50, 50);
  fBuffev = new TH1F("hBuffev", "Buffered events", 100, 50, 50);
  fThread = new TH1F("hThread", "Number of threads", 100, -16, 16);
  fPriority = new TH1F("hPriority", "Priority", 100, -1, 1);
  fSteps = new TH1F("hSteps", "Number steps for learning", 1000, 0, 1000);
  fVector = new TH1F("hVector", "Vector size", 100, 0, 100);
  saved->cd();
}

/*bool HistogramManager::CheckValue(ROOT::Internal::TTreeReaderValueBase *value)
{
  if (value->GetSetupStatus() < 0) {
    std::cerr << "Error " << value->GetSetupStatus() << "setting up reader for "
              << value->GetBranchName() << '\n';
    return false;
  }
  return true;
}

*/
HistogramManager::~HistogramManager() {}

void HistogramManager::HFill(Population *pop, char *file) {
  /*
  TFile* ga = TFile::Open(file);
  TTreeReader reader("Genetic Algorithm TTree", ga);
  TTreeReaderValue<std::vector<Genes>> Population(reader, "Population");
  TTreeReaderArray<Genes> Genes(reader, "Genes");
  if (!CheckValue(Population)) return false;
  if (!CheckValue(Genes)) return false;

  while (reader.Next()) {
    if (reader.GetEntryStatus() == kEntryValid) {
      std::cout << "Loaded entry " << reader.GetCurrentEntry() << '\n';
    } else {
      switch (reader.GetEntryStatus()) {
        kEntryValid:
        // Handled above.
          break;
        kEntryNotLoaded:
          std::cerr << "Error: TTreeReader has not loaded any data yet!\n";
          break;
        kEntryNoTree:
          std::cerr << "Error: TTreeReader cannot find a tree names \"Genetic
Algorithm TTree\"!\n";
          break;
        kEntryNotFound:
          // Can't really happen as TTreeReader::Next() knows when to stop.
          std::cerr << "Error: The entry number doe not exist\n";
          break;
        kEntryChainSetupError:
          std::cerr << "Error: TTreeReader cannot access a chain element, e.g.
file without the tree\n";
          break;
        kEntryChainFileError:
          std::cerr << "Error: TTreeReader cannot open a chain element, e.g.
missing file\n";
          break;
        kEntryDictionaryError:
          std::cerr << "Error: TTreeReader cannot find the dictionary for some
data\n";
           break;
         }
      return false;
    }
if (!Genes.IsEmpty()) {
         float currentWeight = *weight;
         for (int iGene  = 0, nGene = Genes.GetSize(); iGene < nGene; ++iGene) {
            hist->Fill([iGene].fAllev(), hAllev);
            hist->Fill([iGene].fBuffev(), hBuffev);
            hist->Fill([iGene].fThread(), hThread);
            hist->Fill([iGene].fPriority(), hPriority);
            hist->Fill([iGene].fSteps(), hSteps);
            hist->Fill([iGene].fVector(), hVector);
         }
      }
      */
}

ClassImp(HistogramManager)
