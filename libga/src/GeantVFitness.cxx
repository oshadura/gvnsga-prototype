#include "GeantVFitness.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"

#include <algorithm>


ClassImp(GeantVFitness)

    struct CompairMemResident {

  bool operator()(ProcInfo_t lhs, ProcInfo_t rhs) {
    return lhs.fMemResident < rhs.fMemResident;
  }
};
struct CompairMemVirtual {
  bool operator()(ProcInfo_t lhs, ProcInfo_t rhs) {
    return lhs.fMemVirtual < rhs.fMemVirtual;
  }
};

void GeantVFitness::LogMemoryFitness() {
  ProcInfo_t info;
  gSystem->GetProcInfo(&info);
  fMemoryVector.push_back(info);
}

void GeantVFitness::LogTimeFitness() {
}

void GeantVFitness::HistOutputFitness(std::string file) {
  hfile = new TFile(file.c_str(),"RECREATE");
  hfile->mkdir("Fitness");
  hfile->cd("Fitness");
  int numBins = fMemoryVector.size();
  hMemRes = new TH1F ("memory_resident", "Resident memory usage", numBins, 0, numBins);
  hMemRes->GetXaxis()->SetTitle("Gene");
  hMemRes->GetYaxis()->SetTitle("Resident memory (GB)");
  hMemVirt = new TH1F("memory_virtual", "Virtual memory usage", numBins, 0, numBins);
  hMemVirt->GetXaxis()->SetTitle("Gene");
  hMemVirt->GetYaxis()->SetTitle("Virtual memory (GB)");
  int bin = 1;
  for (std::vector<ProcInfo_t>::iterator it = fMemoryVector.begin();
       it != fMemoryVector.end(); it++, bin++) {
    ProcInfo_t info = (*it);
    hMemRes->SetBinContent(bin, info.fMemResident / (1024. * 1024.));
    hMemVirt->SetBinContent(bin, info.fMemVirtual / (1024. * 1024.));
  }
  hMemRes->Write();
  hMemVirt->Write();
  double maxMemResident =
      (std::max_element(fMemoryVector.begin(), fMemoryVector.end(),
                        CompairMemResident()))
          ->fMemResident;
  double maxMemVirtual =
      (std::max_element(fMemoryVector.begin(), fMemoryVector.end(),
                        CompairMemVirtual()))
          ->fMemVirtual;
  std::printf("Maximum resident memory usage:%f\n",
              (maxMemResident / (1024. * 1024.)));
  std::printf("Maximum virtual memory usage:%f\n",
              (maxMemVirtual / (1024. * 1024.)));
  maxMemVirtual = 0;
  maxMemResident = 0;
  fMemoryVector.clean();
  hfile->cd();
  hfile->Write();
  hfile->Close();
}
