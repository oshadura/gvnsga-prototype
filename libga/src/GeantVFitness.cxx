#include "GeantVFitness.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"

#include <algorithm>
#include <iostream>

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

void GeantVFitness::TemporarySolution() {
  ProcInfo_t info;
  gSystem->GetProcInfo(&info);
  fMemoryVector.push_back(info);
}

Double_t GeantVFitness::LogMemoryFitness(std::string file) {
  if (!hfile) {
    hfile = new TFile(file.c_str(), "RECREATE");
    hfile->mkdir("Fitness");
    hfile->cd("Fitness");
  } else {
    hfile = new TFile(file.c_str(), "UPDATE");
  }
  if (fMemorySwitch = true) {
    ProcInfo_t info;
    gSystem->GetProcInfo(&info);
    fMemoryVector.push_back(info);
    sleep(10);
  }
  for (auto &i : fMemoryVector) {
    std::cout << "DEBUG: " << i.fMemResident / (1024. * 1024.) << " GB"
              << std::endl;
  }
  int numBins = fMemoryVector.size();
  if (!hMemRes && !hMemRes) {
    hMemRes = new TH1F("memory_resident", "Resident memory usage", numBins, 0,
                       numBins);
    hMemVirt =
        new TH1F("memory_virtual", "Virtual memory usage", numBins, 0, numBins);
    hMemRes->GetXaxis()->SetTitle("Gene");
    hMemRes->GetYaxis()->SetTitle("Resident memory (GB)");
    hMemVirt->GetXaxis()->SetTitle("Gene");
    hMemVirt->GetYaxis()->SetTitle("Virtual memory (GB)");
  }
  int bin = 1;
  for (std::vector<ProcInfo_t>::iterator it = fMemoryVector.begin();
       it != fMemoryVector.end(); it++, bin++) {
    ProcInfo_t info = (*it);
    hMemRes->SetBinContent(bin, info.fMemResident / (1024. * 1024.));
    hMemVirt->SetBinContent(bin, info.fMemVirtual / (1024. * 1024.));
  }
  hMemVirt->SetLineColor(kBlue);
  hMemVirt->SetLineStyle(kDashed);
  hMemVirt->SetLineColor(kRed);
  hMemVirt->SetLineStyle(kDashed);
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
  return maxMemResident;
}

void GeantVFitness::LogTimeFitness() {}
