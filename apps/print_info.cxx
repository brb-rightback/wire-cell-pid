#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  TString filename = argv[1];
  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("Event/Sim");
  int num = T->GetEntries();
  Float_t elifetime;
  T->SetBranchAddress("elifetime",&elifetime);
  Int_t runNo;
  T->SetBranchAddress("runNo",&runNo);
  if (num >0){
    T->GetEntry(0);
    std::cout << num << " " << runNo << " " << elifetime << std::endl;
  }
}
