#include <iostream>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char* argv[])
{
  TFile *file = new TFile(argv[1]);
  TTree *T_proj_data = (TTree*)file->Get("T_proj_data");
  if (T_proj_data==0)
    std::cout << argv[1] << " " << T_proj_data << std::endl;
}
