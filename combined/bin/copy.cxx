#include <iostream>
#include <exception>
#include <string>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char** argv) {
  
  std::string filename_in(argv[1]);
  std::cout << "Got " << filename_in << " input filename" << std::endl;
  std::string filename_out = "filtered_" + filename_in;
  
  if (argc!=2) throw std::exception();

  TFile *file_in  = new TFile(filename_in.c_str() ,"READ");
  file_in->cd();
  
  TTree *tree_in  = (TTree*) file_in->Get("MiniBooNE_CCQE");
  std::cout << "Read in " << tree_in->GetEntries() << std::endl;

  TFile *file_out = new TFile(filename_out.c_str(),"RECREATE");
  file_out->cd();
  
  TTree *tree_out = tree_in->CopyTree("Weight>0");
  std::cout << "Writing " << tree_out->GetEntries() << std::endl;
  
  tree_out->Write();

  file_in->Close();
  file_out->Close();

  std::cout << "Wrote " << filename_out << std::endl;
  return 0;
}
 

