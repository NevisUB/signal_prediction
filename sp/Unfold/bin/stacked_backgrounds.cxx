#include "Combined/CombinedFunctions.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv) {

  std::string fname(argv[1]);

  auto tf = TFile::Open(fname.c_str(),"READ");

  int NFSP;
  std::vector<int>* FSPType_v = nullptr;

  std::vector<float>* VertexX_v = nullptr;
  std::vector<float>* VertexY_v = nullptr;
  std::vector<float>* VertexZ_v = nullptr;

  std::vector<float>* MomX_v = nullptr;
  std::vector<float>* MomY_v = nullptr;
  std::vector<float>* MomZ_v = nullptr;
  std::vector<float>* MomT_v = nullptr;
  
  int NUANCEChan;
  int NuType;
  int NuParentID;
  bool PassOsc;

  float RecoEnuQE;
  float Weight;

  auto tree = (TTree*) tf->Get("MiniBooNE_CCQE");
    
  std::cout << "Got " << tree->GetEntries() << " @ tree=" << tree << " called " << tree->GetName() << std::endl;
    
  tree->SetBranchAddress("NFSP"     , &NFSP);
  tree->SetBranchAddress("FSPType"  , &FSPType_v);
  tree->SetBranchAddress("VertexX"  , &VertexX_v);
  tree->SetBranchAddress("VertexY"  , &VertexY_v);
  tree->SetBranchAddress("VertexZ"  , &VertexZ_v);

  tree->SetBranchAddress("MomX"  , &MomX_v);
  tree->SetBranchAddress("MomY"  , &MomY_v);
  tree->SetBranchAddress("MomZ"  , &MomZ_v);
  tree->SetBranchAddress("MomT"  , &MomT_v);

  tree->SetBranchAddress("NUANCEChan" , &NUANCEChan);
  tree->SetBranchAddress("NuType"     , &NuType);
  tree->SetBranchAddress("NuParentID" , &NuParentID);
  tree->SetBranchAddress("PassOsc"    , &PassOsc);
  
  tree->SetBranchAddress("RecoEnuQE", &RecoEnuQE);
  tree->SetBranchAddress("Weight"   , &Weight);
  
  
  
  std::vector<TH1D
  
  for(size_t entry = 0; entry < (size_t)tree->GetEntries(); ++entry) {
    std::cout << "@ entry=" << entry << std::endl;
    tree->GetEntry(entry);
    
    if(!PassOsc) continue;

    //
    // check if gamma is pi0
    //

    unsigned isPi0 = sp::Pi0Details(NFSP,*FSPType_v,
				    *VertexX_v,*VertexY_v,*VertexZ_v,
				    *MomX_v,*MomY_v,*MomZ_v,*MomT_v);
    
    //
    // get the background type ID
    //
    sp::StackedBkgdType_t bkg_type = sp::StackHistoBkgd(0,(bool)isPi0,(sp::NuanceType_t)NUANCEChan,(sp::NuType_t)NuType,(sp::GEANT3Type_t)NuParentID);

    //
    // it not a background, move on
    //
    if (bkg_type == sp::kBKGD_INVALID) continue;
    
    //
    // it is a background
    //
    

  }

  return 0;
}
