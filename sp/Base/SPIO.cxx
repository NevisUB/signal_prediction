#ifndef __SPIO_CXX__
#define __SPIO_CXX__

#include "SPIO.h"
#include "SPErr.h"

namespace sp {

  SPIO::SPIO(): 
    _in_tree(nullptr),
    _unfold_file(nullptr),
    _mc_tree_name(""),
    _unfold_file_name("")
  {}

  void SPIO::add_mc_in_file(const std::string& filename) {
    _in_mc_file_v.push_back(filename);
  }
  
  void SPIO::set_mc_tree_name(const std::string& treename) {
    _mc_tree_name = treename;
  }

  void SPIO::set_unfold_in_file(const std::string& filename) {
    _unfold_file_name = filename;
  }
  
  bool SPIO::initialize() {

    if(_unfold_file_name.empty())  {
      std::cout << "Unfold filename not set" << std::endl;
      _unfold_file_name = "unfold_data.root";
      std::cout << "Set unfold filename " << _unfold_file_name << std::endl;
    }
    
    _unfold_file = TFile::Open(_unfold_file_name.c_str(),"UPDATE");

    if (!_in_mc_file_v.empty()) {
      if(_mc_tree_name.empty()) throw sperr("No treename specified");
      
      _in_tree = new TChain(_mc_tree_name.c_str());
      for(const auto& fname : _in_mc_file_v) {
	std::cout << "READ: " << fname << std::endl;
	_in_tree->Add(fname.c_str());
      }

      _in_n_entries = _in_tree->GetEntries();
      if(!_in_n_entries) throw sperr("Empty mc tree specified (NO ENTRIES)");

      TObjArray *branch_list = _in_tree->GetListOfBranches();
      if(!branch_list->GetEntries()) throw sperr("Empty mc tree specified (NO BRANCHES)");
      
      for(size_t branch_id=0;branch_id<branch_list->GetEntries(); branch_id++)
	_in_mc_branch_v.insert(branch_list->At(branch_id)->GetName());
      
    }

    return true;
  }

  void SPIO::dump_branches() {
    std::cout << "tree called " << _mc_tree_name << " @ chain" << _in_tree << " has branches:" << std::endl;;
    std::cout << "[" << std::endl;
    for(const auto& branch : _in_mc_branch_v) 
      std::cout << branch << "," << std::endl;
    std::cout << "]" << std::endl;
  }
  
  // bool SPIO::calculate_migration() {

    

  // }
  
}

#endif
