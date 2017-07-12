#ifndef __SPIO_CXX__
#define __SPIO_CXX__

#include "SPIO.h"
#include "SPErr.h"
#include <sstream>

namespace sp {

  SPIO::SPIO() : 
    _mc_tree_name(""),
    _unfold_file_name(""),
    _in_tree(nullptr),
    _unfold_file(nullptr)
  {std::cout << "S @ " << this << std::endl;}

  SPIO::~SPIO() {
    std::cout << "~S @ " << this << std::endl;
  }
  
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
    
    TH1::AddDirectory(kFALSE);

    if(_unfold_file_name.empty())  {
      std::cout << "Unfold filename not set" << std::endl;
      _unfold_file_name = "unfold_data.root";
      std::cout << "Set unfold filename " << _unfold_file_name << std::endl;
    }

    _unfold_file = TFile::Open(_unfold_file_name.c_str(),"UPDATE");
    std::cout << "UPDATE: " << _unfold_file_name << std::endl;
    
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

    prepare_unfold_file_parameters();
    prepare_unfold_file_responses();
    
    return true;
  }

  void SPIO::dump_branches() {
    std::cout << "tree called " << _mc_tree_name << " @ chain" << _in_tree << " has branches:" << std::endl;;
    std::cout << "[" << std::endl;
    for(const auto& branch : _in_mc_branch_v) 
      std::cout << branch << "," << std::endl;
    std::cout << "]" << std::endl;
  }
  
  
  void SPIO::add_true_parameter(const std::string& true_var_name,
				const std::vector<double>& bin_lo_v,
				const std::vector<double>& bin_hi_v) {
    Parameter in_param(true_var_name,bin_lo_v,bin_hi_v);
    
    const auto unfold_param = search_unfold_parameters(in_param);

    if (unfold_param) {
      _true_parameter_v.emplace_back(*unfold_param);
      _true_parameter_v.back()._from_file = true;
    }
    else {
      _true_parameter_v.emplace_back(std::move(in_param));
    }

    std::cout << "Now P @ " << &_true_parameter_v.back() << std::endl;
    std::cout << "From file? " << _true_parameter_v.back()._from_file << std::endl;
  }
  
  void SPIO::add_reco_parameter(const std::string reco_var_name,
				const std::vector<double>& bin_lo_v,
				const std::vector<double>& bin_hi_v)  {
    Parameter in_param(reco_var_name,bin_lo_v,bin_hi_v);

    const auto unfold_param = search_unfold_parameters(in_param);

    if (unfold_param) {
      _reco_parameter_v.emplace_back(*unfold_param);
      _reco_parameter_v.back()._from_file = true;
    }
    else {
      _reco_parameter_v.emplace_back(std::move(in_param));
    }

    std::cout << "Now P @ " << &_reco_parameter_v.back() << std::endl;
    std::cout << "From file? " << _reco_parameter_v.back()._from_file << std::endl;
  }
  
  const Parameter* SPIO::search_unfold_parameters(const Parameter& in_param) {
    for(size_t pid = 0; pid < _unfold_parameter_ptr_v.size(); ++pid) {
      if (in_param == *(_unfold_parameter_ptr_v[pid])) {
	std::cout << "Found parameter in file @ position " << pid << " return " << _unfold_parameter_ptr_v[pid] << std::endl;
	return _unfold_parameter_ptr_v[pid];
      }
    }
    std::cout << "No parameter found in file!" << std::endl;
    return nullptr;
  }
  

  void SPIO::prepare_unfold_file_parameters() {
    std::cout << std::endl;
    std::cout << "Scan existing paramters" << std::endl;
    
    auto dir = (TDirectory*)_unfold_file->GetDirectory("Parameters");

    if(!dir) { std::cout << "No Parameters folder found" << std::endl; return;}
    dir->cd();

    TList *objs = dir->GetListOfKeys();
    
    std::cout << "GOT: " << objs->GetEntries() << " parameter definitions" << std::endl;
    if(!objs->GetEntries()) return;

    for(size_t param_id=0; param_id < objs->GetEntries(); param_id++) {
      const auto ref = (const Parameter*) dir->Get(objs->At(param_id)->GetName());
      if(!ref) throw sperr("Failed reading parameter!");
      _unfold_parameter_ptr_v.push_back(ref);
      std::cout << "Found parameter @ ("<< param_id << "," << _unfold_parameter_ptr_v.back() << ")" << std::endl;
    }
  }


  void SPIO::prepare_unfold_file_responses() {
    std::cout << std::endl;
    std::cout << "Scan existing responses" << std::endl;
    
    auto dir = (TDirectory*)_unfold_file->GetDirectory("Responses");

    if(!dir) { std::cout << "No Responses folder found" << std::endl; return;}
    dir->cd();

    TList *objs = dir->GetListOfKeys();

    std::cout << "Got: " << objs->GetEntries() << " response definitions" << std::endl;
    if(!objs->GetEntries()) return;
    
    for(size_t resp_id=0; resp_id < objs->GetEntries(); resp_id++) {
      const auto ref = (const Response*) dir->Get(objs->At(resp_id)->GetName());
      if(!ref) throw sperr("Failed reading parameter!");
      _unfold_response_ptr_v.push_back(ref);
      std::cout << "Found response @ ("<< resp_id << "," << _unfold_response_ptr_v.back() << ")" << std::endl;
    }
  }
  
  
  bool SPIO::write_unfold_file() {
    auto dir = (TDirectory*)_unfold_file->GetDirectory("Parameters");
    if (dir) dir->cd();
    else {
      auto new_dir = (TDirectory*)_unfold_file->mkdir("Parameters");
      if(!new_dir) throw sperr("Directory Parameters could not be made");
      new_dir->cd();
    }

    for(const auto& param : _true_parameter_v) {
      if (param._from_file) { std::cout << "SKIP parameter @ " << &param << std::endl; continue;}
      std::cout << "WRITE: " << param._name << std::endl;
      auto res = param.Write(param._name.c_str());
      if(!res) throw sperr("Could not write true parameter object");
    }
        
    for(const auto& param : _reco_parameter_v) {
      if (param._from_file) { std::cout << "SKIP parameter @ " << &param << std::endl; continue;}
      auto res = param.Write(param._name.c_str());
      std::cout << "WRITE: " << param._name << std::endl;
      if(!res) throw sperr("Could not write reco parameter object");
    }

    dir = (TDirectory*)_unfold_file->GetDirectory("Responses");
    if (dir) dir->cd();
    else {
      auto new_dir = (TDirectory*)_unfold_file->mkdir("Responses");
      if(!new_dir) throw sperr("Directory Responses could not be made");
      new_dir->cd();
    }
    
    for(const auto& response : _response_v) {
      if (response.from_file) { std::cout << "SKIP response @ " << &response << std::endl; continue;}
      std::cout << "WRITE: " << response.name << std::endl;
      auto res = response.Write(response.name.c_str());
      if(!res) throw sperr("Could not write response object");
    }
    
    _unfold_file->Write();
    _unfold_file->Close();
    return true;
  }

  bool SPIO::init_response_matrix() {
    std::cout << std::endl;
    std::cout << "Instantiate response matrix" << std::endl;
    
    for(size_t true_id = 0; true_id < _true_parameter_v.size(); ++true_id) {
      auto& true_param = _true_parameter_v[true_id];
      for(size_t reco_id = 0; reco_id < _reco_parameter_v.size(); ++reco_id) {
	auto& reco_param = _reco_parameter_v[reco_id];

	Response response(&true_param,&reco_param);

	auto res = search_unfold_responses(response);
	if (res) {
	  std::cout << "Got from file: " << res->name << std::endl;
	  std::cout << "Filling " << res->name << " response @ " << res << std::endl;
	  _response_v.emplace_back(*res);
	  _response_v.back().from_file = true;
	}
	else  {
	  std::stringstream ss;
	  ss << "response_" << true_id << "_" << reco_id << "_" << response.true_param->_name << "_" << response.reco_param->_name;
	  std::cout << "Made new " << ss.str() << std::endl;
	  response.name = ss.str();
	  std::cout << "Filling " << response.name << " response @ " << &response << std::endl;
	  _response_v.emplace_back(std::move(response));
	}
	auto& this_res = _response_v.back();

	this_res.true_param = &true_param;
	this_res.reco_param = &reco_param;

	this_res.true_param->_response = &this_res;
	this_res.reco_param->_response = &this_res;

	std::cout << "Now R @ " << &this_res << std::endl;	  
      }
    }
    return true;
  }

  const Response* SPIO::search_unfold_responses(const Response& in_response) {
    for(size_t rid = 0; rid < _unfold_response_ptr_v.size(); ++rid) {
      if (in_response == *(_unfold_response_ptr_v[rid])) {
	std::cout << "Found response in file @ position " << rid << " return " << _unfold_response_ptr_v[rid] << std::endl;
	return _unfold_response_ptr_v[rid];
      }
    }
    std::cout << "No response found in file!" << std::endl;
    return nullptr;
  }

  
  bool SPIO::fill_responses() {

    std::vector<Response*> unfilled_response_v;
    unfilled_response_v.reserve(_response_v.size());

    for(auto& res : _response_v) {
      if (res.filled()) continue;
      unfilled_response_v.emplace_back(&res);
    }
    
    if (unfilled_response_v.empty()) 
      { std::cout << "All responses filled!" << std::endl; return true; }

    for(size_t rid = 0; rid< _response_v.size(); ++rid) {
      auto& response = _response_v[rid];
      _in_tree->SetBranchAddress(response.true_param->_name.c_str(),&response.true_param->_data);
      _in_tree->SetBranchAddress(response.reco_param->_name.c_str(),&response.reco_param->_data);
    }

    float weight;
    bool passosc;
    int nutype;
    
    _in_tree->SetBranchAddress("Weight" , &weight);
    _in_tree->SetBranchAddress("PassOsc", &passosc);
    _in_tree->SetBranchAddress("NuType" , &nutype);
    
    std::cout << "READING: " << _in_n_entries << " entries from MC file" << std::endl;
    for(size_t entry = 0; entry < _in_n_entries; ++entry) {
      _in_tree->GetEntry(entry);

      if (entry>=10000) break;
      if (weight<=0) continue;

      for(auto& res : unfilled_response_v)
	res->Fill(weight,passosc,nutype);

    }

    for(auto& res : unfilled_response_v)
      res->Finalize();
    
    return true;
  }


  
}

#endif
