#ifndef __SPIO_CXX__
#define __SPIO_CXX__

#include "SPIO.h"
#include "SPErr.h"
#include <sstream>

#include "Combined/CombinedUtil.h"
#include "Combined/CombinedFunctions.h"

namespace sp {

  SPIO::SPIO() : 
    sp_base("SPIO"),
    _mc_tree_name(""),
    _unfold_file_name(""),
    _in_tree(nullptr),
    _unfold_file(nullptr),
    _model(nullptr)
  {}

  SPIO::~SPIO() {}
  
  void SPIO::add_mc_in_file(const std::string& filename) {
    _in_mc_file_v.push_back(filename);
  }
  
  void SPIO::set_mc_tree_name(const std::string& treename) {
    _mc_tree_name = treename;
  }

  void SPIO::set_unfold_in_file(const std::string& filename) {
    _unfold_file_name = filename;
  }

  std::vector<TH1D> SPIO::gen_background(const std::string& filename, 
					 const std::string& param,
					 const std::vector<double>& bins_lo_v) {

    SP_DEBUG() << "start" << std::endl;
    auto tf = TFile::Open(filename.c_str(),"READ");

    int NFSP;
    std::vector<int>* FSPType_v = nullptr;
    
    std::vector<float>* VertexX_v = nullptr;
    std::vector<float>* VertexY_v = nullptr;
    std::vector<float>* VertexZ_v = nullptr;

    std::vector<float>* MomX_v = nullptr;
    std::vector<float>* MomY_v = nullptr;
    std::vector<float>* MomZ_v = nullptr;
    std::vector<float>* MomT_v = nullptr;
  

    float Weight;    
    float Data;

    bool PassOsc;

    int NUANCEChan;
    int NuType;
    int NuParentID;
    
    auto tree = (TTree*) tf->Get("MiniBooNE_CCQE");
    SP_DEBUG() << "See MiniBooNE_CCQE @ " << tree << std::endl;
    if (tree == nullptr) throw sperr("Bad treename from file");
    
    SP_DEBUG() << "Got " << tree->GetEntries() << " @ tree=" << tree << " called " << tree->GetName() << std::endl;

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

    tree->SetBranchAddress("Weight" , &Weight);
    
    if(1) {}
    else if (param == "RecoEnuQE") {}
    else if (param == "CosTheta") {}
    else if (param == "Energy") {}
    else throw sperr("Did not specify a valid TTree variable");

    SP_DEBUG() << "Set param " << param << " @ " << &Data << std::endl;
    tree->SetBranchAddress(param.c_str(), &Data);
  
    std::vector<TH1D> res_v((size_t)sp::kBKGD_MAX);
    for(size_t res_id = 0; res_id < res_v.size(); ++res_id) {
      auto name = sp::StackedBkgd2String((sp::StackedBkgdType_t)res_id);
      SP_DEBUG() << "Set histo @ " << res_id << " named " << name << std::endl;
      res_v[res_id] = TH1D(name.c_str(),
			   ";;",
			   bins_lo_v.size()-1,
			   bins_lo_v.data());
    }

    for(size_t entry = 0; entry < (size_t)tree->GetEntries(); ++entry) {
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
      sp::StackedBkgdType_t bkg_type = sp::StackHistoBkgd(0,
							  (bool)isPi0,
							  (sp::NuanceType_t)NUANCEChan,
							  (sp::NuType_t)NuType,
							  (sp::GEANT3Type_t)NuParentID);
      
      auto& res = res_v.at((size_t) bkg_type);
      res.Fill(Data,Weight);
    }
  
    tf->Close();
    return res_v;
  }
  

  bool SPIO::initialize() {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Initialize" << std::endl;
    TH1::AddDirectory(kFALSE);

    if(_unfold_file_name.empty())  {
      SP_DEBUG() << "Unfold filename not set" << std::endl;
      _unfold_file_name = "unfold_data.root";
      SP_DEBUG() << "Set unfold filename " << _unfold_file_name << std::endl;
    }

    _unfold_file = TFile::Open(_unfold_file_name.c_str(),"UPDATE");
    SP_DEBUG() << "UPDATE: " << _unfold_file_name << std::endl;
    
    if (!_in_mc_file_v.empty()) {
      if(_mc_tree_name.empty()) throw sperr("No treename specified");
      
      _in_tree = new TChain(_mc_tree_name.c_str());
      for(const auto& fname : _in_mc_file_v) {
	SP_DEBUG() << "READ: " << fname << std::endl;
	_in_tree->Add(fname.c_str());
      }

      _in_n_entries = _in_tree->GetEntries();
      if(!_in_n_entries) throw sperr("Empty mc tree specified (NO ENTRIES)");
      
      TObjArray *branch_list = _in_tree->GetListOfBranches();
      if(!branch_list->GetEntries()) throw sperr("Empty mc tree specified (NO BRANCHES)");
      
      for(size_t branch_id=0;branch_id<(size_t)branch_list->GetEntries(); branch_id++)
	_in_mc_branch_v.insert(branch_list->At(branch_id)->GetName());
    }

    prepare_unfold_file_parameters();
    prepare_unfold_file_responses();
    SP_DEBUG() << std::endl;
    return true;
  }
  
  
  void SPIO::add_true_parameter(const std::vector<std::string>& var_v,
				const std::vector<double>& bin_lo_v,
				Operation_t op) {
  
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Requested to add true parameter" << std::endl;
    Parameter in_param(var_v,bin_lo_v,op);
    
    const auto unfold_param = search_unfold_parameters(in_param);

    if (unfold_param) {
      _true_parameter_v.emplace_back(*unfold_param);
      _true_parameter_v.back()._from_file = true;
    }
    else {
      _true_parameter_v.emplace_back(std::move(in_param));
    }
    SP_DEBUG() << "Now P @ " << &_true_parameter_v.back() << " from file: " << _true_parameter_v.back()._from_file << std::endl;
    SP_DEBUG() << std::endl;
  }
  
  void SPIO::add_reco_parameter(const std::vector<std::string>& var_v,
				const std::vector<double>& bin_lo_v,
				Operation_t op) {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Requested to add reco parameter" << std::endl;
    Parameter in_param(var_v,bin_lo_v,op);

    const auto unfold_param = search_unfold_parameters(in_param);

    if (unfold_param) {
      _reco_parameter_v.emplace_back(*unfold_param);
      _reco_parameter_v.back()._from_file = true;
    }
    else {
      _reco_parameter_v.emplace_back(std::move(in_param));
    }
    SP_DEBUG() << "Now P @ " << &_reco_parameter_v.back() << " from file: " << _reco_parameter_v.back()._from_file << std::endl;
    SP_DEBUG() << std::endl;
  }
  
  const Parameter* SPIO::search_unfold_parameters(const Parameter& in_param) {
    for(size_t pid = 0; pid < _unfold_parameter_ptr_v.size(); ++pid) {
      if (in_param == *(_unfold_parameter_ptr_v[pid])) {
	SP_DEBUG() << "Found parameter in file @ position " << pid << " return " << _unfold_parameter_ptr_v[pid] << std::endl;
	return _unfold_parameter_ptr_v[pid];
      }
    }
    SP_DEBUG() << "No parameter found in file!" << std::endl;
    return nullptr;
  }
  
  
  void SPIO::prepare_unfold_file_parameters() {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Scan existing Paramters" << std::endl;
    
    auto dir = (TDirectory*)_unfold_file->GetDirectory("Parameters");

    if(!dir) { SP_DEBUG() << "No Parameters folder found" << std::endl; return;}
    dir->cd();

    TList *objs = dir->GetListOfKeys();
    
    SP_DEBUG() << "GOT: " << objs->GetEntries() << " parameter definitions" << std::endl;
    if(!objs->GetEntries()) return;

    for(size_t param_id=0; param_id < (size_t)objs->GetEntries(); param_id++) {
      const auto ref = (const Parameter*) dir->Get(objs->At(param_id)->GetName());
      if(!ref) throw sperr("Failed reading parameter!");
      _unfold_parameter_ptr_v.push_back(ref);
      SP_DEBUG() << "Found parameter @ ("<< param_id << "," << _unfold_parameter_ptr_v.back() << ")" << std::endl;
    }
    SP_DEBUG() << std::endl;
  }


  void SPIO::prepare_unfold_file_responses() {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Scan existing Responses" << std::endl;
    
    auto dir = (TDirectory*)_unfold_file->GetDirectory("Responses");

    if(!dir) { SP_DEBUG() << "No Responses folder found" << std::endl; return;}
    dir->cd();

    TList *objs = dir->GetListOfKeys();

    SP_DEBUG() << "Got: " << objs->GetEntries() << " response definitions" << std::endl;
    if(!objs->GetEntries()) return;
    
    for(size_t resp_id=0; resp_id < (size_t)objs->GetEntries(); resp_id++) {
      const auto ref = (const Response*) dir->Get(objs->At(resp_id)->GetName());
      if(!ref) throw sperr("Failed reading parameter!");
      _unfold_response_ptr_v.push_back(ref);
      SP_DEBUG() << "Found response @ ("<< resp_id << "," << _unfold_response_ptr_v.back() << ")" << std::endl;
    }
    SP_DEBUG() << std::endl;
  }
  

  bool SPIO::init_response_matrix() {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Instantiate response matrix" << std::endl;

    for(size_t true_id = 0; true_id < _true_parameter_v.size(); ++true_id) {
      auto& true_param = _true_parameter_v[true_id];
      for(size_t reco_id = 0; reco_id < _reco_parameter_v.size(); ++reco_id) {
	auto& reco_param = _reco_parameter_v[reco_id];

	Response response(&true_param,&reco_param);

	auto res = search_unfold_responses(response);
	if (res) {
	  SP_DEBUG() << "Got from file: " << res->_name << std::endl;
	  SP_DEBUG() << "Filling " << res->_name << " response @ " << res << std::endl;
	  _response_v.emplace_back(*res);
	  _response_v.back()._from_file = true;
	auto& this_res = _response_v.back();

	SP_DEBUG() << "Set true_param for response @ " << &true_param << std::endl;
	SP_DEBUG() << "Set reco_param for response @ " << &reco_param << std::endl;
	this_res._true_param = &true_param;
	this_res._reco_param = &reco_param;

	SP_DEBUG() << "Now R @ " << &this_res << std::endl;	  
      }
    }
    SP_DEBUG() << std::endl;
    return true;
  }

  const Response* SPIO::search_unfold_responses(const Response& in_response) {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Search unfold responses for... " << &in_response << std::endl;
    for(size_t rid = 0; rid < _unfold_response_ptr_v.size(); ++rid) {
      if (in_response == *(_unfold_response_ptr_v[rid])) {
	SP_DEBUG() << "Found response in file @ position " << rid << " return " << _unfold_response_ptr_v[rid] << std::endl;
	return _unfold_response_ptr_v[rid];
      }
    }
    SP_DEBUG() << "No response found in file" << std::endl;
    SP_DEBUG() << std::endl;
    return nullptr;
  }
  
  bool SPIO::fill_responses() {
    SP_DEBUG() << std::endl;
    SP_DEBUG() << "Fill responses" << std::endl;
    std::vector<Response*> unfilled_response_v;
    unfilled_response_v.reserve(_response_v.size());

    for(auto& res : _response_v) {
      if (res.filled()) continue;
      SP_DEBUG() << "... response @ " << &res << " unfilled!"  << std::endl;
      unfilled_response_v.emplace_back(&res);
    }
    
    if (unfilled_response_v.empty())  {
      SP_DEBUG() << "All responses filled, return" << std::endl;
      return true;
    }

    for(size_t rid = 0; rid < unfilled_response_v.size(); ++rid) {
      auto response = unfilled_response_v[rid];
      SP_DEBUG() << std::endl;
      SP_DEBUG() << rid << ") @ response " << response << std::endl;
      for(size_t vid = 0; vid < response->_true_param->_variable_v.size(); ++vid) {
	SP_DEBUG() << "SET: " << response->_true_param 
		  << " true branch " << response->_true_param->_variable_v[vid] 
		  << " @ " << &response->_true_param->_data_v[vid] << std::endl;
	_in_tree->SetBranchAddress(response->_true_param->_variable_v[vid].c_str(),&response->_true_param->_data_v[vid]);
      }

      for(size_t vid = 0; vid < response->_reco_param->_variable_v.size(); ++vid) {
	SP_DEBUG() << "SET: " << response->_reco_param << " reco branch " << response->_reco_param->_variable_v[vid] << " @ " << &response->_reco_param->_data_v[vid] << std::endl;
	_in_tree->SetBranchAddress(response->_reco_param->_variable_v[vid].c_str(),&response->_reco_param->_data_v[vid]);
      }

      float weight;
      bool passosc;
      _in_tree->SetBranchAddress("Weight"  , &weight);
      _in_tree->SetBranchAddress("PassOsc" , &passosc);
      
      if (!_model) throw sperr("No Model specified");

      _in_tree->SetBranchAddress("NuType"    , &_model->NuType);
      _in_tree->SetBranchAddress("NUANCEChan", &_model->NUANCEChan);

      SP_DEBUG() << "READING: " << _in_n_entries << " entries from MC file" << std::endl;
      for(size_t entry = 0; entry < _in_n_entries; ++entry) {
	_in_tree->GetEntry(entry);
	
	if (_model->Valid())
	  response->Fill(weight,passosc);
	
      }

      response->Finalize();
    }

    SP_DEBUG() << std::endl; 
    return true;
  }

  bool SPIO::write_unfold_file() {
    SP_DEBUG() << "start" << std::endl;

    auto dir = (TDirectory*)_unfold_file->GetDirectory("Parameters");
    if (dir) dir->cd();
    else {
      auto new_dir = (TDirectory*)_unfold_file->mkdir("Parameters");
      if(!new_dir) throw sperr("Directory Parameters could not be made");
      new_dir->cd();
    }

    for(const auto& param : _true_parameter_v) {
      if (param._from_file) { SP_DEBUG() << "SKIP parameter @ " << &param << std::endl; continue;}
      SP_DEBUG() << "WRITE: " << param._name << std::endl;
      auto res = param.Write(param._name.c_str());
      if(!res) throw sperr("Could not write true parameter object");
    }
        
    for(const auto& param : _reco_parameter_v) {
      if (param._from_file) { SP_DEBUG() << "SKIP parameter @ " << &param << std::endl; continue;}
      auto res = param.Write(param._name.c_str());
      SP_DEBUG() << "WRITE: " << param._name << std::endl;
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
      if (response._from_file) { SP_DEBUG() << "SKIP response @ " << &response << std::endl; continue;}
      SP_DEBUG() << "WRITE: " << response._name << std::endl;
      auto res = response.Write(response._name.c_str());
      if(!res) throw sperr("Could not write response object");
    }
    
    SP_DEBUG() << "write" << std::endl;
    _unfold_file->Write();
    
    SP_DEBUG() << "close" << std::endl;
    _unfold_file->Close();

    SP_DEBUG() << "done" << std::endl; 
    return true;
  }

  void SPIO::dump_branches() {
    SP_DEBUG() << "tree called " << _mc_tree_name << " @ chain" << _in_tree << " has branches:" << std::endl;;
    SP_DEBUG() << "[" << std::endl;
    for(const auto& branch : _in_mc_branch_v) 
      SP_DEBUG() << branch << "," << std::endl;
    SP_DEBUG() << "]" << std::endl;
  }
  

}

#endif
