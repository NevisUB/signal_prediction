#ifndef __STORAGE_CXX__
#define __STORAGE_CXX__

#include "Storage.h"
#include <iostream>

namespace sp {

  //
  // Parameter
  //
  Parameter::Parameter(const std::string& name,
		       const std::vector<double>& bin_lo_v,
		       const std::vector<double>& bin_hi_v,
		       std::vector<std::string> options_v) :
    _name(name),
    _bin_lo_v(bin_lo_v),
    _bin_hi_v(bin_hi_v),
    _filled(false),
    _hist(nullptr)
  {
    _hist = new TH1D(_name.c_str(),"",(int)(bin_lo_v.size()-1),_bin_lo_v.data());
    this->dump();
  }

  void Parameter::Fill(float weight) {
    _hist->Fill(_data,weight);
  }
  
  void Parameter::dump() {
    std::cout << std::endl;
    std::cout << "Parameter @ " << this << std::endl;

    std::cout << "name : " << _name << std::endl;

    std::cout << "bin_lo_v: [";
    for(auto v : _bin_lo_v) std::cout << v << ",";
    std::cout << "]" << std::endl;
    
    std::cout << "bin_hi_v: [";
    for(auto v : _bin_hi_v) std::cout << v << ",";
    std::cout << "]" << std::endl;

    std::cout << "bin_w_v: [";
    for(auto v : _bin_w_v) std::cout << v << ",";
    std::cout << "]" << std::endl;

    if(_hist) {
      std::cout << "TH1D pointer : " << _hist << std::endl;
      auto copy = TH1D(*_hist);
      std::cout << "TH1D name : " << copy.GetName() << std::endl;
      std::cout << "TH1D n_entries : " << (*_hist).GetEntries() << std::endl;
    }
    
    std::cout << "Filled: " << _filled << std::endl;
    std::cout << std::endl;
  }



  //
  // Response
  //
  Response::Response(Parameter* true_p, Parameter* reco_p) :
    name(""),
    true_param(true_p),
    reco_param(reco_p),
    response_h(nullptr)
  { 
    response_h = new TH2D(name.c_str(),"",
			  (int)(true_param->_bin_lo_v.size()-1),true_param->_bin_lo_v.data(),
			  (int)(reco_param->_bin_lo_v.size()-1),reco_param->_bin_lo_v.data());
    this->dump();
  }

  
  void Response::Fill(float weight, bool passosc, int nutype) {
    response_h->Fill(true_param->_data,reco_param->_data,weight);
    true_param->Fill(weight);
    reco_param->Fill(weight);
  }

  void Response::Finalize() {
    response.ResizeTo(response_h->GetNbinsX()+2,response_h->GetNbinsY()+2);
    response = TMatrixD(response_h->GetNbinsX()+2,response_h->GetNbinsY()+2,response_h->GetArray());
  }

  void Response::dump() {
    std::cout << std::endl;
    std::cout << "Response @ " << this << std::endl;

    std::cout << "name : " << name << std::endl;
    std::cout << "True parameter" << std::endl;
    if(true_param) true_param->dump();
    else std::cout << "NULL" << std::endl;
    
    std::cout << "Reco parameter" << std::endl;
    if(reco_param) reco_param->dump();
    else std::cout << "NULL" << std::endl;

    if(response_h) std::cout << "TH2D n_entries : " << response_h->GetEntries() << std::endl;

    std::cout << std::endl;
  }
  
}

#endif
