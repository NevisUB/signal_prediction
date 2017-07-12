#ifndef __STORAGE_CXX__
#define __STORAGE_CXX__

#include "Storage.h"

namespace sp {

  //
  // Parameter
  //
  Parameter::Parameter() : 
    _name(""), 
    _bin_lo_v(std::vector<double>()),
    _bin_hi_v(std::vector<double>()), 
    _filled(false),
    _from_file(false),
    _hist(TH1D()), 
    _response(nullptr), 
    _data(0.0)
  { std::cout << "P @ " << this << std::endl; }

  Parameter::Parameter(const std::string& name,
		       const std::vector<double>& bin_lo_v,
		       const std::vector<double>& bin_hi_v) :
    _name(name),
    _bin_lo_v(bin_lo_v),
    _bin_hi_v(bin_hi_v),
    _filled(false),
    _from_file(false),
    _response(nullptr),
    _data(0.0)
  {
    std::cout << "P @ " << this << std::endl;
    _hist = TH1D(_name.c_str(),"",(int)(bin_lo_v.size()-1),_bin_lo_v.data());
    _hist.SetDirectory(0);
    this->dump();
  }

  void Parameter::Fill(float weight) {
    _hist.Fill(_data,weight);
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

    std::cout << "TH1D pointer : " << &_hist << std::endl;
    std::cout << "TH1D n_entries : " << _hist.GetEntries() << std::endl;

    std::cout << "From File: " << _from_file << std::endl;
    std::cout << "Filled: " << _filled << std::endl;
    std::cout << std::endl;
  }



  //
  // Response
  //
  Response::Response() :
    name(""),
    true_param(nullptr),
    reco_param(nullptr),
    response(TMatrixD()),
    true_v(TVectorD()),
    reco_v(TVectorD()),
    response_h(TH2D()),
    from_file(false)
  { std::cout << "R @ " << this << std::endl; }


  Response::Response(Parameter* true_p, Parameter* reco_p) :
    name(""),
    true_param(true_p),
    reco_param(reco_p),
    response(TMatrixD()),
    true_v(TVectorD()),
    reco_v(TVectorD()),
    response_h(TH2D()),
    from_file(false)
  { 
    std::cout << "R @ " << this << std::endl;
    response_h = TH2D(name.c_str(),"",
		      (int)(true_param->_bin_lo_v.size()-1),true_param->_bin_lo_v.data(),
		      (int)(reco_param->_bin_lo_v.size()-1),reco_param->_bin_lo_v.data());
    response_h.SetDirectory(0);
    this->dump();
  }

  
  void Response::Fill(float weight, bool passosc, int nutype) {
    response_h.Fill(true_param->_data,reco_param->_data,weight);
    true_param->Fill(weight);
    reco_param->Fill(weight);
  }

  void Response::Finalize() {
    response.ResizeTo(response_h.GetNbinsX()+2,response_h.GetNbinsY()+2);
    response = TMatrixD(response_h.GetNbinsX()+2,response_h.GetNbinsY()+2,response_h.GetArray());

    true_param->_filled = true;
    reco_param->_filled = true;
  }

  bool Response::filled() const {
    if(!true_param) return false;
    if(!reco_param) return false;

    if(!true_param->_filled) return false;
    if(!reco_param->_filled) return false;
    
    return true;
  }

  void Response::dump() const {
    std::cout << std::endl;
    std::cout << "Response @ " << this << std::endl;

    std::cout << "name : " << name << std::endl;
    std::cout << "True parameter" << std::endl;
    if(true_param) true_param->dump();
    else std::cout << "NULL" << std::endl;
    
    std::cout << "Reco parameter" << std::endl;
    if(reco_param) reco_param->dump();
    else std::cout << "NULL" << std::endl;

    std::cout << "TH2D n_entries : " << response_h.GetEntries() << std::endl;

    std::cout << "From file : " << from_file << std::endl;
    std::cout << "Filled : " << this->filled() << std::endl;
    std::cout << std::endl;
  }
  
}

#endif
