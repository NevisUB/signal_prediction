#ifndef __RESPONSE_CXX__
#define __RESPONSE_CXX__

#include "Response.h"

namespace sp {

  //
  // Response
  //
  Response::Response() :
    _name(""),
    _true_param(nullptr),
    _reco_param(nullptr),
    _response(TMatrixD()),
    _response_h(TH2D()),
    _from_file(false)
  {}


  Response::Response(Parameter* true_p, Parameter* reco_p) :
    _name(""),
    _true_param(true_p),
    _reco_param(reco_p),
    _response(TMatrixD()),
    _response_h(TH2D()),
    _from_file(false)
  { 
    std::cout << "sp::Repsonse::Response(true, false) || Creating a response with "<<_true_param->_bin_lo_v.size()<<" bin edges in truth and "<<_reco_param->_bin_lo_v.size()<<" bin edges in reco"<<std::endl;
    _response_h = TH2D(_name.c_str(),"",
		       (int)(_reco_param->_bin_lo_v.size()-1),_reco_param->_bin_lo_v.data(),
		       (int)(_true_param->_bin_lo_v.size()-1),_true_param->_bin_lo_v.data());
    _response_h.SetDirectory(0);
  }


  void Response::Fill(float weight, bool passosc) {

    auto t_res = _true_param->Fill(weight);
    if(passosc){
      auto r_res = _reco_param->Fill(weight);
      _response_h.Fill(r_res, t_res, weight);
    }

  }

  void Response::Finalize() {
    _response.ResizeTo(_response_h.GetNbinsX()+2,_response_h.GetNbinsY()+2);
		
    // Direct using GetArray caused some problems oddly..
    // _response = TMatrixD(_response_h.GetNbinsX()+2,_response_h.GetNbinsY()+2, _response_h.GetArray());
    for(int i=0; i < _response_h.GetNbinsX()+2; i++){
      for(int a=0; a < _response_h.GetNbinsY()+2; a++){
	_response(i,a) = _response_h.GetBinContent(i,a);
      }
    }

    _true_param->_filled = true;
    _reco_param->_filled = true;
  }

  bool Response::filled() const {
    if(!_true_param) return false;
    if(!_reco_param) return false;

    if(!_true_param->_filled) return false;
    if(!_reco_param->_filled) return false;
    return true;
  }

  void Response::dump() const {
    std::cout << std::endl;
    std::cout << "Response @ " << this << std::endl;

    std::cout << "name : " << _name << std::endl;
    std::cout << "True parameter" << std::endl;
    if(_true_param) _true_param->dump();
    else std::cout << "NULL" << std::endl;

    std::cout << "Reco parameter" << std::endl;
    if(_reco_param) _reco_param->dump();
    else std::cout << "NULL" << std::endl;

    std::cout << "TH2D n_entries : " << _response_h.GetEntries() << std::endl;

    std::cout << "From file : " << _from_file << std::endl;
    std::cout << "Filled : " << this->filled() << std::endl;
    std::cout << std::endl;
  }

}
#endif