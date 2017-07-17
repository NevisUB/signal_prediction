#ifndef __PARAMETER_CXX__
#define __PARAMETER_CXX__

#include "Parameter.h"
#include "UnfoldUtil.h"

namespace sp {

  Parameter::Parameter() : 
    _name(""), 
    _variable_v(std::vector<std::string>()),
    _bin_lo_v(std::vector<double>()),
    _operation(kOP_INVALID),
    _filled(false),
    _from_file(false),
    _hist(TH1D()), 
    _data_v(std::vector<float>())
  { 
    std::cout << "P @ " << this << std::endl; 
  }
  
  Parameter::Parameter(const std::vector<std::string>& variable_v,
		       const std::vector<double>& bin_lo_v,
		       Operation_t op) :
    _name(""),
    _variable_v(variable_v),
    _bin_lo_v(bin_lo_v),
    _operation(op),
    _filled(false),
    _from_file(false),
    _data_v(std::vector<float>(variable_v.size(),0.0))
  {
	  
    std::cout << "P @ " << this << std::endl;
    _name = to_name(op,variable_v,bin_lo_v);
	  
    _hist = TH1D(concatenate(variable_v).c_str(),"",(int)(bin_lo_v.size()-1),_bin_lo_v.data());
    _hist.SetDirectory(0);
  }
    
  float Parameter::Fill(float weight) {
    auto res = Operate(_data_v,_operation);
    if (!_filled)
      _hist.Fill(res,weight);
    return res;
  }
  
  void Parameter::dump() {
    std::cout << std::endl;
    std::cout << "Parameter @ " << this << std::endl;

    std::cout << "name :    " << _name << std::endl;
    std::cout << "operation :    " << (unsigned)_operation << std::endl;

    std::cout << "variable_v: [";
    for(auto v : _variable_v) std::cout << v << ",";
    std::cout << "]" << std::endl;

    std::cout << "bin_lo_v: [";
    for(auto v : _bin_lo_v) std::cout << v << ",";
    std::cout << "]" << std::endl;
    
    std::cout << "TH1D pointer : " << &_hist << std::endl;
    std::cout << "TH1D n_entries : " << _hist.GetEntries() << std::endl;

    std::cout << "From File: " << _from_file << std::endl;
    std::cout << "Filled: " << _filled << std::endl;
    std::cout << std::endl;
  }
}
#endif
