#ifndef __UNFOLDUTIL_CXX__
#define __UNFOLDUTIL_CXX__

#include <sstream>

#include "UnfoldUtil.h"
#include "Base/SPErr.h"

namespace sp {

  std::vector<double> Subtract(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) throw sperr("Given vectors are different sizes!");
    std::vector<double> res(v1.size(),0.0);
    for(size_t i=0; i<v1.size(); ++i) res[i] = v1[i] - v2[i];
    return res;
  }

  std::string concatenate(const std::vector<std::string>& string_v) {
    std::stringstream ss;
    for(size_t sid=0; sid<string_v.size(); ++sid) {
      ss << string_v[sid];
      if (sid != string_v.size()-1) ss << "_";
    }
    return ss.str();
  }

  std::string to_name(Operation_t operation,
		      const std::vector<std::string>& prefix_v,
		      const std::vector<double>& bin_lo_v) {
    
    std::stringstream ss;
    ss.precision(3);

    ss << (unsigned)operation << "_";
    for(const auto& v : prefix_v) ss << v << "_";
    for(const auto& v : bin_lo_v) ss << v << "_";
    return ss.str();
  }


  float GeV2MeV(float gev) { return (1000.0 * gev); }
  float MeV2GeV(float mev) { return (mev / 1000.0); }

  // for now it has to be the same was as function definition...
  float Operate(const std::vector<float>& data, Operation_t operation) {

    switch(operation) {
      
    case kOP_INVALID : 
      { 
	return data.front(); 
      }
    case kOP_EQE : 
      { 
	float res = CombinedFit_EnuQE_ryan(data.at(0),data.at(1),kNUE);
	res *= 1000.0;
	return res;
      }

    case kOP_Q2 : 
      { 
	return CombinedFit_Q2_ryan(data.at(0),data.at(1),kNUE);
      }
    case kOP_MEV2GEV :
      {
	return MeV2GeV(data.front());
      }
    case kOP_GEV2MEV :
      {
	return GeV2MeV(data.front());
      }
    default:
      {
	throw sperr("Unknown operation observed");
      }
    }
    return kINVALID_FLOAT;
  }
  
  TVectorD hist2TVec(TH1D * hist){
    int dim = hist->GetNbinsX()+2;
    TVectorD tmp(dim);

    for(int i = 0; i<dim; i++){
      tmp[i] = hist->GetBinContent(i);
    }

    return tmp;	
  }




}

#endif 
