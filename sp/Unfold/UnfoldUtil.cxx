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
		      const std::vector<double>& bin_lo_v,
		      const std::vector<double>& bin_hi_v) {
    
    std::stringstream ss;
    ss.precision(3);

    ss << (unsigned)operation << "_";
    for(const auto& v : prefix_v) ss << v << "_";
    for(const auto& v : bin_lo_v) ss << v << "_";
    for(const auto& v : bin_hi_v) ss << v << "_";
    return ss.str();
  }

  // for now it has to be the same was as function definition...
  float Operate(const std::vector<float>& data, Operation_t operation) {

    switch(operation) {
      
    case kOP_INVALID : 
      { 
	return data.front(); 
      }
    case kOP_EQE : 
      { 
	float res = CombinedFit_EnuQE_ryan(data[0],data[1],kNUE);
	res *= 1000.0;
	//std::cout << "EQE: (" << data[0] << "," << data[1] << ") = " << res << std::endl;
	return res;
      }

    case kOP_Q2 : 
      { 
	return CombinedFit_Q2_ryan(data[0],data[1],kNUE);
      }
    default:
      {
	throw sperr("Unknown operation observed");
      }
    }
    return kINVALID_FLOAT;
  }


}

#endif 
