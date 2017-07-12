#ifndef __UNFOLDUTIL_H__
#define __UNFOLDUTIL_H__

#include <vector>
#include "Combined/CombinedFunctions.h"
#include "Base/SPErr.h"
#include "UnfoldTypes.h"

namespace sp {

  std::vector<double> Subtract(const std::vector<double>& v1, const std::vector<double>& v2);

  std::string concatenate(const std::vector<std::string>& string_v);

  std::string to_name(Operation_t operation,
		      const std::vector<std::string>& prefix_v,
		      const std::vector<double>& bin_lo_v,
		      const std::vector<double>& bin_hi_v);

  // for now it has to be the same was as function definition...
  float Operate(const std::vector<float>& data, Operation_t operation);
  
}

#endif 
