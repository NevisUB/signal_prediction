#ifndef __SPUTIL_CXX__
#define __SPUTIL_CXX__

#include "SPUtil.h"
#include "SPErr.h"

namespace sp {

  std::vector<double> Subtract(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) throw sperr("Given vectors are different sizes!");
    std::vector<double> res(v1.size(),0.0);
    for(size_t i=0; i<v1.size(); ++i) res[i] = v1[i] - v2[i];
    return res;
  }


}

#endif 
