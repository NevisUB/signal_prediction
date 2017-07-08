#ifndef __UNFOLDUTIL_CXX__
#define __UNFOLDUTIL_CXX__

#include "UnfoldUtil.h"

namespace sp {

  Eigen::VectorXf to_vector_eigen(std::vector<float> vec) {
    return Eigen::Map<Eigen::VectorXf>(vec.data(), vec.size());
  }

  std::vector<float> to_vector_std(const Eigen::VectorXf vec) {
    return std::vector<float>(vec.data(), vec.data() + vec.size());
  }
}

#endif
