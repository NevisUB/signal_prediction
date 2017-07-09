#ifndef __UNFOLDUTIL_CXX__
#define __UNFOLDUTIL_CXX__

#include "UnfoldUtil.h"

namespace sp {

  class LoadUnfold {
  public:
    LoadUnfold(){}
    ~LoadUnfold(){}
  };

  
  Eigen::VectorXf to_vector_eigen(const std::vector<float>& vec) {
    return Eigen::Map<const Eigen::VectorXf>(vec.data(), vec.size());
  }

  Eigen::VectorXd to_vector_eigen(std::vector<double>& vec) {
    return Eigen::Map<const Eigen::VectorXd>(vec.data(), vec.size());
  }

  Eigen::VectorXf to_vector_eigen(const TH1F& th) {
    return Eigen::Map<const Eigen::VectorXf>(th.GetArray(), th.GetSize());
  }

  Eigen::VectorXd to_vector_eigen(const TH1D& th) {
    return Eigen::Map<const Eigen::VectorXd>(th.GetArray(), th.GetSize());
  }
  
  std::vector<float> to_vector_std(const Eigen::VectorXf& vec) {
    return std::vector<float>(vec.data(), vec.data() + vec.size());
  }

  std::vector<double> to_vector_std(const Eigen::VectorXd& vec) {
    return std::vector<double>(vec.data(), vec.data() + vec.size());
  }
  
  
}

#endif
