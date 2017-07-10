#ifndef __UNFOLDUTIL_CXX__
#define __UNFOLDUTIL_CXX__

#include "UnfoldUtil.h"

namespace sp {

  //
  // Eigen::VectorXx
  //
  Eigen::VectorXf to_vector_eigen(const std::vector<float>& vec) {
    return Eigen::Map<const Eigen::VectorXf>(vec.data(), vec.size());
  }

  Eigen::VectorXd to_vector_eigen(std::vector<double>& vec) {
    return Eigen::Map<const Eigen::VectorXd>(vec.data(), vec.size());
  }

  Eigen::VectorXf to_vector_eigen(const TH1F& th) {
    return Eigen::Map<const Eigen::VectorXf>(th.GetArray() + 1, th.GetSize() - 2);
  }

  Eigen::VectorXd to_vector_eigen(const TH1D& th) {
    return Eigen::Map<const Eigen::VectorXd>(th.GetArray() + 1, th.GetSize() - 2);
  }


  //
  // std::vector
  //
  std::vector<float> to_vector_std(const Eigen::VectorXf& vec) {
    return std::vector<float>(vec.data(), vec.data() + vec.size());
  }

  std::vector<double> to_vector_std(const Eigen::VectorXd& vec) {
    return std::vector<double>(vec.data(), vec.data() + vec.size());
  }
  
  std::vector<float> to_vector_std(const TH1F& th) {
    return std::vector<float>(th.GetArray() + 1, (th.GetArray() + 1) + (th.GetSize() - 2) );
  }
  
  std::vector<double> to_vector_std(const TH1D& th) {
    return std::vector<double>(th.GetArray() + 1, (th.GetArray() + 1) + (th.GetSize() - 2) );
  }


  //
  // Eigen::MatrixXx
  // 

  Eigen::MatrixXf to_mat_eigen(const TH2F& th) {
    auto mat = Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>( th.GetArray(),
													th.GetNbinsY()+2,
													th.GetNbinsX()+2);
    return mat.block(1,1,th.GetNbinsY(),th.GetNbinsX());
  }
  
  Eigen::MatrixXd to_mat_eigen(const TH2D& th) {
    auto mat = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>( th.GetArray(),
													 th.GetNbinsY()+2,
													 th.GetNbinsX()+2);
    return mat.block(1,1,th.GetNbinsY(),th.GetNbinsX());
  }
  
}

#endif
