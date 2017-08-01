#ifndef __ALGORITHMUTIL_CXX__
#define __ALGORITHMUTIL_CXX__

#include "AlgorithmUtil.h"

namespace sp {


  //
  //
  //
  Eigen::MatrixXf Cov(const Eigen::VectorXf& vec) {
    Eigen::MatrixXf centered = vec.rowwise() - vec.colwise().mean();
    float denom = centered.rows() > 1  ? centered.rows() - 1 : 1;
    Eigen::MatrixXf cov = (centered.adjoint() * centered) / denom;
    return cov;
  }

  Eigen::MatrixXd Cov(const Eigen::VectorXd& vec) {
    Eigen::MatrixXd centered = vec.rowwise() - vec.colwise().mean();
    double denom = centered.rows() > 1  ? centered.rows() - 1 : 1;
    Eigen::MatrixXd cov = (centered.adjoint() * centered) / denom;
    return cov;
  }

  
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

  //
  //
  //
  Eigen::MatrixXf Invert(const Eigen::MatrixXf& mat){
    return mat.inverse();
  }
  Eigen::MatrixXd Invert(const Eigen::MatrixXd& mat) {
    return mat.inverse();
  }

  //
  //
  //
  Eigen::MatrixXf Transpose(const Eigen::MatrixXf& mat){
    return mat.transpose();
  }
  Eigen::MatrixXd Transpose(const Eigen::MatrixXd& mat) {
    return mat.transpose();
  }

  //
  //
  //
  Eigen::MatrixXf Multiply(const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2){
    return mat1 * mat2;
  }

  Eigen::MatrixXd Multiply(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2){
    return mat1 * mat2;
  }

  Eigen::MatrixXf Multiply(const Eigen::MatrixXf& mat1, const Eigen::VectorXf& vec1){
    return mat1 * vec1;
  }

  Eigen::MatrixXd Multiply(const Eigen::MatrixXd& mat1, const Eigen::VectorXd& vec1){
    return mat1 * vec1;
  }

  //
  //
  //
  Eigen::VectorXf CMultiply(const Eigen::VectorXf& vec1, const Eigen::VectorXf& vec2){
    return vec1.cwiseProduct(vec2);
  }

  Eigen::VectorXd CMultiply(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2){
    return vec1.cwiseProduct(vec2);
  }

  Eigen::VectorXf CInvert(const Eigen::VectorXf& vec1) {
    return vec1.cwiseInverse();
  }

  Eigen::VectorXd CInvert(const Eigen::VectorXd& vec1) {
    return vec1.cwiseInverse();
  }

  //
  //
  //
  Eigen::MatrixXf Add(const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2){
    return mat1 + mat2;
  }

  Eigen::MatrixXd Add(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2){
    return mat1 + mat2;
  }

  //
  //
  //
  Eigen::MatrixXf Subtract(const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2){
    return mat1 - mat2;
  }

  Eigen::MatrixXd Subtract(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2){
    return mat1 - mat2;
  }

  
}

#endif
