#ifndef __ALGORITHMUTIL_H__
#define __ALGORITHMUTIL_H__

#include <Eigen/Dense>
#include <vector>

#include "TH1F.h"
#include "TH1D.h"

#include "TH2F.h"
#include "TH2D.h"


namespace sp {


  //
  // Useful functions
  //
  Eigen::MatrixXf Cov(const Eigen::VectorXf& vec);
  Eigen::MatrixXd Cov(const Eigen::VectorXd& vec);
  
  
  //
  // Conversion wrapper for python
  //
  Eigen::VectorXf to_vector_eigen(const std::vector<float>& vec);
  Eigen::VectorXf to_vector_eigen(const TH1F& th);

  Eigen::VectorXd to_vector_eigen(const std::vector<double>& vec);
  Eigen::VectorXd to_vector_eigen(const TH1D& th);
  
  std::vector<float>  to_vector_std(const Eigen::VectorXf& vec);
  std::vector<float> to_vector_std(const TH1F& th);

  std::vector<double> to_vector_std(const Eigen::VectorXd& vec);
  std::vector<double> to_vector_std(const TH1D& th);

  Eigen::MatrixXf to_mat_eigen(const TH2F& th);
  Eigen::MatrixXd to_mat_eigen(const TH2D& th);


  //
  // Arithmetic convenience wrappers for python
  //

  Eigen::MatrixXf Invert(const Eigen::MatrixXf& mat);
  Eigen::MatrixXd Invert(const Eigen::MatrixXd& mat);

  Eigen::MatrixXf Transpose(const Eigen::MatrixXf& mat);
  Eigen::MatrixXd Transpose(const Eigen::MatrixXd& mat);

  Eigen::MatrixXf Multiply(const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2);
  Eigen::MatrixXd Multiply(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);

  Eigen::MatrixXf Multiply(const Eigen::MatrixXf& mat1, const Eigen::VectorXf& vec1);
  Eigen::MatrixXd Multiply(const Eigen::MatrixXd& mat1, const Eigen::VectorXd& vec1);

  Eigen::VectorXf CMultiply(const Eigen::VectorXf& vec1, const Eigen::VectorXf& vec2);
  Eigen::VectorXd CMultiply(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2);

  Eigen::VectorXf CInvert(const Eigen::VectorXf& vec1);
  Eigen::VectorXd CInvert(const Eigen::VectorXd& vec1);

  Eigen::MatrixXf Add(const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2);
  Eigen::MatrixXd Add(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);
  
  Eigen::MatrixXf Subtract(const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2);
  Eigen::MatrixXd Subtract(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);

  
}



#endif
