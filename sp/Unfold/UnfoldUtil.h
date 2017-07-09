#ifndef __UNFOLDUTIL_H__
#define __UNFOLDUTIL_H__

#include <Eigen/Dense>
#include <vector>

#include "TH1F.h"
#include "TH1D.h"

namespace sp {


  Eigen::VectorXf to_vector_eigen(const std::vector<float>& vec);
  Eigen::VectorXf to_vector_eigen(const TH1F& th);

  Eigen::VectorXd to_vector_eigen(const std::vector<double>& vec);
  Eigen::VectorXd to_vector_eigen(const TH1D& th);
  
  std::vector<float>  to_vector_std(const Eigen::VectorXf& vec);
  std::vector<float> to_vector_std(const TH1F& th);

  std::vector<double> to_vector_std(const Eigen::VectorXd& vec);
  std::vector<double> to_vector_std(const TH1D& th);
  
}

#endif
