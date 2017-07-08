#ifndef __UNFOLDUTIL_H__
#define __UNFOLDUTIL_H__

#include <Eigen/Dense>
#include <vector>

namespace sp {

  // why forced to make a copy?
  Eigen::VectorXf to_vector_eigen(std::vector<float> vec);

  std::vector<float> to_vector_std(const Eigen::VectorXf vec);

  
}

#endif
