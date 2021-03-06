#ifndef __TIKONOVSVD_H__
#define __TIKONOVSVD_H__

#include "Base/SPTypes.h"
#include <Eigen/Dense>

namespace sp {

  class TikhonovSVD {

  public:
    TikhonovSVD() : _epsilon(1e-4) {}
    ~TikhonovSVD(){}

    void Initialize(const Eigen::MatrixXf& A);
    Eigen::VectorXf Unfold(const Eigen::VectorXf& b, size_t k = kINVALID_SIZE);
    size_t Rank() const { return _s.size(); }

    float _epsilon;
    
  private:

    Eigen::MatrixXf _C;
    Eigen::MatrixXf _C_inv;
    Eigen::MatrixXf _A;
    Eigen::MatrixXf _A_C_inv;
    Eigen::MatrixXf _U;
    Eigen::MatrixXf _V;
    Eigen::VectorXf _s;
    Eigen::VectorXf _d;

  public:
        
    const Eigen::MatrixXf& A() const { return _A; }
    const Eigen::MatrixXf& C() const { return _C; }
    const Eigen::MatrixXf& C_inv() const { return _C_inv; }
    const Eigen::MatrixXf& A_C_inv() const { return _A_C_inv; }
    const Eigen::MatrixXf& U() const { return _U; }
    const Eigen::MatrixXf& V() const { return _V; }
    const Eigen::VectorXf& s() const { return _s; }
    const Eigen::VectorXf& d() const { return _d; }    

  };
  

}

#endif
