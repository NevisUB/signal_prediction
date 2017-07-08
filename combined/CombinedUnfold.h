#ifndef __COMBINEDUNFOLD_H__
#define __COMBINEDUNFOLD_H__

#include "CombinedTypes.h"
#include "Eigen/Dense"

namespace sp {

  class TikhonovSVD {
  public:
    TikhonovSVD() : _tau(1e-4) {}
    ~TikhonovSVD(){}

    void Unfold(const Eigen::MatrixXf& A);
    
    const Eigen::MatrixXf& A() const { return _A; }
    const Eigen::MatrixXf& C() const { return _C; }
    const Eigen::MatrixXf& C_inv() const { return _C_inv; }
    const Eigen::MatrixXf& A_C_inv() const { return _A_C_inv; }
    const Eigen::MatrixXf& U() const { return _U; }
    const Eigen::MatrixXf& V() const { return _V; }
    const Eigen::VectorXf& s() const { return _s; }

    float _tau;
    
  private:
    Eigen::MatrixXf _C;
    Eigen::MatrixXf _C_inv;
    Eigen::MatrixXf _A;
    Eigen::MatrixXf _A_C_inv;
    Eigen::MatrixXf _U;
    Eigen::MatrixXf _V;
    Eigen::VectorXf _s;
    
  };
  

}

#endif
