#ifndef __COMBINEDUNFOLD_CXX__
#define __COMBINEDUNFOLD_CXX__
#include <iostream>
#include "CombinedUnfold.h"
#include <Eigen/SVD>

namespace sp {
  
  void TikhonovSVD::Unfold(const Eigen::MatrixXf& A) {

    //
    // Store A
    //
    _A = A;

    //
    // Initialize C
    //
    const size_t c_sz = _A.cols();

    Eigen::VectorXf diag_v(c_sz);
    Eigen::VectorXf sub_diag_v(c_sz-1);

    for(size_t c_i = 0; c_i < c_sz; ++c_i) diag_v[c_i] = -2.0 + _tau;
    
    diag_v[0]      = -1.0 + _tau;
    diag_v[c_sz-1] = -1.0 + _tau;

    for(size_t c_i = 0; c_i < c_sz-1; ++c_i) sub_diag_v[c_i] = 1.0;

    _C.resize(c_sz,c_sz);
    _C.setZero(c_sz,c_sz);
    
    _C.diagonal()     = diag_v;
    _C.diagonal< 1>() = sub_diag_v;
    _C.diagonal<-1>() = sub_diag_v;

    //
    // Invert C
    //

    _C_inv = _C.inverse();

    //
    // Calculate A* C_inv
    //
    _A_C_inv = _A * _C_inv;

    //
    // SVD
    //
    Eigen::BDCSVD<Eigen::MatrixXf> bdc(_A_C_inv, Eigen::ComputeFullU | Eigen::ComputeFullV);

    _U = bdc.matrixU();
    _V = bdc.matrixV();
    _s = bdc.singularValues();
  }
  
}

#endif
