#ifndef __TIKONOVSVD_CXX__
#define __TIKONOVSVD_CXX__

#include "TikonovSVD.h"
#include <Eigen/SVD>
#include <iostream>

namespace sp {
  
  void TikhonovSVD::Initialize(const Eigen::MatrixXf& A) {

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

    for(size_t c_i = 0; c_i < c_sz; ++c_i) diag_v[c_i] = -2.0 + _epsilon;
    
    diag_v[0]      = -1.0 + _epsilon;
    diag_v[c_sz-1] = -1.0 + _epsilon;

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


  Eigen::VectorXf TikhonovSVD::Unfold(const Eigen::VectorXf& b, size_t k) {
    
    _d = _U.transpose() * b;

    if (k==kINVALID_SIZE) k = b.size() - 1;
    
    auto tau = _s(k) * _s(k);
    
    auto top = _d.cwiseProduct(_s);
    auto bot = _s.cwiseProduct(_s) + tau*Eigen::VectorXf::Ones(_s.size());
    auto z = top.cwiseProduct(bot.cwiseInverse());
    
    return _C_inv * ( _V * z);
  }
  
}

#endif
