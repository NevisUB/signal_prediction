#ifndef __UNFOLDALGOSVD_CXX__
#define __UNFOLDALGOSVD_CXX__

#include "UnfoldAlgoSVD.h"
#include <assert.h> 

namespace sp {

  UnfoldAlgoSVD::UnfoldAlgoSVD() : UnfoldAlgoBase("SVD") {}

  void UnfoldAlgoSVD::_Initialize_() 
  {
    name = "SVD";

    SP_DEBUG() << "start n_t = " << n_t << std::endl;

    C.ResizeTo(n_t,n_t);
    inv_C.ResizeTo(n_t,n_t);

    double xi = 1e-5;

    for(int i=0; i<n_t;i++){
      for(int j=0; j<n_t;j++){
	if(i == j){
	  if(i == 0 or i == n_t-1){
	    C(i,j) = -1.0+xi;
	  }else{
	    C(i,j) = -2.0+xi;
	  }
	}else{
	  C(i,j)=0.0;
	}
      }
    }

    for(int i=0; i<n_t-1; i++) {
      C(i  ,i+1) = 1;
      C(i+1,i  ) = 1;
    }

    TDecompSVD svd_C(C);

    const auto& cs = svd_C.GetSig();
    TMatrixD c_inv_s(n_t,n_t);

    for(int i=0; i<n_t; i++){
      for(int j=0; j<n_t; j++){

	if(i==j) {
	  if (cs(i) == 0.0)  {
	    SP_CRITICAL() << "cs @ " << i << " is zero" << std::endl;
	    throw sperr();
	  }
	  c_inv_s(i,j) = 1.0 / cs(i);
	}
	else {
	  c_inv_s(i,j) = 0.0;
	}
      }
    }
    
    const auto& cU = svd_C.GetU();
    const auto& cV = svd_C.GetV();
    auto cUT = cU;
    cUT.T();

    inv_C = cV * c_inv_s * cUT;

    R.ResizeTo(n_r, n_r);
    T.ResizeTo(n_t, n_t);
    
    for(int j=0; j<n_t; j++)
      R(j,j) = std::max(1.0, std::sqrt( r(j) ));
    
    for(int a=0; a<n_t; a++)
      T(a,a) = std::max(1.0, std::sqrt( t(a) ));
    
    SP_DEBUG() << "end" << std::endl;
  }
  

  void UnfoldAlgoSVD::Unfold(){
    SP_INFO() << "Unfold" << std::endl;

    TMatrixD tilde_A(n_r,n_t);
    TVectorD tilde_d(n_r);
    
    SP_DEBUG() << "Decomposing D " << std::endl;
    
    TDecompSVD svd_D(D);
    const auto& r_D = svd_D.GetSig();
    const auto& Q_D = svd_D.GetU();

    SP_DEBUG() << "Return r_D ==> " << std::endl;
    if(this->logger().level() == msg::kDEBUG) 
      r_D.Print();

    SP_DEBUG() << "Return Q_D ==> " << std::endl;
    if(this->logger().level() == msg::kDEBUG) 
      Q_D.Print();

    SP_DEBUG() << "Rotate & Rescale A"<<std::endl;	
    rotate_rescale(tilde_A, r_D, Q_D, A);
    
    SP_DEBUG() << "Rotate & Rescale d"<<std::endl;	
    rotate_rescale(tilde_d, r_D, Q_D, d);

    auto tilde_A_inv_C = tilde_A * inv_C;
    
    SP_DEBUG() << "Decomposing A.C^{-1}" << std::endl;	

    TDecompSVD svd_taic(tilde_A_inv_C);
    s_taic.ResizeTo(n_r);
    s_taic = svd_taic.GetSig();
    const auto& U_taic = svd_taic.GetU();
    const auto& V_taic = svd_taic.GetV();

    auto UT_taic = U_taic; 
    UT_taic.T();
    
    UT_td.ResizeTo(n_t);
    UT_td = UT_taic * tilde_d;
    
    SP_DEBUG() << "di ==> " << std::endl;
    if (this->logger().level() == msg::kDEBUG)
      UT_td.Print();

    SP_DEBUG() << "Apply regularization" << std::endl;
    
    assert(regularization >= 0);
    assert(regularization < s_taic.GetNrows());
    double tau = s_taic(regularization) * s_taic(regularization);
    
    z.ResizeTo(n_r);
    Z.ResizeTo(n_r,n_r);
    
    for(int i=0; i<n_r; i++){
      double num = UT_td(i)  * s_taic(i);
      double s2  = s_taic(i) * s_taic(i);
      double den = s2 + tau;
      if (den==0.0) {
	SP_CRITICAL() << "calculation of z contains invalid entry @ i="<<i<<std::endl;
	throw sperr();
      }
      z(i)  =  num / den;
      Z(i,i) = s2 / (den * den);
    }

    // unfolded spectrum & covariance
    w.ResizeTo(n_t);
    W.ResizeTo(n_t,n_t);

    w = inv_C * (V_taic * z);

    auto inv_CT  = inv_C;
    auto VT_taic = V_taic;
    
    inv_CT.T();
    VT_taic.T();

    W = inv_C * V_taic * Z * VT_taic * inv_CT;
    
    u.ResizeTo(n_t);
    U.ResizeTo(n_t,n_t);
    
    for(int i=0; i<n_t; i++){
      u(i) = t(i) * w(i);
      
      for(int j=0; j<n_t; j++)
	U(i,j) = t(i) * W(i,j) * t(j);
    }

    SP_DEBUG() << "Unfolded u ==> " << std::endl;
    if(this->logger().level() == msg::kDEBUG)
      u.Print();

  }

  void UnfoldAlgoSVD::rotate_rescale(TMatrixD& tilde_Ain, 
				     const TVectorD& rin,
				     const TMatrixD& Q,
				     const TMatrixD& Ain){
    
    auto QAin = Q * Ain;

    assert(QAin.GetNrows() == tilde_Ain.GetNrows());
    assert(QAin.GetNcols() == tilde_Ain.GetNcols());

    for(int i=0; i < Ain.GetNrows(); i++){
      for(int a=0; a < Ain.GetNcols(); a++){
	double denom = rin(i);
	if (denom == 0.0) {
	  SP_CRITICAL() << "zero denominator" << std::endl;
	  throw sperr();
	}
	double num =  QAin(i,a);

	tilde_Ain(i,a) = num / denom;
      }
    }

  }

  void UnfoldAlgoSVD::rotate_rescale(TVectorD& tilde_b,
				     const TVectorD& rin,
				     const TMatrixD& Q,
				     const TVectorD& bin ){
	
    tilde_b = Q * bin;

    for(int i=0; i < rin.GetNrows(); i++){
      if (rin(i) == 0.0) {
	SP_CRITICAL() << "zero denominator" << std::endl;
	throw sperr();
      }
      
      tilde_b(i) *= 1.0 / rin(i);
    }
  }

}

#endif
