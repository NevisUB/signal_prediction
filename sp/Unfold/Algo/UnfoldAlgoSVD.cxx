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

    R.ResizeTo(n_t,n_t);
    T.ResizeTo(n_t,n_t);

    for(int j=0; j<n_t; j++){
      R(j,j) = std::max(1.0, std::sqrt( r(j) ));
      T(j,j) = std::max(1.0, std::sqrt( t(j) ));
    }
    
    SP_DEBUG() << "end" << std::endl;
  }
  

  void UnfoldAlgoSVD::Unfold(){
    SP_INFO() << "Unfold" << std::endl;

    TMatrixD tilde_A(n_t,n_t);
    TVectorD tilde_d(n_t);
    
    SP_DEBUG() << "Decomposing D " << std::endl;
    
    TDecompSVD svd_D(D);
    const auto& r_D = svd_D.GetSig();
    const auto& Q_D = svd_D.GetU();

    SP_DEBUG() << "Rotate & Rescale A"<<std::endl;	
    rotate_rescale(tilde_A, r_D, Q_D, A);
    
    SP_DEBUG() << "Rotate & Rescale d"<<std::endl;	
    rotate_rescale(tilde_d, r_D, Q_D, d);

    auto tilde_A_inv_C = tilde_A * inv_C;

    SP_DEBUG() << "Decomposing A.C^{-1}" << std::endl;	

    TDecompSVD svd_taic(tilde_A_inv_C);
    const auto& s_taic = svd_taic.GetSig();
    const auto& U_taic = svd_taic.GetU();
    const auto& V_taic = svd_taic.GetV();

    auto UT_taic = U_taic; 
    UT_taic.T();

    auto UT_td = UT_taic * tilde_d;

    //
    // Apply regularization 
    //

    assert(regularization < s_taic.GetNrows());
    double tau = s_taic(regularization) * s_taic(regularization);
    
    z.ResizeTo(n_t);
    Z.ResizeTo(n_t,n_t);
    
    for(int i=0; i<n_t; i++){
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
    u.ResizeTo(n_t);
    U.ResizeTo(n_t,n_t);

    u = inv_C * (V_taic * z);

    auto inv_CT  = inv_C;
    auto VT_taic = V_taic;
    
    inv_CT.T();
    VT_taic.T();

    U = inv_C * V_taic * Z * VT_taic * inv_CT;

    xtau.ResizeTo(n_t);
    Xtau.ResizeTo(n_t,n_t);

    for(int i=0; i<n_t; i++){
      xtau(i) = t(i) * u(i);
      
      for(int j=0; j<n_t; j++)
	Xtau(i,j) = t(i) * U(i,j) * t(j);
    }


    //
    // Write out
    //
    /*

    std::vector<double> x(nt,0.0);
    std::vector<double> y(nt,0.0);
    std::vector<double> y2(nt,0.0);

    for (int i=0; i<nt; i++) {
      x[i]  = i;
      y[i]  = fabs(d(i)) ;
      y2[i] = s(i);
    }



    TCanvas c;
    c.Divide(2,1);

    TPad *p1 = (TPad *)(c->cd(1)); 
    p1->SetLogy();

    TGraph gr1(nt, x, &(y[0]));
    gr1.SetTitle("|d_{i}|");

    TPad *p2 = (TPad *)(c->cd(2)); 
    p2->SetLogy();

    TGraph gr2(nt, x, &(y2[0]));
    gr2.SetTitle("Single Values");
    gr2.SetMinimum(1e-3);
    gr2.Draw("AB");
    */
  }

  void UnfoldAlgoSVD::rotate_rescale(TMatrixD& tilde_Ain, 
				     const TVectorD& rin,
				     const TMatrixD& Q,
				     const TMatrixD& Ain){
    
    auto QAin = Q * Ain;

    assert(QAin.GetNrows() == tilde_Ain.GetNrows());
    assert(QAin.GetNcols() == tilde_Ain.GetNcols());

    for(int i=0; i < Ain.GetNrows(); i++){
      for(int j=0; j < Ain.GetNrows(); j++){
	double denom = rin(i);
	if (denom == 0.0) {
	  SP_CRITICAL() << "zero denominator" << std::endl;
	  throw sperr();
	}
	double num =  QAin(i,j);

	tilde_Ain(i,j) = num / denom;
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
