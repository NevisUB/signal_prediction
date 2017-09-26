#ifndef __UNFOLDALGOSVD_CXX__
#define __UNFOLDALGOSVD_CXX__

#include "UnfoldAlgoSVD.h"
#include <assert.h> 

namespace sp {

	UnfoldAlgoSVD::UnfoldAlgoSVD() : UnfoldAlgoBase("SVD") {
	direct_regularization = false;
	}

	void UnfoldAlgoSVD::_Initialize_() 
	{
		name = "SVD";

		SP_DEBUG() << "Starting SVD initilization, n_t  = " << n_t <<" and n_r = "<<n_r<<std::endl;

		assert(n_r>=n_t);

		//resizing up some matricies
		dudd.ResizeTo(n_t,n_r);
		C.ResizeTo(n_t,n_t);
		inv_C.ResizeTo(n_t,n_t);


		//Making the curvature matrix C, xi is a small off diagonal to invert
		double xi = 1e-6;
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


		//This is just to invert C, via SVD
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

		//And save the value of the inverted C
		inv_C = cV * c_inv_s * cUT;



		//R and T are the uncertanties on R and T in MC, unsure if they are used.
		R.ResizeTo(n_r, n_r);
		T.ResizeTo(n_t, n_t);

		for(int j=0; j<n_r; j++){
			R(j,j) = std::max(1.0, std::sqrt( r(j) ));
		}

		for(int a=0; a<n_t; a++){
			T(a,a) = std::max(1.0, std::sqrt( t(a) ));
		}


		SP_DEBUG() << "Finished intializing SVD algorithm" << std::endl;
	}


	void UnfoldAlgoSVD::Unfold(){
		SP_INFO() << "Unfold() Begin" << std::endl;

		//The tilde variables are to keep the rotated and rescaled
		TMatrixD tilde_A(n_r,n_t);
		TVectorD tilde_d(n_r);

		//We decompose D (the covariance matrix of the observed Data) for the rescaling
		SP_DEBUG() << "Decomposing D ==> " << std::endl;
		
		//This should be true
		for(int i=0; i<n_r; i++){
			assert(D(i,i)>0);

		}


		TDecompSVD svd_D(D);
		const auto& r_D = svd_D.GetSig();//get significant values
		const auto& Q_D = svd_D.GetU();

		SP_DEBUG() << "Return r_D ==> " << std::endl;
		if(this->logger().level() == msg::kDEBUG) 
			r_D.Print();

		SP_DEBUG() << "Return Q_D ==> " << std::endl;
		if(this->logger().level() == msg::kDEBUG) 
			Q_D.Print();

		SP_DEBUG() << "Rotate & Rescale A"<<std::endl;	
		rotate_rescale(tilde_A, r_D, Q_D, N);

		SP_DEBUG() << "Rotate & Rescale d"<<std::endl;	
		rotate_rescale(tilde_d, r_D, Q_D, d);

		auto tilde_A_inv_C = tilde_A * inv_C;

		SP_DEBUG() << "Decomposing A.C^{-1}" << std::endl;	

		TDecompSVD svd_taic(tilde_A_inv_C);



		//After decomposition  we have A(n_r, n_t) = U(n_r,n_r) S(n_r, n_t) Vt (n_t,n_t)
		//But as S is diagonal we just keep its diagonal which is the ssmaller one, n_t as n_r>=n_t 
		s_taic.ResizeTo(n_t);
		s_taic = svd_taic.GetSig();

		const auto& U_taic = svd_taic.GetU();//nr x nr
		const auto& V_taic = svd_taic.GetV();//nt x nt

		auto UT_taic = U_taic; 
		UT_taic.T();

		UT_td.ResizeTo(n_r);
		//SP_DEBUG() << "3 " <<n_t<<" "<<n_r<<" UT: "<<UT_taic.GetNcols()<<" "<<UT_taic.GetNrows()<< " VT: "<<V_taic.GetNcols()<<" "<<V_taic.GetNrows()<<std::endl;	
		UT_td = UT_taic * tilde_d;

		SP_DEBUG() << "di ==> " << std::endl;
		if (this->logger().level() == msg::kDEBUG)
			UT_td.Print();

		SP_DEBUG() << "Apply regularization" << std::endl;

		assert(regularization >= 0);
		//assert(regularization < s_taic.GetNrows());
		std::cout<<"STAIC: "<<s_taic(0)<<" "<<s_taic(n_t-1)<<std::endl;
		
		//Either pass in a integer or a direct tau
		double tau;
		if(direct_regularization){
			tau = regularization;
		}else{
			tau = s_taic(regularization) * s_taic(regularization);
		}

		// z is the final unknown we are solving for, a function of w which is itself t (all n_t)
		z.ResizeTo(n_t);
		Z.ResizeTo(n_t,n_t);


		for(int a=0; a<n_t; a++){
			double num = UT_td(a)  * s_taic(a);
			double s2  = s_taic(a) * s_taic(a);
			double den = s2 + tau;
			if (den==0.0) {
				SP_CRITICAL() << "calculation of z contains invalid entry @ a="<<a<<std::endl;
				throw sperr();
			}
			z(a)  =  num / den;
			Z(a,a) = s2 / (den * den);
		}

		// unfolded spectrum & covariance
		w.ResizeTo(n_t);
		W.ResizeTo(n_t,n_t);

		w = inv_C * V_taic * z;

		//transpose a bit
		auto inv_CT  = inv_C;
		auto VT_taic = V_taic;
		inv_CT.T();
		VT_taic.T();


		W = inv_C * V_taic * Z * VT_taic * inv_CT;

		u.ResizeTo(n_t);
		U.ResizeTo(n_t,n_t);
		u.Zero();
		U.Zero();


		for(int a=0; a<n_t; a++){
			u(a) = t(a) * w(a);



			for(int b=0; b<n_t; b++){
				U(a,b) = t(a) * W(a,b) * t(b);

				//if(a==b)std::cout<<"LBIN: "<<a<<" t "<<t(a)<<" w "<<w(a)<<" b "<<b<<" t "<<t(b)<<" w "<<w(b)<<" W "<<W(a,b)<<" "<<U(a,b)<<" Err :" <<sqrt(U(a,b))<<std::endl;
			}
		}
	


		SP_DEBUG() << "Unfolded u ==> " << std::endl;
		if(this->logger().level() == msg::kDEBUG)
			u.Print();


		//Marks quick tests of matrix formulism.
		//Note S is not necessarily matrix of signular values
		TMatrixD S(n_t, n_r);
		S.Zero();

		for(int i=0; i< n_r; i++){
			for(int a=0; a<n_t; a++){
				if(a==i){
					S(a,i) = s_taic(a)/(s_taic(a)*s_taic(a)+tau );
				}

			}
		}

		TVectorD w_other(n_t);
	
		//	nt x nt * nt_nt * n_t x n_r  * nr x nr * n_r	
		w_other = inv_C * V_taic* S * UT_taic *tilde_d;

		for(int a=0; a<n_t;a++){
			std::cout<<"TESTA1: "<<a<<" Old "<<w(a)<<" New: "<<w_other(a)<<" "<<w(a)-w_other(a)<<std::endl;
		}	


		SP_DEBUG()<<"Begininning calculation of bias, using derivative + taylor approx from Cowan. Calc on w and btilde, then propagate."<<std::endl;
	
		//Useful to define a n_t unit matrix here	
		TMatrixD unit(n_t,n_t);
		unit.Zero();
		for(int a=0;a<n_t;a++){
			unit(a,a)=1.0;
		} 


		//Error propagation matrix, dw/d(tilde_d), explicitly assumes that V_taic and S are NOT dependent on tilde_b, which isn't super true as tilde_A has scaled out by r_i
		TMatrixD dw_dtilded(n_t, n_r);
		dw_dtilded = inv_C * V_taic* S * UT_taic;


		for(int a =0; a<n_t;a++){
			b(a) = 0;

			for(int j=0; j<n_r; j++){
				double refold_j = 0;
				
				for(int g=0; g<n_t; g++){
					refold_j += tilde_A(j,g)*w(g);
				}

				//Bias is currently on  w(a)
				b(a) += dw_dtilded(a,j)*(refold_j-tilde_d(j));
			}
		}
		
		//and move from bias on w to bias on u by removing t-scaling	
		for(int a =0; a<n_t;a++){
			b(a) = b(a)*t(a);
		}



		SP_DEBUG()<<"Begininning calculation of Covariance on the bias, ignoring dn/dd variance."<<std::endl;
		
		TMatrixD CRmI(n_t, n_t);
		//CRmI = dw_dtilded*tilde_A-unit;
		for(int a=0; a<n_t; a++){
			for(int g=0; g< n_t; g++){
				CRmI(a,g)=0;

				for( int i=0; i< n_r; i++){
					CRmI(a,g) += t(a)*dw_dtilded(a,i)*tilde_A(i,g);
				}
			}
		}
		CRmI = CRmI-unit;
	

		TMatrixD CRmI_T = CRmI;
		CRmI_T.T();
		
		B.ResizeTo(n_t, n_t);
		B = CRmI*W*CRmI_T;
	
		for(int a=0;a<n_t; a++){
			std::cout<<"BIN: "<<a<<" B "<<B(a,a)<<" W "<<W(a,a)<<std::endl;

		}

		//Bit of comparing to wiener

		SP_DEBUG()<<"Finished Unfold() call"<<std::endl;

	}

	void UnfoldAlgoSVD::SetDirectRegularization(double reg){
		SP_DEBUG()<<"Setting Direct regularization parameter to "<<reg<<std::endl;
		direct_regularization = true;
		regularization = reg;
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

				double denom = sqrt( rin(i) );
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

			tilde_b(i) *= 1.0 / sqrt( rin(i));
		}
	}

}

#endif
