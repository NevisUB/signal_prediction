#ifndef __UNFOLDALGODAGNOSTINI_CXX__
#define __UNFOLDALGODAGNOSTNI_CXX__

#include "UnfoldAlgoDAgnostini.h"

namespace sp {


	void UnfoldAlgoDAgnostini::Unfold(){

		D.ResizeTo(n_r,n_r);
		D.Zero();
		for(int i=0; i<n_r; i++){
			D(i,i)=d(i);
		}

		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Beginning unfolding with Algorithm: "<<name<<std::endl;

		u.ResizeTo(n_t);
		//what is the initial guess for bayes theorem;
		u=t;
		
		double ntruth =0;
		for(int i=0; i<n_r; i++){
			ntruth+=t(i);
		}
	

		for(int a=0; a<n_t; a++){
			//u(a) = ntruth/((double)n_t); 
		}

		TVectorT<double> ep(n_t);
		ep.Zero();
		for(int b=0; b<n_t; b++){
			for(int l=0; l<n_r; l++){
				ep(b) += A(l,b) ;
			}
		}


		TMatrixT<double> meas_error_prop(n_t, n_r);
		TVectorT<double> u_last(n_t);

		for(int k=0; k<regularization; k++){
			std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || On Iteration "<<k<<" of "<<regularization<<std::endl;


			TMatrixT<double> Munfold(n_t, n_r);	

			for(int a=0; a< n_t; a++){ //loop over the u elemements

				//Calculate "unfolding matrix"
				for(int i=0; i< n_r; i++){
					double de =0;

					for(int b=0;b<n_t; b++){
						de += (A)(i,b)*u(b);
					}

					Munfold(a,i) += A(i,a)*u(a)/(ep(a)*de) ;

					if(ep(a)==0){
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || WARNING efficiency is: "<<ep(a)<<" leading to NAN's"<<std::endl;
					}

					if(de==0){
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || WARNING denominator is: "<<de <<" leading to NAN's"<<std::endl;
					}
				}

			}//end of u loop

			//Unfold!
			u = Munfold*d;

			/***************************************
			 *	For each iteration, need to keep track of covariance and error matrix better
			 * ************************************/
			TMatrixT<double> meas_error_prop_tmp(n_t, n_r);
			meas_error_prop_tmp = Munfold;

			if(k !=0){
				for(int a=0; a<n_t; a++){
					for(int i =0; i<n_r; i++){
						//already have munfold in
						meas_error_prop_tmp(a,i) += u(a)/u_last(a)*meas_error_prop(a,i);
						for(int j=0; j< n_r; j++){
							for(int b=0; b< n_t; b++){


								meas_error_prop_tmp(a,i)-= d(j)*ep(b)/u_last(b)*Munfold(a,j)*Munfold(b,j)*meas_error_prop(b,i);


							}
						}


					}
				}


			}


			/***********************************
			 *	Then calculate U due to D propagation
			 * ********************************/
			U.Zero();
			for(int a=0; a<n_t; a++){
				for(int b=0; b<n_t; b++){
			
					for(int i=0; i<n_r; i++){
						for(int j=0; j<n_r; j++){
							U(a,b) += meas_error_prop(a,i)*D(i,j)*meas_error_prop(b,j);
						}
					}

				}
			}



			meas_error_prop = meas_error_prop_tmp;
			u_last=u;

		}//end of iterative loop


	}

}

#endif
