#ifndef __UNFOLDALGODAGNOSTINI_CXX__
#define __UNFOLDALGODAGNOSTINI_CXX__

#include "UnfoldAlgoDAgnostini.h"

namespace sp {


	void UnfoldAlgoDAgnostini::Unfold(){
		TMatrixT<double> meas_error_prop(n_t, n_r);
		TVectorT<double> u_last(n_t);

		std::vector<std::vector<std::vector<double>>> dudA_last(n_t, std::vector<std::vector<double>>(n_r, std::vector<double>(n_t,0)));

		std::vector<std::vector<std::vector<double>>> covA(n_r, std::vector<std::vector<double>>(n_t, std::vector<double>(n_r,0)));//thre indicies, as we have ignored inter truth correlations i think
		//This is necessary to propagate the errors on A through..
		for(int r=0;r<n_r; r++){
			for(int a=0;a<n_t; a++){
				for(int s=0;s<n_r; s++){
					if(r==s){
						covA.at(r).at(a).at(s) = 1.0/t(a) *A(a,r)*(1-A(a,r)); 
					}else{
						covA.at(r).at(a).at(s) = - 1.0/t(a)*A(a,r)*A(a,s); 
					}

				}
			}
		}




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

		u_last = u;
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

			std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Begininning calculation of covariance matrix due to measurement uncertainty D."<<std::endl;
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

			TMatrixT<double> tt(n_t,n_r);
			TMatrixT<double> tr(n_r,n_t);
			tr.Transpose(meas_error_prop_tmp);

			for(int a =0; a<n_t; a++){
				for(int i=0; i<n_r; i++){
					//		std::cout<<"M "<<a<<" "<<i<<" "<<meas_error_prop_tmp(a,i)<<" "<<D(a,i)<<" "<<tt(a,i)<<std::endl;

				}
			}


			U.Zero();
			U = meas_error_prop_tmp*D*tr;


			/*
			   for(int a=0; a<n_t; a++){
			   for(int b=0; b<n_t; b++){
			   for(int i=0; i<n_r; i++){
			   for(int j=0; j<n_r; j++){
			   U(a,b) += meas_error_prop(a,i)*D(i,j)*meas_error_prop(b,j);
			   }
			   }

			   }
			   }
			   */

			//What is f!!


			/***********************************
			 *	Then calculate U due to Underlying A
			 * ********************************/
			std::vector<std::vector<std::vector<double>>> dudA(n_t, std::vector<std::vector<double>>(n_r, std::vector<double>(n_t,0)));

			TVectorD f(n_r);
			f = A*u_last;

			std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Begininning calculation of covariance matrix due to unertainty in response matrix A."<<std::endl;

			for(int a =0; a<n_t; a++){
				for(int j =0; j<n_r; j++){
					for(int b=0; b<n_t; b++){
						double dab=0;
						if(a==b) dab = 1.0;

						double t1 =1.0/ep(a)*(u_last(a)*d(j)/f(j)-u(a))*dab;			
						double t2 = -u_last(b)*d(j)/f(j)*Munfold(a,j);
						double t3 = u(a)/u_last(a)*dudA_last.at(a).at(j).at(b);

						double t4 =0;
						for(int l=0; l<n_r; l++){
							for(int g=0;g<n_t; g++){
								t4+=d(l)*Munfold(a,l)*Munfold(g,l)*dudA_last.at(g).at(j).at(b);

							}
						}
						t4 = -t4*ep(a)/u_last(a);



						dudA.at(a).at(j).at(b) = t1+t2+t3+t4;

					}
				}
			}


			//these two are outer indicied
			for(int a=0; a<n_t; a++){			
				for(int b=0; b<n_t; b++){			
					UA(a,b) = 0;


					for(int j=0; j<n_r; j++){
						for(int i=0; i<n_r; i++){
							for(int g=0; g<n_t; g++){
								UA(a,b)+= dudA.at(a).at(j).at(g)*covA.at(j).at(g).at(i)*dudA.at(b).at(i).at(g);
							}
						}


					}//end of interior loop
				}
			}



				dudA_last = dudA;
				meas_error_prop = meas_error_prop_tmp;
				u_last=u;



			}//end of iterative loop


		}

	}

#endif
