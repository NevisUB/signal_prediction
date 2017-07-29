#ifndef __UNFOLDALGODAGNOSTINI_CXX__
#define __UNFOLDALGODAGNOSTINI_CXX__

#include "UnfoldAlgoDAgnostini.h"

namespace sp {


	void UnfoldAlgoDAgnostini::Unfold(){
		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Beginning unfolding with Algorithm: "<<name<<std::endl;

		TMatrixT<double> dudd(n_t, n_r);
		TVectorT<double> u_last(n_t);

		std::vector<std::vector<std::vector<double>>> dudA_last(n_t, std::vector<std::vector<double>>(n_r, std::vector<double>(n_t,0)));
		std::vector<std::vector<std::vector<double>>> covA(n_r, std::vector<std::vector<double>>(n_t, std::vector<double>(n_r,0)));//thre indicies, as we have ignored inter truth correlations i think


		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Calculating covariance on Response. MUST CHECK WARNING."<<std::endl;
		//This is necessary to propagate the errors on A through..
		for(int r=0;r<n_r; r++){
			for(int a=0;a<n_t; a++){
				for(int s=0;s<n_r; s++){
					if(r==s){
						covA.at(r).at(a).at(s) = 1.0/t(a) *A(r,a)*(1-A(r,a)); 
					}else{
						covA.at(r).at(a).at(s) = - 1.0/t(a)*A(r,a)*A(s,a); 
					}

				}
			}
		}






		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Setting initial guess."<<std::endl;
		u.ResizeTo(n_t);
		//what is the initial guess for bayes theorem;
		u=t;

		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Truth is : ";
		double ntruth =0;
		for(int a=0; a<n_t; a++){
			ntruth+=t(a);
			std::cout<<t(a)<<" ";
		}
		std::cout<<std::endl;

		u_last = u;
		for(int a=0; a<n_t; a++){
			// UNCOMMENT  this for flat prior 
			//u(a) = ntruth/((double)n_t); 
		}

	

		for(int k=0; k<regularization; k++){
			std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || On Iteration "<<k<<" of "<<regularization<<std::endl;


			TMatrixT<double> Munfold(n_t, n_r);	

			for(int a=0; a< n_t; a++){ //loop over the u elemements

				//Calculate "unfolding matrix"
				for(int i=0; i< n_r; i++){
					double de =0;

					for(int b=0;b<n_t; b++){
						de += A(i,b)*u(b);
					}

					Munfold(a,i) += A(i,a)*u(a)/(ep(a)*de) ;

				
					if(de==0){
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || ERROR denominator is: "<<de <<" leading to NAN's"<<std::endl;
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || ERROR This means the current reconstructed variable in this iteration is 0 for bin "<<i<<std::endl;
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || ERROR suggestion is to merge the last two reco bins, and repeat"<<std::endl;
						exit(EXIT_FAILURE);
					}
				}

			}//end of u loop

			//Unfold!
			u = Munfold*d;

			/***************************************
			 *	For each iteration, need to keep track of covariance and error matrix better
			 * ************************************/
			TMatrixT<double> dudd_cur(n_t, n_r);
			dudd_cur = Munfold;

			std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Begininning calculation of covariance matrix due to measurement uncertainty D."<<std::endl;
			if(k !=0){
				for(int a=0; a<n_t; a++){
					for(int i =0; i<n_r; i++){
						//already have munfold in

						if(u_last(a)==0){
							std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || WARNING: u_last(a) is zero, leading to infinity in dudd_cur: "<<dudd_cur(a,i)<<" @ a = "<<a<<". Previous dudd(a,i)= "<<dudd(a,i)<<" and u(i): "<<u(a)<<std::endl;
						}	


							dudd_cur(a,i) += u(a)/u_last(a)*dudd(a,i);


						for(int j=0; j< n_r; j++){
							for(int b=0; b< n_t; b++){


								dudd_cur(a,i)-= d(j)*ep(b)/u_last(b)*Munfold(a,j)*Munfold(b,j)*dudd(b,i);


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
			tr.Transpose(dudd_cur);

			for(int a =0; a<n_t; a++){
				for(int i=0; i<n_r; i++){
					//		std::cout<<"M "<<a<<" "<<i<<" "<<dudd_cur(a,i)<<" "<<D(a,i)<<" "<<tt(a,i)<<std::endl;

				}
			}


			U.Zero();
			U = dudd_cur*D*tr;




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



			std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Begininning calculation of bias!"<<std::endl;
			// And bias! This may not be taking it into account properly as of yet. Using bias approximation from cowans book. should be sufficient. 
			for(int a =0; a<n_t;a++){

				b(a) = 0;
				for(int j=0; j<n_r; j++){
					double vj = 0;
					for(int b=0; b<n_t; b++){
						vj+= A(j,b)*u(b);
					}
					b(a) += dudd_cur(a,j)*(vj-d(j));
					//b(a) += Munfold(a,j)*(vj-d(j));
				}
			}




			dudA_last = dudA;
			dudd = dudd_cur;
			u_last=u;



		}//end of iterative loop


	}

}

#endif
