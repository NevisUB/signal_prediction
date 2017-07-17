#ifndef __UNFOLDALGODAGNOSTINI_CXX_
#define __UNFOLDALGODAGNOSTNI_CXX__

#include "UnfoldAlgoDAgnostini.h"

namespace sp {


	void UnfoldAlgoDAgnostini::Unfold(const TVectorD * data, const TMatrixD * Din){
		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || Beginning unfolding with Algorithm: "<<name<<std::endl;

		u.ResizeTo(n_t);
		u=t;

		U.ResizeTo(n_t,n_t);
		U.Zero();

		d.ResizeTo(n_r);
		d = *data;

		for(int k=0; k<regularization; k++){
		std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || On Iteration "<<k<<" of "<<regularization<<std::endl;

			for(int a=0; a< n_t; a++){ //loop over the u elemements

				double tml =0;
				double ep = 0;

				//Calculate epsilon
				for(int i=0; i<n_t; i++){
					ep += A(i,a) ;
				}
				//Calculate "unfolding matrix"
				for(int i=0; i< n_r; i++){
					double de =0;

					for(int b=0;b<n_t; b++){
						de += (A)(i,b)*u(b);
					}
					tml += A(i,a)*u(a)*d(i)/(ep*de) ;

					if(ep==0){
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || WARNING efficiency is: "<<ep<<" leading to NAN's"<<std::endl;
					}
	
					if(de==0){
						std::cout<<"sp::UnfoldAlgoDAgnostini::Unfold || WARNING denominator is: "<<de <<" leading to NAN's"<<std::endl;
					}
				}
				
				u(a) = tml;
			}//end of u loop

		}//end of iterative loop

	}

	void UnfoldAlgoDAgnostini::Unfold(const TVectorD * data){
		d.Zero();
		d = *data;

		this->Unfold();

	}

	void UnfoldAlgoDAgnostini::Unfold(){
		this->Unfold(&d,&D);



	}

}

#endif
