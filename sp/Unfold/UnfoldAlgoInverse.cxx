#ifndef __UNFOLDALGOINVERSE_CXX__
#define __UNFOLDALGOINVERSE_CXX__

#include "UnfoldAlgoInverse.h"

namespace sp {


	void UnfoldAlgoInverse::Unfold(){

		std::cout<<"sp::UnfoldAlgoInverse::Unfold || Beginning unfolding with Algorithm: "<<name<<std::endl;
		if(n_r != n_t){
			std::cout<<"sp::UnfoldAlgoInverse::Unfold || ERROR: Inverse is only possible if n_r == n_t, i.e square matrix!"<<std::endl;
			exit(EXIT_FAILURE);
		}	

		TMatrixD invA(n_t,n_r);
		double det = 0;

		std::cout<<"sp::UnfoldAlgoInverse::Unfold || Inverting Response Matrix."<<std::endl;
		TMatrixD tmpA = A;
		invA = tmpA.Invert(&det);
		std::cout<<"sp::UnfoldAlgoInverse::Unfold || Inverted, Determinant = "<<det<<std::endl;
		std::cout<<"sp::UnfoldAlgoInverse::Unfold || Unfolding."<<std::endl;
		u = invA*d;	

		std::cout<<"sp::UnfoldAlgoInverse::Unfold || Calculating Covariance Matrix"<<std::endl;
		TMatrixD TrInvA = invA;
		TrInvA.Transpose(invA);
		U = invA*D*TrInvA;

		std::cout<<"sp::UnfoldAlgoInverse::Unfold || Done."<<std::endl;


	}
}

#endif
