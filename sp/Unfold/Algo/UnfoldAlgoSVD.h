#ifndef __UNFOLDALGOSVD_H__
#define __UNFOLDALGOSVD_H__

#include "Unfold/Core/UnfoldAlgoBase.h"
#include <TDecompSVD.h>
namespace sp {

  class UnfoldAlgoSVD : public UnfoldAlgoBase {

  private:
  public:
		
    UnfoldAlgoSVD();
		
    TMatrixD C;
    TMatrixD invC;

    double xi;



    void Unfold();

    int rotateRescale(TMatrixT<double> * tilde_A, TVectorT<double> *r, TMatrixT<double> * Q, TMatrixT<double> * A );
    int rotateRescale(TVectorT<double> * tilde_b, TVectorT<double> *r, TMatrixT<double> * Q, TVectorT<double> * b );


  };


}

#endif