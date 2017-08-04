#ifndef __UNFOLDALGOSVD_H__
#define __UNFOLDALGOSVD_H__

#include "Unfold/Core/UnfoldAlgoBase.h"
#include <TDecompSVD.h>

namespace sp {

  class UnfoldAlgoSVD : public UnfoldAlgoBase {
  public:
    UnfoldAlgoSVD();
    ~UnfoldAlgoSVD(){}
    
    TMatrixD C;
    TMatrixD inv_C;

    TVectorD UT_td;
    TVectorD s_taic;

    TVectorD z;
    TMatrixD Z;

    TVectorD w;
    TMatrixD W;

    TMatrixD dudd; //nt x nr

    void Unfold();
    
  protected:
    
    void _Initialize_();
      
  private:

    void rotate_rescale(TMatrixD& tilde_A, 
			const TVectorD& r, 
			const TMatrixD& Q, 
			const TMatrixD& A);
    
    void rotate_rescale(TVectorD& tilde_b,
			const TVectorD& r,
			const TMatrixD& Q,
			const TVectorD& b);
  };

}

#endif
