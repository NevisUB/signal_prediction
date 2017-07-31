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

    TVectorD z;
    TMatrixD Z;

    TVectorD xtau;
    TMatrixD Xtau;

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
