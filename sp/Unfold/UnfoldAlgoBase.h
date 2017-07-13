#ifndef __UNFOLDALGOBASE_H__
#define __UNFOLDALGOBASE_H__

#include "Response.h"

namespace sp {

  class UnfoldAlgoBase { 
  public:
    UnfoldAlgoBase() {}
    virtual ~UnfoldAlgoBase() {}

    //
    // Pure virtual mothods called by Manager
    // 
    virtual void Initialize();
    virtual void Unfold(const TMatrixD& A) = 0;

  };
  

}

#endif
