#ifndef __UNFOLDALGOINVERSE_H__
#define __UNFOLDALGOINVERSE_H__

#include "Unfold/Core/UnfoldAlgoBase.h"

namespace sp {

  class UnfoldAlgoInverse : public UnfoldAlgoBase {

  private:
  public:
    UnfoldAlgoInverse() : UnfoldAlgoBase() {
      name = "Inverse";
    }

    void Unfold();

  };


}

#endif
