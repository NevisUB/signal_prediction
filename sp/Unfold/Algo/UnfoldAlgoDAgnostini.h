#ifndef __UNFOLDALGODAGNOSTINI_H__
#define __UNFOLDALGODAGNOSTINI_H__

#include "Unfold/Core/UnfoldAlgoBase.h"

namespace sp {

  class UnfoldAlgoDAgnostini : public UnfoldAlgoBase {

  public:
    UnfoldAlgoDAgnostini() : UnfoldAlgoBase("UnfoldAlgoDAgostini") {}
    ~UnfoldAlgoDAgnostini() {}

    void Unfold();

  };


}

#endif
