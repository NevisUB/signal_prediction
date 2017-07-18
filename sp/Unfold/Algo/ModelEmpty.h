#ifndef __MODELEMPTY_H__
#define __MODELEMPTY_H__

#include "Unfold/Core/SPModelBase.h"

namespace sp {

  class ModelEmpty : public SPModelBase {
    
  public:
    ModelEmpty() 
      : SPModelBase() 
    {}
    ~ModelEmpty() {}

    bool Valid() { return true; }
  };
}

#endif
