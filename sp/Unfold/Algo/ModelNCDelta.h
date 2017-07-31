#ifndef __MODELNCDELTA_H__
#define __MODELNCDELTA_H__

#include "Unfold/Core/SPModelBase.h"

namespace sp {

  class ModelNCDelta : public SPModelBase {
    
  public:
    ModelNCDelta() 
      : SPModelBase() 
    {}
    ~ModelNCDelta() {}

    bool Valid();
  };
}

#endif
