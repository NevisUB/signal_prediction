#ifndef __MODELNUECCQE_H__
#define __MODELNUECCQE_H__

#include "Unfold/Core/SPModelBase.h"

namespace sp {

  class ModelNueCCQE : public SPModelBase {
    
  public:
    ModelNueCCQE() 
      : SPModelBase() 
    {}
    ~ModelNueCCQE() {}

    bool Valid();
  };
}

#endif
