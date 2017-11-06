#ifndef __MODELNCPI0_H__
#define __MODELNCPI0_H__

#include "Unfold/Core/SPModelBase.h"

namespace sp {

  class ModelNCpi0 : public SPModelBase {
    
  public:
    ModelNCpi0() 
      : SPModelBase() 
    {}
    ~ModelNCpi0() {}

    bool Valid();
    float Operate(std::string);

  };
}

#endif
