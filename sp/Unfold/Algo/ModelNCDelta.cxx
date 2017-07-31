#ifndef __MODELNCDELTA_CXX__
#define __MODELNCDELTA_CXX__

#include "ModelNCDelta.h"

namespace sp {

  bool ModelNCDelta::Valid() {

    if (NUANCEChan == 94)
      return true;

    return false;
  }

}

#endif
