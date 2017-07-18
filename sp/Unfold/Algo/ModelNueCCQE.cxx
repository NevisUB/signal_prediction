#ifndef __MODELNUECCQE_CXX__
#define __MODELNUECCQE_CXX__

#include "ModelNueCCQE.h"

namespace sp {

  bool ModelNueCCQE::Valid() {

    if (NuType == 3 or NuType == 4) 
      return true;

    return false;
  }

}

#endif
