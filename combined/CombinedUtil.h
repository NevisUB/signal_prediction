#ifndef __COMBINEDUTIL_H__
#define __COMBINEDUTIL_H__

#include "CombinedTypes.h"
#include <string>

namespace sp {
  
  std::string Bkgd2String(const BkgdType_t bkgd);
  BkgdType_t String2Bkgd(const std::string& name);

  std::string GEANT32String(const GEANT3Type_t particle);
  GEANT3Type_t String2GEANT3(const std::string& name);
  
}

#endif
