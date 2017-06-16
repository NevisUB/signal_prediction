#ifndef __COMBINEDFUNCTIONS_H__
#define __COMBINEDFUNCTIONS_H__

#include "CombinedTypes.h"
#include <string>

namespace sp {


  std::string Bkgd2String(BkgdType_t bkgd);
  BkgdType_t String2Bkgd(const std::string& name);

  BkgdType_t CombinedFit_bkgd_type(NuanceType_t evwt, NuType_t inno);

  // The Energy formula used by Mike S. in the Run Plan
  double CombinedFit_nue_qe(const double E, const double UZ);
}

#endif
