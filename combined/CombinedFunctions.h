#ifndef __COMBINEDFUNCTIONS_H__
#define __COMBINEDFUNCTIONS_H__

#include "CombinedTypes.h"
#include <string>

namespace sp {

  

  std::string Bkgd2String(BkgdType_t bkgd);
  BkgdType_t String2Bkgd(const std::string& name);


  //
  // Functions ported from EventLevelRoutines.F
  // https://goo.gl/N24LXp
  //
  
  BkgdType_t CombinedFit_bkgd_type(NuanceType_t evwt, NuType_t inno);

  // The Energy formula used by Mike S. in the Run Plan
  double CombinedFit_nue_qe(const double E, const double UZ);

  //Ryan's calculation of EnuQE
  //Use OneTrack_E and OneTrack_UZ
  double CombinedFit_EnuQE_ryan(const double evis_MeV, const double costheta, const NuType_t lepton_type);
  
  //Use OneTrack_E and OneTrack_UZ to calculate Q2
  double CombinedFit_Q2_ryan(const double evis_MeV, const double costheta, const NuType_t lepton_type);
}

#endif
