#ifndef __COMBINEDFUNCTIONS_H__
#define __COMBINEDFUNCTIONS_H__

#include "CombinedTypes.h"
#include "PyCombined.h"

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

  
  //
  // Functions ported from  PhysicsFunctions.F
  // https://goo.gl/ksYaoH

  unsigned Pi0Details(const int nfsp,
		      const std::vector<int>& ipfs,
		      const std::vector<float>& vrtx_x,
		      const std::vector<float>& vrtx_y,
		      const std::vector<float>& vrtx_z,
		      const std::vector<float>& pfsp_x,
		      const std::vector<float>& pfsp_y,
		      const std::vector<float>& pfsp_z,
		      const std::vector<float>& pfsp_t,
		      int& ldalitz,
		      std::vector<float>& ppi,
		      std::vector<float>& pgam_x,
		      std::vector<float>& pgam_y,
		      std::vector<float>& pgam_z,
		      std::vector<float>& pgam_t,
		      std::vector<float>& vtx);

  
  unsigned Pi0Details(const int nfsp,
		      const std::vector<int>& ipfs,
		      const std::vector<float>& vrtx_x,
		      const std::vector<float>& vrtx_y,
		      const std::vector<float>& vrtx_z,
		      const std::vector<float>& pfsp_x,
		      const std::vector<float>& pfsp_y,
		      const std::vector<float>& pfsp_z,
		      const std::vector<float>& pfsp_t);


  unsigned Pi0Details(const int nfsp,
		      PyObject* ipfs,
		      PyObject* vrtx_x,
		      PyObject* vrtx_y,
		      PyObject* vrtx_z,
		      PyObject* pfsp_x,
		      PyObject* pfsp_y,
		      PyObject* pfsp_z,
		      PyObject* pfsp_t);

  

}

#endif
