#ifndef __COMBINEDFUNCTIONS_H__
#define __COMBINEDFUNCTIONS_H__

#include "Base/SPTypes.h"
#include "CombinedTypes.h"

#include <string>
#include <vector>

namespace sp {

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

  //
  // Functions ported from CombinedFit_MakeStackHistograms.F
  // https://goo.gl/3oWC45
  //
  
  StackedBkgdType_t StackHistoBkgd(const bool Event_is_dirt,
				   const bool Event_is_pi0, 
				   const NuanceType_t evwt, 
				   const NuType_t nutype,   
				   const GEANT3Type_t parent_id);

  // Reweight of montecarlo based on number of pions, kaon, and background from fit?
  float StackedHistoBkdgWeight( const bool Event_is_dirt,   // tuple var 11
			        const bool Event_is_pi0,    // tuple var 10
				const NuanceType_t evwt,    // tuple var 3
				const NuType_t nutype,      // tuple var 4
				const GEANT3Type_t parent_id, // tuple var 8
				const float in_weight,
				const float Npi,
				const float Nka,
				const float Nbkg);
  
  
  
}

#endif
