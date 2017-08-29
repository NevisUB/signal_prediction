#ifndef __MODELNUECCQE_CXX__
#define __MODELNUECCQE_CXX__

#include "ModelNueCCQE.h"

namespace sp {

	bool ModelNueCCQE::Valid() {

		if (NuType == 3 or NuType == 4) 
			return true;

		return false;
	}

	float ModelNueCCQE::Operate(std::string reco_or_truth){

		if(reco_or_truth == "reco"){
			return RecoEnuQE*1000;		

		}else if (reco_or_truth == "truth"){

			//return NuMomT*1000;
			return NuMomT*1000;

		}
		return -99;
	}

}

#endif
