#ifndef __MODELNCPI0_CXX__
#define __MODELNCPI0_CXX__

#include "ModelNCpi0.h"

namespace sp {

	bool ModelNCpi0::Valid() {

		if (NUANCEChan == 6 || NUANCEChan == 8 || NUANCEChan == 96 || NUANCEChan == 13 || NUANCEChan == 15) 
			return true;
		return false;
	}

	float ModelNCpi0::Operate(std::string reco_or_truth){

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
