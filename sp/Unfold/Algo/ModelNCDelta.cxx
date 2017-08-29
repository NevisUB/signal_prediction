#ifndef __MODELNCDELTA_CXX__
#define __MODELNCDELTA_CXX__

#include "ModelNCDelta.h"

namespace sp {

	bool ModelNCDelta::Valid() { 


		if (NUANCEChan == 94)
			return true;




		return false;
	}


	float ModelNCDelta::Operate(std::string reco_or_truth){

		if(reco_or_truth == "reco"){
			return RecoEnuQE*1000;		

		}else if (reco_or_truth == "truth"){

			float gamma_max=-999;
			for(int i=0; i< NFSP; i++){
				if(FSPType->at(i)==1){
					if(MomT->at(i)>gamma_max){
						gamma_max= MomT->at(i);

					}
				}

			}
			//return NuMomT*1000;
			return gamma_max*1000;

		}
		return -99;
	}

}
#endif
