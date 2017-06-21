#ifndef BOONE_FNC_H_
#define BOONE_FNC_H_

#include <vector>
#include <string>

/**************************
 * DEFINE some  MiniBooNE codes
 **************************/

#define NOTNU 0
#define NUMU 1
#define NUMUBAR 2
#define NUE 3
#define NUEBAR 4

/**************************
 * DEFINE some  nuance codes
 * **************************/

#define NUA_CCQE 1
#define NUA_NCQE_PI 2
#define NUA_CCPIP_P 3
#define NUA_CCPIP_N 5
#define NUA_NCPI0_P 6
#define NUA_NCPI0_N 8
#define NUA_NCPIP 7
#define NUA_NCPIM 9
#define NUA_NCPI0_COH 96
#define NUA_CCPI0_COH 97

#define NUA_CCDELTA 93
#define NUA_NCDELTA 94

#define NUABAR_CCPIP_P 10
#define NUABAR_CCPIP_N 12
#define NUABAR_NCPI0_P 13
#define NUABAR_NCPI0_N 15

/**************************
 * DEFINE some  IBKGD codes
 **************************/
#define BCCQE 0
#define BCCPIP 1
#define BNCPI0 2
#define BCHPI0 3
#define BDELTA 4

std::string stringCalcIBKGD(int nu_type, int nuance_type ){

	if( nu_type == NUE || nu_type == NUMU){

		switch(nuance_type){
			case NUA_CCQE:
				if(nu_type == NUE){
					return "CCQE_nue";
				}else
				{
					return "CCQE_numu";

				}
			case NUA_NCQE_PI:
				return "NCQE_PI";
			case NUA_CCPIP_P:
			case NUA_CCPIP_N:
				return "CCPIP";
			case NUA_NCPI0_P:
			case NUA_NCPI0_N:
				return "NCPI0";
			case NUA_NCPI0_COH:
			case NUA_CCPI0_COH:
				return "CHPI0";
			case NUA_CCDELTA:
			case NUA_NCDELTA:
				return "DELTA";
			default:
				//std::cout<<nuance_type<<std::endl;
				return "Other";
		}	




	}else if (nu_type ==  NUEBAR || nu_type == NUMUBAR  ){



		switch(nuance_type){
			case NUA_CCQE:
				if(nu_type == NUEBAR){
					return "CCQE_nue";
				}else
				{
					return "CCQE_numu";

				}
			case NUA_NCQE_PI:
				return "NCQE_PI";
			case NUABAR_CCPIP_P:
			case NUABAR_CCPIP_N:
				return "CCPIP";
			case NUABAR_NCPI0_P:
			case NUABAR_NCPI0_N:
				return "NCPI0";
			case NUA_NCPI0_COH:
			case NUA_CCPI0_COH:
				return "CHPI0";
			case NUA_CCDELTA:
			case NUA_NCDELTA:
				return "DELTA";
			default:
				//std::cout<<nuance_type<<std::endl;
				return "Other";
		}	







	}


	return "Other";
}







int calcIBKGD(int nu_type, int nuance_type ){
	int ans = -999;


	if( nu_type == NUE || nu_type == NUMU){

		switch(nuance_type){
			case NUA_CCQE:
				return BCCQE;
			case NUA_CCPIP_P:
			case NUA_CCPIP_N:
				return BCCPIP;
			case NUA_NCPI0_P:
			case NUA_NCPI0_N:
				return BNCPI0;
			case NUA_NCPI0_COH:
			case NUA_CCPI0_COH:
				return BCHPI0;
			case NUA_CCDELTA:
			case NUA_NCDELTA:
				return BDELTA;
			default:
				return -99;
		}	




	}else if (nu_type ==  NUEBAR || nu_type == NUMUBAR  ){



		switch(nuance_type){
			case NUA_CCQE:
				return BCCQE;
			case NUABAR_CCPIP_P:
			case NUABAR_CCPIP_N:
				return BCCPIP;
			case NUABAR_NCPI0_P:
			case NUABAR_NCPI0_N:
				return BNCPI0;
			case NUA_NCPI0_COH:
			case NUA_CCPI0_COH:
				return BCHPI0;
			case NUA_CCDELTA:
			case NUA_NCDELTA:
				return BDELTA;
			default:
				return -99;
		}	







	}


	return -99;
}





#endif
