#ifndef __COMBINEDFUNCTIONS_CXX__
#define __COMBINEDFUNCTIONS_CXX__

#include "CombinedFunctions.h"

#include <iostream>
#include <cmath>

namespace sp { 

  BkgdType_t CombinedFit_bkgd_type(const NuanceType_t evwt, const NuType_t inno) {

    BkgdType_t ibkgd = kBINVALID;
    
    if ((inno == kNUMU) or (inno == kNUE)) {
      if (evwt == kCCQE) ibkgd = kBCCQE;
      else if (evwt == kCC1p1pi  or evwt == kCC1n1pi)  ibkgd = kBCCPIP;
      else if (evwt == kNC1p1pi0 or evwt == kNC1n1pi0) ibkgd = kBNCPI0;
      else if (evwt == kNC1pi01a or evwt == kCC1pi1a)  ibkgd = kBCHPI0;
      else if (evwt == kCC1pNg   or evwt == kNC1pNg)   ibkgd = kBDELTA;
    }
    else  {
      if ( evwt == kCCQE ) ibkgd = kBCCQE;
      else if (evwt == kCCbar1p1pi  or evwt == kCCbar1n1pi)  ibkgd = kBCCPIP;
      else if (evwt == kNCbar1p1pi0 or evwt == kNCbar1n1pi0) ibkgd = kBNCPI0;
      else if (evwt == kNC1pi01a    or evwt == kCC1pi1a)     ibkgd = kBCHPI0;
      else if (evwt == kCC1pNg      or evwt == kNC1pNg)      ibkgd = kBDELTA;
    }
    
    return ibkgd;
  }


  double CombinedFit_nue_qe(const double E, const double UZ) {

    double MPROT = 0.938;
    double MEL = 0.0005;
    double PEL = -0.22281e-2 + 0.10068e-2*E;
    double EEL = std::sqrt(PEL*PEL+MEL*MEL);
    double UZCORR = 0.80111e-2 + 0.99295*UZ;
    double ENUQE = 0.5 * (2.*MPROT*EEL - MEL*MEL) / (MPROT - EEL+PEL*UZCORR);

    double res = 0.73420E-01 + 0.92557*ENUQE;

    return res;
  }

  double CombinedFit_EnuQE_ryan(const double evis_MeV, const double costheta, const NuType_t lepton_type) {

    double Me = 0.00051099906;
    double Mm = 0.105658389;
    double Mp = 0.93827231;
    double Mn = 0.93956563;
    
    double Ml;
    if (lepton_type == kNUE or lepton_type == kNUEBAR)
      Ml = Me;
    else if (lepton_type == kNUMU or lepton_type == kNUMUBAR)
      Ml = Mm;
    else {
      std::cerr << "Unknown lepton type" << std::endl;
      throw std::exception();
    }
    
    double Etot  = (evis_MeV/1000.) + Ml;
    double EnuQE = 0.5*(2.*Mn*Etot + Mp*Mp - Mn*Mn - Ml*Ml) / (Mn - Etot + costheta*std::sqrt(Etot*Etot - Ml*Ml));

    return EnuQE;
  }

  double CombinedFit_Q2_ryan(const double evis_MeV, const double costheta, const NuType_t lepton_type) {
    
    double Me = 0.00051099906;
    double Mm = 0.105658389;
    double Mp = 0.93827231;
    double Mn = 0.93956563;

    double Ml;
    if (lepton_type == kNUE or lepton_type == kNUEBAR)
      Ml = Me;
    else if (lepton_type == kNUMU or lepton_type == kNUMUBAR)
      Ml = Mm;
    else {
      std::cerr << "Unknown lepton type" << std::endl;
      throw std::exception();
    }

    double Etot  = (evis_MeV/1000.) + Ml;
    double EnuQE = 0.5*(2.*Mn*Etot + Mp*Mp - Mn*Mn - Ml*Ml) / (Mn - Etot + costheta*std::sqrt(Etot*Etot - Ml*Ml));
    double Q2 = 2.0*EnuQE*Etot*(1.0-costheta);
    
    return Q2;
  }

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
		      std::vector<float>& vtx) {

    
    // returns Event_num_pi0s
    
    // nfsp = number of final state particles
    // ipfs = final state geant3 ID
    // vrtx = vertex XYZ of each final state particle
    // pfsp = 4 momentum of each final state particle
    // ldalitz = is dalitz decay?
    // ppi = Event pi0 momentum
    // pgram = Event momentum of each gamma?
    // vtx = Pi0 decay point?
    // pg = 4 momentum of gamma?
    
    bool ld[127];
    unsigned nnpi0;
    GEANT3Type_t it, it2;
    
    float vrtx[3][127];
    float pfsp[4][127];
    float pgam[4][2];

    for(unsigned i=0; i<nfsp; ++i){
      vrtx[0][i] = vrtx_x[i];
      vrtx[1][i] = vrtx_y[i];
      vrtx[2][i] = vrtx_z[i];

      pfsp[0][i] = pfsp_x[i];
      pfsp[1][i] = pfsp_y[i];
      pfsp[2][i] = pfsp_z[i];
      pfsp[3][i] = pfsp_t[i];
    }
    
    unsigned ng;
    unsigned mg;
    float pg[4][127];
    float mtemp;

    float mpi0;
    float me;
    mpi0 = 0.1349766;
    me   = 0.511e-3;

    unsigned npi0 = 0;
    ng = 0;
    mg = 0;
    ldalitz = 0;

    // it = GEANT3 particle ID number
    // i loops over final state particles
    for(unsigned i = 0; i < nfsp; ++i) {
      it = (GEANT3Type_t)ipfs[i];
      if (it == kPION0) {
	npi0 = npi0 + 1;
	for(unsigned j = 0; j<3; ++j){
	  vtx[j] = vrtx[j][i];
	}
      }
      
      if (it == kGAMMA) {
	ng = ng + 1;
	ld[ng] = false;
	pg[3][ng] = pfsp[3][i]; // store this gamma energy

	for(unsigned j=0; j<3; ++j) {
	  pg[j][ng] = pg[3][ng] * pfsp[j][i]; // store the gamma vertex
	}
	
	if(pg[3][ng] > 0.025) mg = mg + 1; // if above certain threshold increment mg (?)
      }

      if (it == kPOSITRON or it == kELECTRON) {
	for(unsigned j = i+1; j < nfsp; ++j) {
	  it2 = (GEANT3Type_t)ipfs[j];
	  if ((it2 == kPOSITRON or it2 == kELECTRON) and it != it2) {
	    // one if positron, one is electron, its a gamma?
	    ng = ng + 1;
	    ld[ng] = true;
	    pg[3][ng] = std::sqrt(pfsp[3][i]*pfsp[3][i] + me*me) + std::sqrt(pfsp[3][j]*pfsp[3][j] + me*me);
	    for(unsigned k=0; k<3; ++k) {
	      pg[k][ng] = pfsp[3][i]*pfsp[k][i] + pfsp[3][j]*pfsp[k][j];
	    }
	    if(pg[3][ng] > 0.025) mg = mg + 1;
	  } 
	}
      }
      
    } // end loop over final state particle

    nnpi0 = 0;
    //pairwise loop over gamma?
    for(unsigned i=0; i < ng; ++i) {
      for(unsigned j=i+1; j < ng; ++j) {

	// m = sqrt((E1+E2)^2 - (p)^2)
	
	mtemp = std::sqrt(std::pow(pg[3][i]+pg[3][j],2) - std::pow(pg[0][i]+pg[0][j],2) -
			  std::pow(pg[1][i]+pg[1][j],2) - std::pow(pg[2][i]+pg[2][j],2));
	  
	  if(std::abs(mtemp - mpi0) < 0.000001) {
	    
            nnpi0 = nnpi0 + 1;

	    if(nnpi0 == 1 or (ppi[3] < (pg[3][i] + pg[3][j]))) {
	      for(unsigned k=0; k<4; ++k) {
		ppi[k] = pg[k][i] + pg[k][j];
	      }
		
	      if (ld[i] or  ld[j]) {
		ldalitz = 1;
	      }
		
	      if (pg[3][i] > pg[3][j]) {
		for(unsigned k=0; k<4; ++k) {
		  pgam[k][0] = pg[k][i];
		  pgam[k][1] = pg[k][j];
		}
	      } else {
		for(unsigned k=0; k<4; ++k) {
		  pgam[k][0] = pg[k][j];
		  pgam[k][1] = pg[k][i];
		}
	      }
	    }
	  } // if pizero mass is satisfied
      } // end that gamma
    } // end this gamma

    for(unsigned i=0; i<2; ++i) {
      pgam_x[i] = pgam[0][i];
      pgam_y[i] = pgam[1][i];
      pgam_z[i] = pgam[2][i];
      pgam_t[i] = pgam[3][i];
    }      

    npi0 = std::max(npi0,nnpi0);
    
    return npi0;
  }

  unsigned Pi0Details(const int nfsp,
  		      const std::vector<int>& ipfs,
  		      const std::vector<float>& vrtx_x,
  		      const std::vector<float>& vrtx_y,
  		      const std::vector<float>& vrtx_z,
  		      const std::vector<float>& pfsp_x,
  		      const std::vector<float>& pfsp_y,
  		      const std::vector<float>& pfsp_z,
  		      const std::vector<float>& pfsp_t) {

    static int ldalitz;
    static std::vector<float> ppi(4,0);
    static std::vector<float> pgam_x(2,0);
    static std::vector<float> pgam_y(2,0);
    static std::vector<float> pgam_z(2,0);
    static std::vector<float> pgam_t(2,0);
    static std::vector<float> vtx(3,0);
    
    return Pi0Details(nfsp,
  		      ipfs,
  		      vrtx_x,vrtx_y,vrtx_z,
  		      pfsp_x,pfsp_y,pfsp_z,pfsp_t,
  		      ldalitz,
  		      ppi,
  		      pgam_x,pgam_y,pgam_z,pgam_t,
  		      vtx);
  }

  unsigned Pi0Details(const int nfsp,
		      PyObject* ipfs,
		      PyObject* vrtx_x,
		      PyObject* vrtx_y,
		      PyObject* vrtx_z,
		      PyObject* pfsp_x,
		      PyObject* pfsp_y,
		      PyObject* pfsp_z,
		      PyObject* pfsp_t) {

    return Pi0Details(nfsp,
    		      as_vector_int32(ipfs),
    		      as_vector_float32(vrtx_x),
    		      as_vector_float32(vrtx_y),
    		      as_vector_float32(vrtx_z),
    		      as_vector_float32(pfsp_x),
    		      as_vector_float32(pfsp_y),
    		      as_vector_float32(pfsp_z),
    		      as_vector_float32(pfsp_t));
  }
  

  StackedBkgdType_t StackHistoBkgd(const bool Event_is_dirt, // tuple var 11
				   const bool Event_is_pi0,  // tuple var 10
				   const NuanceType_t evwt,  // tuple var 3
				   const NuType_t nutype,    // tuple var 4
				   const GEANT3Type_t parent_id) // tuple var 8

  {
    StackedBkgdType_t bkg = kBKGD_INVALID;
    
    if      (Event_is_dirt) return kBKGD_DIRT;
    else if (Event_is_pi0)  return kBKGD_PI0;
    else if (evwt == kCC1pNg or evwt == kNC1pNg) return kBKGD_DELTA;
    else if (parent_id == kPIONP and nutype == kNUE) return kBKGD_NUEPIP;
    else if (parent_id == kKAONP and nutype == kNUE) return kBKGD_NUEKP;
    else if (parent_id == kKAON0LONG and nutype == kNUE) return kBKGD_NUEK0;
    else return kBKGD_OTHER;
    
    return bkg;
  }
  


  float StackedHistoBkdgWeight( const bool Event_is_dirt,   // tuple var 11
			        const bool Event_is_pi0,    // tuple var 10
				const NuanceType_t evwt,    // tuple var 3
				const NuType_t nutype,      // tuple var 4
				const GEANT3Type_t parent_id, // tuple var 8
				const float in_weight,
				const float Npi,
				const float Nka,
				const float Nbkg) 
  {
    float out_weight = in_weight;
    
    //out_weight *= FLUX_WEIGHT(INU_COCKT) * Nbkg;

    if (!(Event_is_pi0 or Event_is_dirt or evwt == kNC1pNg)) {

      if (parent_id == kPIONP) out_weight *= Npi;
      if (parent_id == kKAONP) out_weight *= Nka;
	  
    }

    return out_weight;
  }
  

}

#endif
