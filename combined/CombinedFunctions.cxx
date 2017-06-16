#ifndef __COMBINEDFUNCTIONS_CXX__
#define __COMBINEDFUNCTIONS_CXX__

#include "CombinedFunctions.h"

#include <iostream>
#include <cmath>
namespace sp { 

  std::string Bkgd2String(const BkgdType_t bkgd) {

    switch(bkgd) {
    case kBINVALID: return "BINVALID";
    case kBCCQE:    return "BCCQE";
    case kBCCPIP:   return "BCCPIP";
    case kBNCPI0:   return "BNCPI0";
    case kBCHPI0:   return "BCHPI0";
    case kBDELTA:   return "BDELTA";
    default:        return "BINVALID";
    }
    
    return "BINVALID";
  }

  BkgdType_t String2Bkgd(const std::string& name) {

    if (name == "BINVALID") return kBINVALID;
    if (name == "BCCQE")    return kBCCQE;
    if (name == "BCCPIP")   return kBCCPIP;
    if (name == "BNCPI0")   return kBNCPI0;
    if (name == "BCHPI0")   return kBCHPI0;
    if (name == "BDELTA")   return kBDELTA;
    
    return kBINVALID;
  }

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
      std::cout << "Unknown lepton type" << std::endl;
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
      std::cout << "Unknown lepton type" << std::endl;
      throw std::exception();
    }

    double Etot  = (evis_MeV/1000.) + Ml;
    double EnuQE = 0.5*(2.*Mn*Etot + Mp*Mp - Mn*Mn - Ml*Ml) / (Mn - Etot + costheta*std::sqrt(Etot*Etot - Ml*Ml));
    double Q2 = 2.0*EnuQE*Etot*(1.0-costheta);
    
    return Q2;
  }
}

#endif
