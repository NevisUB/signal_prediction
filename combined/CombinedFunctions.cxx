#ifndef __COMBINEDFUNCTIONS_CXX__
#define __COMBINEDFUNCTIONS_CXX__

#include "CombinedFunctions.h"

#include <iostream>
namespace sp { 

  std::string Bkgd2String(BkgdType_t bkgd) {

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

  BkgdType_t CombinedFit_bkgd_type(NuanceType_t evwt, NuType_t inno) {

    BkgdType_t ibkgd = kBINVALID;
    
    if ((inno == kNUMU) or (inno == kNUE)) {
      if (evwt == kCCQE) ibkgd = kBCCQE;
      else if (evwt == kCC1p1pi  or evwt == kCC1n1pi)  ibkgd = kBCCPIP;
      else if (evwt == kNC1p1pi0 or evwt == kNC1n1pi0) ibkgd = kBCHPI0;
      else if (evwt == kNC1pi01a or evwt == kCC1pi1a)  ibkgd = kBCHPI0;
      else if (evwt == kCC1pNg   or evwt == kCC1pNg)   ibkgd = kBDELTA;
    }
    else  {
      if ( evwt == kCCQE ) ibkgd = kBCCQE;
      else if (evwt == kCCbar1p1pi  or evwt == kCCbar1n1pi)  ibkgd = kBCCPIP;
      else if (evwt == kNCbar1p1pi0 or evwt == kNCbar1n1pi0) ibkgd = kBCHPI0;
      else if (evwt == kNC1pi01a    or evwt == kCC1pi1a)     ibkgd = kBCHPI0;
      else if (evwt == kCC1pNg      or evwt == kCC1pNg)      ibkgd = kBDELTA;
    }
    
    return ibkgd;
  }

}

#endif
