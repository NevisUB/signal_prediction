#ifndef __COMBINEDUTIL_CXX__
#define __COMBINEDUTIL_CXX__

#include "CombinedUtil.h"

namespace sp {

  //
  // Convenience functions
  //

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
 

  std::string StackedBkgd2PrettyString(const StackedBkgdType_t bkgd) {
  //std::vector<std::string> temp_names = {"Other","#nu_{e} from K^{0}","#nu_{e} from K^{+}","#nu_{e} from K^{0}","#Delta #rightarrow N #gamma","#pi^{0} misid","Dirt"};
    switch(bkgd) {
    case kBKGD_INVALID: return "INVALID";
    case kBKGD_DIRT:    return "Dirt";
    case kBKGD_PI0:     return "#pi^{0} misid";
    case kBKGD_DELTA:   return "#Delta #rightarrow N #gamma";
    case kBKGD_NUEPIP:  return "#nu_{e} from #mu";
    case kBKGD_NUEKP:   return "#nu_{e} from K^{+}";
    case kBKGD_NUEK0:   return "#nu_{e} from K^{0}";
    case kBKGD_OTHER:   return "Other";
    case kBKGD_MAX:     return "MAX";
    default:            return "INVALID";
    }
    
    return "INVALID";
  }


std::string StackedBkgd2String(const StackedBkgdType_t bkgd) {

    
    //std::vector<std::string> temp_names = {"Other","#nu_{e} from K^{0}","#nu_{e} from K^{+}","#nu_{e} from K^{0}","#Delta #rightarrow N #gamma","#pi^{0} misid","Dirt"};
    switch(bkgd) {
    case kBKGD_INVALID: return "INVALID";
    case kBKGD_DIRT:    return "DIRT";
    case kBKGD_PI0:     return "PI0";
    case kBKGD_DELTA:   return "DELTA";
    case kBKGD_NUEPIP:  return "NUEPIP";
    case kBKGD_NUEKP:   return "NUEKP";
    case kBKGD_NUEK0:   return "NUEK0";
    case kBKGD_OTHER:   return "OTHER";
    case kBKGD_MAX:     return "MAX";
    default:            return "INVALID";
    }
    
    return "INVALID";
  }

  StackedBkgdType_t String2StackedBkgd(const std::string& name) {

    if (name == "INVALID") return kBKGD_INVALID;
    if (name == "DIRT")    return kBKGD_DIRT;
    if (name == "PI0")     return kBKGD_PI0;
    if (name == "DELTA")   return kBKGD_DELTA;
    if (name == "NUEPIP")  return kBKGD_NUEPIP;
    if (name == "NUEKP")   return kBKGD_NUEKP;
    if (name == "NUEK0")   return kBKGD_NUEK0;
    if (name == "OTHER")   return kBKGD_OTHER;
    if (name == "MAX")     return kBKGD_MAX;
    
    return kBKGD_INVALID;
  }


  std::string GEANT32String(const GEANT3Type_t particle) {

    switch(particle) {
    case kINVALID :    return "INVALID";
    case kGAMMA :      return "GAMMA";
    case kPOSITRON :   return "POSITRON";
    case kELECTRON :   return "ELECTRON";
    case kNEUTRINO :   return "NEUTRINO";
    case kMUONP :      return "MUONP";
    case kMUONM :      return "MUONM";
    case kPION0 :      return "PION0";
    case kPIONP :      return "PIONP";
    case kPIONM :      return "PIONM";
    case kKAON0LONG :  return "KAON0LONG";
    case kKAONP :      return "KAONP";
    case kKAONM :      return "KAONM";
    case kNEUTRON :    return "NEUTRON";
    case kPROTON :     return "PROTON";
    case kANTIPROTON : return "ANTIPROTON";
    case kKAON0SHORT : return "KAON0SHORT";
    case kETA :        return "ETA";
    case kLAMBDA :     return "LAMBDA";
    case kSIGMAP :     return "SIGMAP";
    case kSIGMA0 :     return "SIGMA0";
    case kSIGMAM :     return "SIGMAM";
    case kXI0 :        return "XI0";
    case kXIM :        return "XIM";
    case kOMEGAM :     return "OMEGAM";
    case kANTINEUTRON :return "ANTINEUTRON";
    case kANTILAMBDA : return "ANTILAMBDA";
    case kANTISIGMAM : return "ANTISIGMAM";
    case kANTISIGMA0 : return "ANTISIGMA0";
    case kANTISIGMAP : return "ANTISIGMAP";
    case kANTIXI0 :    return "ANTIXI0";
    case kANTIXIP :    return "ANTIXIP";
    case kANTIOMEGAP : return "ANTIOMEGAP";
    case kDEUTERON :   return "DEUTERON";
    case kTRITON :     return "TRITON";
    case kALPHA :      return "ALPHA";
    case kGEANTINO :   return "GEANTINO";
    case kHE3 :        return "HE3";
    case kCERENKOV :   return "CERENKOV";
    default: return "INVALID";
    }

    return "INVALID";
  }


  GEANT3Type_t String2GEANT3(const std::string& name) {

    if( name == "INVALID" ) return kINVALID;
    if( name == "GAMMA" ) return kGAMMA;
    if( name == "POSITRON" ) return kPOSITRON;
    if( name == "ELECTRON" ) return kELECTRON;
    if( name == "NEUTRINO" ) return kNEUTRINO;
    if( name == "MUONP" ) return kMUONP;
    if( name == "MUONM" ) return kMUONM;
    if( name == "PION0" ) return kPION0;
    if( name == "PIONP" ) return kPIONP;
    if( name == "PIONM" ) return kPIONM;
    if( name == "KAON0LONG" ) return kKAON0LONG;
    if( name == "KAONP" ) return kKAONP;
    if( name == "KAONM" ) return kKAONM;
    if( name == "NEUTRON" ) return kNEUTRON;
    if( name == "PROTON" ) return kPROTON;
    if( name == "ANTIPROTON" ) return kANTIPROTON;
    if( name == "KAON0SHORT" ) return kKAON0SHORT;
    if( name == "ETA" ) return kETA;
    if( name == "LAMBDA" ) return kLAMBDA;
    if( name == "SIGMAP" ) return kSIGMAP;
    if( name == "SIGMA0" ) return kSIGMA0;
    if( name == "SIGMAM" ) return kSIGMAM;
    if( name == "XI0" ) return kXI0;
    if( name == "XIM" ) return kXIM;
    if( name == "OMEGAM" ) return kOMEGAM;
    if( name == "ANTINEUTRON" ) return kANTINEUTRON;
    if( name == "ANTILAMBDA" ) return kANTILAMBDA;
    if( name == "ANTISIGMAM" ) return kANTISIGMAM;
    if( name == "ANTISIGMA0" ) return kANTISIGMA0;
    if( name == "ANTISIGMAP" ) return kANTISIGMAP;
    if( name == "ANTIXI0" ) return kANTIXI0;
    if( name == "ANTIXIP" ) return kANTIXIP;
    if( name == "ANTIOMEGAP" ) return kANTIOMEGAP;
    if( name == "DEUTERON" ) return kDEUTERON;
    if( name == "TRITON" ) return kTRITON;
    if( name == "ALPHA" ) return kALPHA;
    if( name == "GEANTINO" ) return kGEANTINO;
    if( name == "HE3" ) return kHE3;
    if( name == "CERENKOV" ) return kCERENKOV;

    return kINVALID;
  }
}

#endif
