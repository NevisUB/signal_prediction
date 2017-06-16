#ifndef __COMBINEDTYPES_H__
#define __COMBINEDTYPES_H__

namespace sp {

  // pyROOT doesn't like enum class...
  enum NuType_t
    {
      kNUINVALID=0,
      kNUMU,
      kNUMUBAR,
      kNUE,
      kNUEBAR
    };

  enum BkgdType_t
    {
      kBINVALID=0,
      kBCCQE,
      kBCCPIP,
      kBNCPI0,
      kBCHPI0,
      kBDELTA
    };

  enum NuanceType_t
    {
      kUNKNOWN=0,
      kCCQE=1,
      kNCQE=2,

      // Resonant single pion production
      kCC1p1pi=3,
      kCC1n1pi=5,
      kNC1p1pi0=6,
      kNC1n1pi0=8,
      
      kCCbar1p1pi=10,
      kCCbar1n1pi=12,
      kNCbar1p1pi0=13,
      kNCbar1n1pi0=15,

      // Delta radiative
      kCC1pNg=93,
      kNC1pNg=94,
      kNC1pi01a=96,
      kCC1pi1a=97
    };



}


#endif 
