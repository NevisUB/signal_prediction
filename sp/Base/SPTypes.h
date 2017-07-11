#ifndef __SPTYPES_H__
#define __SPTYPES_H__

#include <limits>

namespace sp {

  /// Used as an invalid value identifier for long long
  const long long          kINVALID_LONGLONG  = std::numeric_limits< long long          >::max();
  /// Used as an invalid value identifier for unsigned long long
  const unsigned long long kINVALID_ULONGLONG = std::numeric_limits< unsigned long long >::max();
  /// Used as an invalid value identifier for size_t
  const size_t             kINVALID_SIZE      = std::numeric_limits< size_t             >::max();
  /// Used as an invalid value identifier for int
  const int                kINVALID_INT       = std::numeric_limits< int                >::max();
  /// Used as an invalid value identifier for unsigned int
  const unsigned int       kINVALID_UINT      = std::numeric_limits< unsigned int       >::max();
  /// Used as an invalid value identifier for unsigned short
  const short              kINVALID_SHORT     = std::numeric_limits< short              >::max();
  /// Used as an invalid value identifier for unsigned unsigned short
  const unsigned short     kINVALID_USHORT    = std::numeric_limits< unsigned short     >::max();
  /// Used as an invalid value idnetifier for single-point precision  
  const float              kINVALID_FLOAT     = std::numeric_limits< float              >::max();
  /// Used as an invalid value idnetifier for double-point precision
  const double             kINVALID_DOUBLE    = std::numeric_limits< double             >::max();
  
}


#endif 
