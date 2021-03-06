#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <iostream>
#include <string>
#include <vector>

#include "TObject.h"
#include "TH1D.h"
#include "UnfoldTypes.h"
#include "TVector.h"
#include "SPModelBase.h"

namespace sp {

  class Parameter : public TObject {
  public:

    Parameter();

    Parameter(const std::vector<std::string>& variable_v,
	      const std::vector<double>& bin_lo_v,
	       SPModelBase *inmodel, Operation_t op = kOP_INVALID);
    
    ~Parameter(){}
    
    std::string _name;
    std::vector<std::string> _variable_v;
    std::vector<double> _bin_lo_v;
    SPModelBase * _model;
    Operation_t _operation;

    bool _filled;
    bool _from_file;

    TH1D _hist;
   
    // Not serialized
    std::vector<float> _data_v;
    
    inline bool operator==(const Parameter& rhs) const {
      if (_variable_v != rhs._variable_v)  return false;
      if (_bin_lo_v != rhs._bin_lo_v) return false;
      return true;
    }

    inline bool operator!=(const Parameter& rhs) const {
      if ((*this) == rhs) return false;
      return true;
    }

    float Fill(float weight, std::string );
    void dump();

    ClassDef(sp::Parameter,1); 
  };

}
#endif
