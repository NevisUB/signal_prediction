#ifndef __STORAGE_H__
#define __STORAGE_H__

#include <string>
#include <vector>

#include "TObject.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>
    

namespace sp {

  class Response;
  
  class Parameter : public TObject {
  public:

    Parameter();

    Parameter(const std::string& name,
	      const std::vector<double>& bin_lo_v,
	      const std::vector<double>& bin_hi_v);
    
    ~Parameter(){ std::cout << "~P @ "<< this  << std::endl;};
    
    std::string _name;
    std::vector<double> _bin_lo_v;
    std::vector<double> _bin_hi_v;

    bool _filled;
    bool _from_file;

    TH1D _hist;
    
    // Not serialized
    Response* _response; //!
    float _data; //!
    
    inline bool operator==(const Parameter& rhs) const {
      if (_name     != rhs._name)     return false;
      if (_bin_lo_v != rhs._bin_lo_v) return false;
      if (_bin_hi_v != rhs._bin_hi_v) return false;
      return true;
    }

    inline bool operator!=(const Parameter& rhs) const {
      if ((*this) == rhs) return false;
      return true;
    }

    void Fill(float weight);
    void dump();

    ClassDef(sp::Parameter,1); 
  };

  
  class Response : public TObject{

  public:

    Response();
    Response(Parameter* true_p, Parameter* reco_p);
    ~Response(){ std::cout << "~R @ " << this << std::endl;}
    
    std::string name;
    
    Parameter* true_param;
    Parameter* reco_param;

    TMatrixD response;
    TVectorD true_v;
    TVectorD reco_v;
    
    TH2D response_h;
    
    bool from_file;

    inline bool operator==(const Response& rhs) const {
      if ( (*true_param) != *(rhs.true_param) ) return false;
      if ( (*reco_param) != *(rhs.reco_param) ) return false;
      return true;
    }

    void Fill(float weight, bool passosc, int nutype);
    void Finalize();
    
    bool filled() const;
    void dump() const;

    ClassDef(sp::Response,1); 
  };
  
}

#endif
