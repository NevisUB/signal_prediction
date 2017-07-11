#ifndef __STORAGE_H__
#define __STORAGE_H__

#include <string>
#include <vector>

#include "TObject.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"

namespace sp {

  class Response;
  
  class Parameter : public TObject {
  public:
    Parameter() : _hist(nullptr) {}
    Parameter(const std::string& name,
	      const std::vector<double>& bin_lo_v,
	      const std::vector<double>& bin_hi_v,
	      std::vector<std::string> options_v = std::vector<std::string>());
    
    ~Parameter(){ if(_hist) delete _hist; _hist=nullptr;}
    
    std::string _name;
    std::vector<double> _bin_lo_v;
    std::vector<double> _bin_hi_v;
    std::vector<double> _bin_w_v;
    std::vector<double> _bin_v;
    std::vector<std::string> options_v;
    
    bool _filled;
    
    TH1D* _hist;

    // Not serialized
    Response* _response; //!
    float _data; //!
    
    inline bool operator==(const Parameter& rhs) const {
      if (_name     != rhs._name)     return false;
      if (_bin_lo_v != rhs._bin_lo_v) return false;
      if (_bin_hi_v != rhs._bin_hi_v) return false;
      if (_bin_v    != rhs._bin_v)    return false;
      if (_filled   != rhs._filled)   return false;
      return true;
    }

    inline bool operator!=(const Parameter& rhs) const {
      if ((*this) == rhs) return false;
      return true;
    }

    void Fill(float weight);
    
    void dump();

    ClassDef(Parameter,1); 
  };

  
  class Response : public TObject{

  public:

    Response() :
      name(""),
      true_param(nullptr),
      reco_param(nullptr),
      response_h(nullptr)
    {}

    Response(Parameter* true_p, Parameter* reco_p);
    
    ~Response(){ if(response_h) delete response_h; response_h=nullptr;}
    
    std::string name;
    
    Parameter* true_param;
    Parameter* reco_param;

    TMatrixD response;
    TVectorD true_v;
    TVectorD reco_v;
    
    TH2D* response_h;
    
    inline bool operator==(const Response& rhs) const {
      if ( (*true_param) != *(rhs.true_param) ) return false;
      if ( (*reco_param) != *(rhs.reco_param) ) return false;
      return true;
    }

    inline bool filled() const { return (true_param->_filled && reco_param->_filled); }
    
    void Fill(float weight, bool passosc, int nutype);
    void Finalize();
    
    void dump();

    ClassDef(Response,1); 
  };
  
}

#endif
