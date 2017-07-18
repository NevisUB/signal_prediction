#ifndef __RESPONSE_H__
#define __RESPONSE_H__

#include "Parameter.h"

#include "TH2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"

namespace sp {

  class Response : public TObject{

  public:

    Response();
    Response(Parameter* true_p, Parameter* reco_p);
    ~Response(){}
    
    std::string _name;
    
    Parameter* _true_param;
    Parameter* _reco_param;

    TMatrixD _response;
    TH2D _response_h;
    
    bool _from_file;

    inline bool operator==(const Response& rhs) const {
      if ( (*_true_param) != *(rhs._true_param) ) return false;
      if ( (*_reco_param) != *(rhs._reco_param) ) return false;
      return true;
    }

    void Fill(float weight, bool passosc);
    void Finalize();
    
    bool filled() const;
    void dump() const;

    ClassDef(sp::Response,1); 
  };
  
}

#endif
