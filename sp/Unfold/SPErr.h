#ifndef __SPERR_H__
#define __SPERR_H__

#include <iostream>
#include <exception>

namespace sp {

  class sperr : public std::exception {
    
  public:
    sperr(std::string msg="") : std::exception()
    {
      _msg = "\033[93m";
      _msg += msg;
      _msg += "\033[00m";
    }

    virtual ~sperr() throw() {}
    virtual const char* what() const throw()
    { return _msg.c_str(); }

  private:
    std::string _msg;
    
  };
}

#endif
