#ifndef __SPMODELBASE_H__
#define __SPMODELBASE_H__

namespace sp {
  class SPModelBase {
  public:
    SPModelBase() {}
    virtual ~SPModelBase() {}
    virtual bool Valid() = 0;
    
    int NuType;
    int NUANCEChan;

  protected:
    
  };
}

#endif
