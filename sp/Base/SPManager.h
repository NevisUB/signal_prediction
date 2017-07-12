#ifndef __SPMANAGER_H__
#define __SPMANAGER_H__

#include "SPIO.h"
#include "UnfoldAlgoBase.h"

namespace sp {

  class SPManager {
  public:
    SPManager() {}
    ~SPManager() {}

    void Initialize();
    void Process();
    void Finalize();

    void AddAlgo(UnfoldAlgoBase* algo);
    
  private:
    SPIO _spio;

    std::vector<UnfoldAlgoBase*> _algo_v;
  };


}

#endif
