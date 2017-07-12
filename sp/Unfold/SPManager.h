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

  private:
    SPIO _spio;
    std::vector<UnfoldAlgoBase*> _algo_v;
    //std::vector<std::tuple<Response*,

  public:

    // Setters
    void AddAlgo(UnfoldAlgoBase* algo);

    // Getters
    SPIO& get_io() { return _spio; } 

  };
}

#endif
