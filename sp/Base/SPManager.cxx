#ifndef __SPMANAGER_CXX_
#define __SPMANAGER_CXX__

#include "SPManager.h"

namespace sp {

  void SPManager::Initialize() {
    for(auto algo : _algo_v)
      algo->Initialize();
  }

  void SPManager::Process() {
    for(const auto& response : _spio.Responses())
      for(auto algo : _algo_v) 
	algo->Unfold(response.response);
  } 

  void SPManager::Finalize() {}

}

#endif
