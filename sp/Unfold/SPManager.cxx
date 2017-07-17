#ifndef __SPMANAGER_CXX_
#define __SPMANAGER_CXX__

#include "SPManager.h"

namespace sp {

  void SPManager::Initialize() {
	  for(const auto &response : _spio.Responses())
	    for(auto algo : _algo_v)
      		algo->Initialize(&response);
  }

  void SPManager::Process() {
      //for(auto algo : _algo_v) 
	//algo->Unfold();
  } 

  void SPManager::Finalize() {}

}

#endif
