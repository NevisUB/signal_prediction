#ifndef __PYCOMBINED_H__
#define __PYCOMBINED_H__

#include "Base/PySP.h"
#include "CombinedFunctions.h"

namespace sp {

  unsigned Pi0Details(const int nfsp,
		      PyObject* ipfs,
		      PyObject* vrtx_x,
		      PyObject* vrtx_y,
		      PyObject* vrtx_z,
		      PyObject* pfsp_x,
		      PyObject* pfsp_y,
		      PyObject* pfsp_z,
		      PyObject* pfsp_t);

}

#endif
