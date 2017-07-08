#ifndef __PYCOMBINED_CXX__
#define __PYCOMBINED_CXX__

#import "PyCombined.h"

namespace sp {

  unsigned Pi0Details(const int nfsp,
		      PyObject* ipfs,
		      PyObject* vrtx_x,
		      PyObject* vrtx_y,
		      PyObject* vrtx_z,
		      PyObject* pfsp_x,
		      PyObject* pfsp_y,
		      PyObject* pfsp_z,
		      PyObject* pfsp_t) {

    return Pi0Details(nfsp,
    		      as_vector_int32(ipfs),
    		      as_vector_float32(vrtx_x),
    		      as_vector_float32(vrtx_y),
    		      as_vector_float32(vrtx_z),
    		      as_vector_float32(pfsp_x),
    		      as_vector_float32(pfsp_y),
    		      as_vector_float32(pfsp_z),
    		      as_vector_float32(pfsp_t));
  }

}

#endif
