#ifndef __PYCOMBINED_H__
#define __PYCOMBINED_H__

struct _object;
typedef _object PyObject;

#ifndef __CLING__
#ifndef __CINT__
#include <Python.h>
#endif
#endif
#include <vector>


namespace sp {

  void SetPyUtil();
  std::vector<int> as_vector_int32(PyObject* pyarray);
  std::vector<float> as_vector_float32(PyObject* pyarray);


}

#endif
