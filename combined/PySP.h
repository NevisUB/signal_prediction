#ifndef __PYSP_H__
#define __PYSP_H__

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
  PyObject* as_array_float32(const std::vector<float>& vec);

}

#endif
