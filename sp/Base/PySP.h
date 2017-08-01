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
#include "TMatrixD.h"
#include "TVectorD.h"

namespace sp {

  void SetPyUtil();

  std::vector<int>    as_vector_int32  (PyObject* pyarray);
  std::vector<float>  as_vector_float32(PyObject* pyarray);
  std::vector<double> as_vector_float64(PyObject* pyarray);
  std::vector<double> as_vector_double (PyObject* pyarray);
  
  PyObject* as_array_float32(const std::vector<float>& vec);
  PyObject* as_array_float64(const std::vector<double>& vec);
  PyObject* as_array_double (const std::vector<double>& vec);

  PyObject* as_array_float64(const TVectorD& vec);
  PyObject* as_array_double (const TVectorD& vec);

  PyObject* as_array_float64(const TMatrixD& mat);
  PyObject* as_array_double (const TMatrixD& mat);

}

#endif
