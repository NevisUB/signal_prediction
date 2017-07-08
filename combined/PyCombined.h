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
#include <Eigen/Dense>

namespace sp {

  void SetPyUtil();

  std::vector<int> as_vector_int32(PyObject* pyarray);

  std::vector<float> as_vector_float32(PyObject* pyarray);
  PyObject* as_array_float32(const std::vector<float>& vec);

  Eigen::MatrixXf as_mat_float32(PyObject* pyarray);
  PyObject* as_array_float32(const Eigen::MatrixXf& mat);
}

#endif
