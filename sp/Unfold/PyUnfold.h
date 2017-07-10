#ifndef __PYUNFOLD_H__
#define __PYUNFOLD_H__

#include "Base/PySP.h"

#include <Eigen/Dense>

namespace sp {

  void InitPyUnfold();
  
  Eigen::MatrixXf as_mat_float32(PyObject* pyarray);
  PyObject* as_array_float32(const Eigen::MatrixXf& mat);
  
}

#endif
