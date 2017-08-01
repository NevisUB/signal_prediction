#ifndef __PYALGORITHM_H__
#define __PYALGORITHM_H__

#include "Base/PySP.h"

#include <Eigen/Dense>

namespace sp {

  void InitPyAlgorithm();
  
  Eigen::MatrixXf as_mat_float32(PyObject* pyarray);
  PyObject* as_array_float32(const Eigen::MatrixXf& mat);
  PyObject* as_array_float32(const Eigen::VectorXf& vec);
}

#endif
