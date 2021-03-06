#ifndef __PYALGORITHM_CXX__
#define __PYALGORITHM_CXX__

#include "PyAlgorithm.h"
#include "AlgorithmUtil.h"

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

#include <numpy/ndarrayobject.h>

#include <exception>
#include <iostream>

namespace sp {

  void InitPyAlgorithm() {
    static bool once=false;
    if(!once) { import_array(); once=true; }
  }
  
  Eigen::MatrixXf as_mat_float32(PyObject* pyarray) {
    InitPyAlgorithm();
    
    float **carray;
    // Create C arrays from numpy objects:
    const int dtype = NPY_FLOAT;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp dims[2];
    if (PyArray_AsCArray(&pyarray, (void **)&carray, dims, 2, descr) < 0) {
      std::cerr << "cannot convert to 2D C-array" << std::endl;
      throw std::exception();
    }

    auto mat = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>( carray[0], dims[0], dims[1] );

    PyArray_Free(pyarray,(void *)carray);

    return mat;
  }


  PyObject* as_array_float32(const Eigen::MatrixXf& mat) {
    InitPyAlgorithm();
        
    npy_intp dims[2];
    dims[0] = mat.rows();
    dims[1] = mat.cols();
    auto array = (PyArrayObject*)(PyArray_ZEROS(2,dims,NPY_FLOAT,0));
    
    float* data_ptr = (float*) PyArray_DATA(array);

    float* eigen_ptr = (float*) std::malloc(sizeof(float)*mat.rows()*mat.cols());

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >( eigen_ptr, dims[0], dims[1] ) = mat;

    memcpy(data_ptr, eigen_ptr, sizeof(float)*(mat.rows()*mat.cols()));

    std::free(eigen_ptr);
    
    return PyArray_Return(array);
  }

  PyObject* as_array_float32(const Eigen::VectorXf& vec) {
    return as_array_float32(to_vector_std(vec));
  }
  
}

#endif
