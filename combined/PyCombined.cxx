#ifndef __PYCOMBINED_CXX__
#define __PYCOMBINED_CXX__


#import "PyCombined.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#include <exception>
#include <iostream>

namespace sp {

  void SetPyUtil() {
    static bool once=false;
    if(!once) { import_array(); once=true; }
  }
  
  std::vector<int> as_vector_int32(PyObject* pyarray) {
    SetPyUtil();

    int *carray;
    const int dtype = NPY_INT32;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp dims[1];
    if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 1, descr) < 0) {
      std::cerr << "cannot convert to 2D C-array" << std::endl;
      throw std::exception();
    }
    
    return std::vector<int>(carray,carray + dims[0]);
  }

  std::vector<float> as_vector_float32(PyObject* pyarray) {
    SetPyUtil();

    float *carray;
    const int dtype = NPY_FLOAT;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp dims[1];
    if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 1, descr) < 0) {
      std::cerr << "cannot convert to 2D C-array" << std::endl;
      throw std::exception();
    }
    
    return std::vector<float>(carray,carray + dims[0]);
  }

  PyObject* as_array_float32(const std::vector<float>& vec) {
    SetPyUtil();
        
    npy_intp dims[1];
    dims[0] = vec.size();
    auto array = (PyArrayObject*)(PyArray_ZEROS(1,dims,NPY_FLOAT,0));
    
    float* data_ptr = (float*) PyArray_DATA(array);

    memcpy(data_ptr, vec.data(), sizeof(float)*(vec.size()));

    return PyArray_Return(array);
  }
  
  
  Eigen::MatrixXf as_mat_float32(PyObject* pyarray) {
    SetPyUtil();

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
    SetPyUtil();
        
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
  
}

#endif
