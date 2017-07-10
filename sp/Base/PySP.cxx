#ifndef __PYSP_CXX__
#define __PYSP_CXX__

#import "PySP.h"
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

  std::vector<double> as_vector_float64(PyObject* pyarray) {
    SetPyUtil();

    double *carray;
    const int dtype = NPY_DOUBLE;
    PyArray_Descr *descr = PyArray_DescrFromType(dtype);
    npy_intp dims[1];
    if (PyArray_AsCArray(&pyarray, (void *)&carray, dims, 1, descr) < 0) {
      std::cerr << "cannot convert to 2D C-array" << std::endl;
      throw std::exception();
    }
    
    return std::vector<double>(carray,carray + dims[0]);
  }

  std::vector<double> as_vector_double(PyObject* pyarray) {
    return as_vector_float64(pyarray);
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


  PyObject* as_array_float64(const std::vector<double>& vec) {
    SetPyUtil();
        
    npy_intp dims[1];
    dims[0] = vec.size();
    auto array = (PyArrayObject*)(PyArray_ZEROS(1,dims,NPY_DOUBLE,0));
    
    double* data_ptr = (double*) PyArray_DATA(array);

    memcpy(data_ptr, vec.data(), sizeof(double)*(vec.size()));

    return PyArray_Return(array);
  }


  PyObject* as_array_double(const std::vector<double>& vec) {
    return as_array_float64(vec);
  }

}

#endif
