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


}

#endif
