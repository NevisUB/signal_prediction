me="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $me

export SP_TOPDIR=$me
export SP_LIBDIR=$me/lib
export EIGEN_INCDIR=/usr/local/include/eigen3

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SP_LIBDIR
