me="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $me

export COMBINED_BASEDIR=$me
export COMBINED_LIBDIR=$me/lib
export EIGEN_INCDIR=/usr/local/include/eigen3

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$COMBINED_LIBDIR
