me="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $me

export COMBINED_BASEDIR=$me
export COMBINED_LIBDIR=$me/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$COMBINED_LIBDIR
