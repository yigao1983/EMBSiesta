#!/bin/bash

export FC="gfortran"
export MPIFC="${HOME}/Programs/openmpi-1.6/build_Jul-16-2012_gfortran/local/bin/mpif90"
export FCFLAGS="-O3 -ffast-math"
export FPPFLAGS="-DGRID_DP"
export SIESTA_PREFIX="${HOME}/Programs/siesta-3.2-hyb-2013-May-17/Obj/"
export BLAS_LIB="-L${HOME}/Programs/lapack-3.3.1/ -lblas"
export LAPACK_LIB="-L${HOME}/Programs/lapack-3.3.1/ -llapack"
export BLACS_LIB="-L${HOME}/Programs/BLACS/LIB/ -lblacsCinit -lblacsF77init -lblacs"
export SCALAPACK_LIB="-L${HOME}/Programs/scalapack-1.8.0/ -lscalapack"

../Src/configure --prefix="$SIESTA_PREFIX" \
                 --with-blas="$BLAS_LIB" \
                 --with-lapack="$LAPACK_LIB" \
                 --with-blacs="$BLACS_LIB" \
                 --with-scalapack="$SCALAPACK_LIB" \
                 --enable-mpi

# end of script
