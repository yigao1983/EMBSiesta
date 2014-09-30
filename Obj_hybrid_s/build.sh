#!/bin/bash

export FC="lf95"
# export MPIFC="${HOME}/Programs/openmpi-1.6/build_Jul-13-2012_lahey/local/bin/mpif90"
export FCFLAGS="-O3"
export FPPFLAGS="-DGRID_DP"
export SIESTA_PREFIX="${HOME}/Programs/siesta-3.2-hyb-2013-May-17/Obj_hybrid_s/"
export BLAS_LIB="-L${HOME}/Programs/lapack-3.4.2/ -lblas"
export LAPACK_LIB="-L${HOME}/Programs/lapack-3.4.2/ -llapack"
# export BLACS_LIB="-L${HOME}/Programs/BLACS-lahey/LIB/ -lblacsCinit -lblacsF77init -lblacs"
# export SCALAPACK_LIB="-L${HOME}/Programs/scalapack-1.8.0-lahey/ -lscalapack"

../Src_hybrid/configure --prefix="$SIESTA_PREFIX" \
                        --with-blas="$BLAS_LIB" \
                        --with-lapack="$LAPACK_LIB"
#                         --with-blacs="$BLACS_LIB" \
#                         --with-scalapack="$SCALAPACK_LIB" \
#                         --enable-mpi

# end of script
