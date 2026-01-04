#!/usr/bin/bash
# shellcheck disable=SC1091

export CMAKE_Fortran_COMPILER=ifx
source "${ONEAPI_ROOT}/compiler/latest/env/vars.sh"
source "${ONEAPI_ROOT}/mkl/latest/env/vars.sh"
source "${ONEAPI_ROOT}/mpi/latest/env/vars.sh"

cmake -B "build/debug" \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_Fortran_COMPILER=ifx

cmake --build "build/debug"

chmod +x tshell.sh && cp tshell.sh build