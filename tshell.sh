#!/bin/bash
# shellcheck disable=SC1091
ulimit -s 524288  # stack size limit to 512 MB
export MKL_THREADING_LAYER=INTEL

if [ $# -lt 1 ]; then
  echo "tshell: 0 argument received" >&2
  exit 1
fi

# release or debug
case "$1" in
  -[dD])
    shift
    prog='tkerneld'
    ;;
  *)
    prog='tkernel'
    ;;
esac

# source oneAPI environment
__saved_args=("$@")
source "${ONEAPI_ROOT}/compiler/latest/env/vars.sh"
source "${ONEAPI_ROOT}/mkl/latest/env/vars.sh"
source "${ONEAPI_ROOT}/mpi/latest/env/vars.sh"
set -- "${__saved_args[@]}"
unset __saved_args

case "$1" in
  #----------------------------------------------------------------------------
  # electronic structure calculation
  [!-]*)
    inp=$1
    if [ ! -f "$inp" ]; then
      echo "tshell: File $inp not found!" >&2
      exit 1
    fi
    thr=$(sed -n "{/threads=/ {s/.*threads= *\([0-9]\+\).*/\1/p}}" "$inp")
    if [ -z "$thr" ]; then
      echo "tshell: no threads information, set 8 threads" >&2
      thr=8
    fi
    export KMP_AFFINITY="granularity=fine,compact,1,0"
    export MKL_NUM_THREADS=$thr
    export OMP_STACKSIZE=256M
    "$prog" "$@"
    # send message to nfty
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
      curl --retry 5 --retry-delay 3 -d \
      "TRESC: Normal termination of job $1" ntfy.sh/batchannel
    elif [ $exit_code -eq 1 ]; then
      curl --retry 5 --retry-delay 3 -d \
      "TRESC: Error termination of job $1" ntfy.sh/batchannel
    fi
    ;;
  #----------------------------------------------------------------------------
  # input hints
  "-h"|"-H"|"--help"|"--hints")
    "$prog" "-v"
    echo " when you call tshell:" >&2
    echo " user -> |tshell -> tkernel|" >&2
    echo "         |------TRESC------|" >&2
    echo " here are some hints for you" >&2
    echo " tshell" >&2
    echo "        input.tre      electronic structure calculation" >&2
    echo "        -v             program version" >&2
    echo "        -h             input hints" >&2
    echo "        -mog1c/-mog2c  generate 1c/2c unstructured MO grid" >&2
    echo "        -cub1c/-cub2c  generate 1c/2c structured MO grid" >&2
    echo "                       in real space" >&2
    echo "        -pro1c/-pro2c  generate 1c/2c structured MO grid" >&2
    echo "                       in real space and momentum space" >&2
    echo "        -d             debug (need tkerneld)" >&2
    echo "        -td            vtune threading analysis" >&2
    echo "        -hp            vtune hotspots analysis" >&2
    echo " For more: https://github.com/Dirac4pi/TRESC" >&2
    ;;
  #----------------------------------------------------------------------------
  # program version
  "-v"|"-V"|"--version"|"--VERSION")
    "$prog" "-v"
    ;;
  #----------------------------------------------------------------------------
  # generate structured/unstructured grid in real/momentum space for orbitals
  "-mog1c"|"-MOG1C"|"-mog2c"|"-MOG2C"|"-cub1c"|"-CUB1C"|"-cub2c"|"-CUB2C"|\
  "-pro1c"|"-PRO1C"|"-pro2c"|"-PRO2C")
    export KMP_AFFINITY="granularity=fine,compact,1,0"
    export MKL_NUM_THREADS=16
    export OMP_STACKSIZE=256M
    "$prog" "$@"
    ;;
  #----------------------------------------------------------------------------
  # vtune hotspots analysis
  "-hotspots"|"-hp"|"-HP"|"-HOTSPOTS")
    source "$ONEAPI_ROOT/vtune/latest/vtune-vars.sh"
    echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
    inp=$2
    if [ ! -f "$inp" ]; then
      echo "tshell: File $inp not found!" >&2
      exit 1
    fi
    thr=$(sed -n "{/threads=/ {s/.*threads= *\([0-9]\+\).*/\1/p}}" "$inp")
    if [ -z "$thr" ]; then
      echo "tshell: no threads information, set 8 threads" >&2
      thr=8
    fi
    export KMP_AFFINITY="granularity=fine,compact,1,0"
    export MKL_NUM_THREADS=$thr
    export OMP_STACKSIZE=256M
    vtune -collect hotspots -result-dir vtune_hotspots "$prog" "$inp"
    ;;
  #----------------------------------------------------------------------------
  # vtune threading analysis
  "-threading"|"-THREADING"|"-td"|"-TD"|"-th"|"-TH")
    source "$ONEAPI_ROOT/vtune/latest/vtune-vars.sh"
    echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
    inp=$2
    if [ ! -f "$inp" ]; then
      echo "tshell: File $inp not found!" >&2
      exit 1
    fi
    thr=$(sed -n "{/threads=/ {s/.*threads= *\([0-9]\+\).*/\1/p}}" "$inp")
    if [ -z "$thr" ]; then
      echo "tshell: no threads information, set 8 threads" >&2
      thr=8
    fi
    export KMP_AFFINITY="granularity=fine,compact,1,0"
    export MKL_NUM_THREADS=$thr
    export OMP_STACKSIZE=256M
    vtune -collect threading -result-dir vtune_threading "$prog" "$inp"
    ;;
esac
