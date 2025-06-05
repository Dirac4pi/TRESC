#!/bin/bash
ulimit -s 262144  # stack size limit to 256 MB
ulimit -c unlimited  # enable core generation
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
    exec "$prog" "$@"
    ;;
  #----------------------------------------------------------------------------
  # running hints
  "-h"|"-H"|"--help"|"--hints")
    "$prog" "-v"
    echo " when you call tshell:" >&2
    echo " user -> |tshell -> tkernel|" >&2
    echo "         |------TRESC------|" >&2
    echo " here are some hints for you" >&2
    echo " tshell" >&2
    echo "        inputfilename.xyz      electronic structure calculation" >&2
    echo "        -v                     program version" >&2
    echo "        -h                     running hints" >&2
    echo "        -1c                    scalar orbital visualization" >&2
    echo "        -2c                    spinor orbital visualization" >&2
    echo "        -d                     debug (need tkerneld)" >&2
    echo "        -td                    vtune threading analysis" >&2
    echo "        -hp                    vtune hotspots analysis" >&2
    echo "        -cd                    core dump (need tkerneld)" >&2
    echo " For more: https://github.com/Dirac4pi/TRESC" >&2
    ;;
  #----------------------------------------------------------------------------
  # program version
  "-v"|"-V"|"--version"|"--VERSION")
    exec "$prog" "-v"
    ;;
  #----------------------------------------------------------------------------
  # orbital visualization
  "-1c"|"-1C"|"-2c"|"-2C")
    export KMP_AFFINITY="granularity=fine,compact,1,0"
    export MKL_NUM_THREADS=8
    export OMP_STACKSIZE=256M
    exec "$prog" "$@"
    ;;
  #----------------------------------------------------------------------------
  # core dump
  "-cd"|"-CD"|"-coredump"|"-COREDUMP")
    export FOR_DUMP_CORE_FILE=TRUE
    echo "core.%p" | sudo tee /proc/sys/kernel/core_pattern
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
    exec "tkerneld" "$inp"
    ;;
  #----------------------------------------------------------------------------
  # vtune hotspots analysis
  "-hotspots"|"-hp"|"-HP"|"-HOTSPOTS")
    source $ONEAPI_ROOT/vtune/latest/vtune-vars.sh
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
    source $ONEAPI_ROOT/vtune/latest/vtune-vars.sh
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
