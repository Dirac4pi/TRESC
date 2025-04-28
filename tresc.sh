#!/bin/bash
ulimit -s 65536
# check arguments exists
if [ $# -lt 1 ]; then
  echo "TRESC.sh: 0 argument received"
  exit 1
fi
if [[ "$1" == "-d" || "$1" == "-D" ]]; then
  shift
  prog='thomasdebug'
else
  prog='thomas'
fi

case "$1" in
  [!-]*)
    # electronic structure calculation
    inp=$1
    if [ ! -f "$inp" ]; then
      echo "TRESC.sh: File $inp not found!"
      exit 1
    fi
    # extract 'threads' from input file
    thr=$(grep -A1 "%Hamiltonian" "$inp" \
          | grep "threads=" \
          | awk -F'=' '{print $2}')
    if [ -z "$thr" ]; then
      echo "TRESC.sh: no threads information, set 8 threads"
      thr=8
    fi
    export KMP_AFFINITY="granularity=fine,compact,1,0"
    export MKL_NUM_THREADS=$thr
    export OMP_STACKSIZE=32M
    exec "$prog" "$@"
    ;;
  "-v"|"-V"|"--version")
    # program version
    exec "$prog" "$@"
    ;;
  "-1c"|"-1C"|"-2c"|"-2C")
    # orbital visualization
    cores=$(nproc)
    export MKL_NUM_THREADS=$cores
    exec "$prog" "$@"
    ;;
esac