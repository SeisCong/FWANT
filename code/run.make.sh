#!/bin/bash

if [ $# -eq 1 ]; then
   target=$1
else
   target='all'
fi

hr="-"; while [ ${#hr} -lt 70 ]; do hr=${hr}"-"; done

function compile_code {
make ${target} 2>&1 | tee make.${target}.log
ierr=${PIPESTATUS[0]}
}

time compile_code

echo -e "\n"
echo $hr
if [ $ierr -eq 0 ]; then
   echo "\"make ${target}\" succeeded"
else
   echo "\"make ${target}\" failed, please check make.${target}.log"
fi
echo $hr
