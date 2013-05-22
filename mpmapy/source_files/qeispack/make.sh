#!/bin/bash
# compile wrapper_qeispack
CC=gfortran
Flags="-fPIC -c -O1"
libs="-lm"
list=(*.o)
test -e ${list} && rm *.o

fname=wrapper_qeispack
dir=../../shared
test -e ${fname}.so  && rm ${fname}.so 
gver=`$CC --version | head -1`
test=`echo "${gver: -5:3} > 4.6" | bc`
if [ `echo $test` == 0 ]; then
    echo "Warning: require gfortran 4.7 or later"
#    exit 1
else
    echo "$CC $Flags ${fname}.f90 $libs"
    echo "$CC $Flags qeispack.f90 $libs"
    echo "$CC -shared -o wrapper_qeispack.so wrapper_qeispack.o qeispack.o -lm"        
    $CC $Flags ${fname}.f90 $libs && \
    $CC $Flags qeispack.f90 $libs && \
    $CC -shared -o wrapper_qeispack.so wrapper_qeispack.o qeispack.o -lm
    
    if [ `echo $?` -ne 0 ]; then
        echo "Error: compile error"
        exit 1
    fi
    test -e ${dir} || mkdir ${dir}
    cp ${fname}.so ${dir}/
    echo "Success: make ${fname}.so in ${dir/..\//}/ directroy."
    exit 0
fi

