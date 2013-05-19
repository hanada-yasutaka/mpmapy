#!/bin/bash
# compile wrapper_qeispack

list=(*.o)
test -e ${list} && rm *.o

fname=wrapper_qeispack
dir=../../shared
test -e ${fname}.so  && rm ${fname}.so 
gver=`gfortran --version | head -1`
echo "${gver}"
test=`echo "${gver: -5:3} > 4.6" | bc`
if [ `echo $test` == 0 ]; then
    echo "Error: require gfortran 4.7 or later"
    exit 1
else
    echo "gfortran -W -fPIC -c ${fname}.f90"
    echo "gfortran -W -fPIC -c qeispack.f90"
    gfortran -W -fPIC -c ${fname}.f90 -lm && \
    gfortran -W -fPIC -c qeispack.f90 -lm && \
    gfortran -shared -o wrapper_qeispack.so wrapper_qeispack.o qeispack.o -lm
    echo "gfortran -shared -o wrapper_qeispack.so wrapper_qeispack.o qeispack.o"    
    if [ `echo $?` -ne 0 ]; then
        echo "Error: compile error"
        exit 1
    fi
    test -e ${dir} || mkdir ${dir}
    cp ${fname}.so ${dir}/
    echo "Success: make ${fname}.so in ${dir/..\//}/ directroy."
    exit 0
fi

