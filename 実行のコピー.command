#!/bin/bash

cd ~/fortran/コード
gfortran -c para.f90
#gfortran -c cortra_jyusin.f90
gfortran -c calcu.f90
gfortran -c record.f90
gfortran -o test arpt.f90 para.o calcu.o record.o #cortra_jyusin.o
./test
