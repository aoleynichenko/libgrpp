#!/bin/bash

test_libgrpp_c.x --grpp > test.out

libgrpp_compare_matrices.x libgrpp_c_arep.txt arep.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x.txt so_x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y.txt so_y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z.txt so_z.txt
if [ $? -ne 0 ]; then exit 1; fi

#libgrpp_compare_matrices.x libgrpp_c_overlap.txt overlap.txt
#if [ $? -ne 0 ]; then exit 1; fi
#libgrpp_compare_matrices.x libgrpp_c_nucattr.txt coulomb.txt
#if [ $? -ne 0 ]; then exit 1; fi


#test_libgrpp_f90.x

#libgrpp_compare_matrices.x libgrpp_f90_arep.txt arep.txt
#if [ $? -ne 0 ]; then exit 1; fi
#libgrpp_compare_matrices.x libgrpp_f90_so_x.txt so_x.txt
#if [ $? -ne 0 ]; then exit 1; fi
#libgrpp_compare_matrices.x libgrpp_f90_so_y.txt so_y.txt
#if [ $? -ne 0 ]; then exit 1; fi
#libgrpp_compare_matrices.x libgrpp_f90_so_z.txt so_z.txt
#if [ $? -ne 0 ]; then exit 1; fi

exit 0


