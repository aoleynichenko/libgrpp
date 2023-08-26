#!/bin/bash

test_libgrpp_c.x --coulomb-point --coulomb-point-num --coulomb-ball --coulomb-gauss --coulomb-fermi > test.out

libgrpp_compare_matrices.x libgrpp_c_nucattr_point.txt nucattr_point.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_nucattr_point_numerical.txt nucattr_point_numerical.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_nucattr_charged_ball.txt nucattr_charged_ball.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_nucattr_gaussian.txt nucattr_gaussian.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_nucattr_fermi.txt nucattr_fermi.txt
if [ $? -ne 0 ]; then exit 1; fi

exit 0


