#!/bin/bash

test_libgrpp_f90.x > test.out

#
# semilocal RPP integrals
#
libgrpp_compare_matrices.x libgrpp_f90_arep.txt arep.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_x.txt so_x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_y.txt so_y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_z.txt so_z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# gradients: averaged part
#
libgrpp_compare_matrices.x libgrpp_f90_arep_grad_0x.txt arep_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_arep_grad_0y.txt arep_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_arep_grad_0z.txt arep_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_arep_grad_1x.txt arep_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_arep_grad_1y.txt arep_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_arep_grad_1z.txt arep_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# gradients: spin-orbit part, x component
#
libgrpp_compare_matrices.x libgrpp_f90_so_x_grad_0x.txt so_x_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_x_grad_0y.txt so_x_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_x_grad_0z.txt so_x_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_x_grad_1x.txt so_x_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_x_grad_1y.txt so_x_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_x_grad_1z.txt so_x_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# gradients: spin-orbit part, y component
#
libgrpp_compare_matrices.x libgrpp_f90_so_y_grad_0x.txt so_y_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_y_grad_0y.txt so_y_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_y_grad_0z.txt so_y_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_y_grad_1x.txt so_y_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_y_grad_1y.txt so_y_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_y_grad_1z.txt so_y_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# gradients: spin-orbit part, z component
#
libgrpp_compare_matrices.x libgrpp_f90_so_z_grad_0x.txt so_z_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_z_grad_0y.txt so_z_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_z_grad_0z.txt so_z_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_z_grad_1x.txt so_z_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_z_grad_1y.txt so_z_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_f90_so_z_grad_1z.txt so_z_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi



















