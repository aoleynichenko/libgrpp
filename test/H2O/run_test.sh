#!/bin/bash

test_libgrpp_c.x --grpp --overlap --coulomb-point --kinetic --momentum --overlap-grad --grpp-grad > test.out

#
# pseudopotential
#

libgrpp_compare_matrices.x libgrpp_c_arep.txt arep.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x.txt so_x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y.txt so_y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z.txt so_z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# other one-electron operators
#

libgrpp_compare_matrices.x libgrpp_c_overlap.txt overlap.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_nucattr_point.txt coulomb.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_kinetic_energy.txt kinetic_energy.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_momentum_x.txt momentum_x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_momentum_y.txt momentum_y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_momentum_z.txt momentum_z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# overlap gradients
#

libgrpp_compare_matrices.x libgrpp_c_overlap_grad_0x.txt overlap_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_overlap_grad_0y.txt overlap_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_overlap_grad_0z.txt overlap_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_overlap_grad_1x.txt overlap_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_overlap_grad_1y.txt overlap_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_overlap_grad_1z.txt overlap_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_overlap_grad_2x.txt overlap_grad_2x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_overlap_grad_2y.txt overlap_grad_2y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_overlap_grad_2z.txt overlap_grad_2z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# grpp (arep) gradients
#

libgrpp_compare_matrices.x libgrpp_c_arep_grad_0x.txt arep_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_arep_grad_0y.txt arep_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_arep_grad_0z.txt arep_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_arep_grad_1x.txt arep_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_arep_grad_1y.txt arep_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_arep_grad_1z.txt arep_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_arep_grad_2x.txt arep_grad_2x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_arep_grad_2y.txt arep_grad_2y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_arep_grad_2z.txt arep_grad_2z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# grpp (so-x) gradients
#

libgrpp_compare_matrices.x libgrpp_c_so_x_grad_0x.txt so_x_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x_grad_0y.txt so_x_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x_grad_0z.txt so_x_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_so_x_grad_1x.txt so_x_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x_grad_1y.txt so_x_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x_grad_1z.txt so_x_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_so_x_grad_2x.txt so_x_grad_2x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x_grad_2y.txt so_x_grad_2y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_x_grad_2z.txt so_x_grad_2z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# grpp (so-y) gradients
#

libgrpp_compare_matrices.x libgrpp_c_so_y_grad_0x.txt so_y_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y_grad_0y.txt so_y_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y_grad_0z.txt so_y_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_so_y_grad_1x.txt so_y_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y_grad_1y.txt so_y_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y_grad_1z.txt so_y_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_so_y_grad_2x.txt so_y_grad_2x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y_grad_2y.txt so_y_grad_2y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_y_grad_2z.txt so_y_grad_2z.txt
if [ $? -ne 0 ]; then exit 1; fi

#
# grpp (so-z) gradients
#

libgrpp_compare_matrices.x libgrpp_c_so_z_grad_0x.txt so_z_grad_0x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z_grad_0y.txt so_z_grad_0y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z_grad_0z.txt so_z_grad_0z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_so_z_grad_1x.txt so_z_grad_1x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z_grad_1y.txt so_z_grad_1y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z_grad_1z.txt so_z_grad_1z.txt
if [ $? -ne 0 ]; then exit 1; fi

libgrpp_compare_matrices.x libgrpp_c_so_z_grad_2x.txt so_z_grad_2x.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z_grad_2y.txt so_z_grad_2y.txt
if [ $? -ne 0 ]; then exit 1; fi
libgrpp_compare_matrices.x libgrpp_c_so_z_grad_2z.txt so_z_grad_2z.txt
if [ $? -ne 0 ]; then exit 1; fi


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


