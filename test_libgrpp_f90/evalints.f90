!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2023 Alexander Oleynichenko
!

module evalints

contains

    subroutine calculate_ecp_integrals(arep_matrix, so_matrices)

        use ecp
        use basis
        use xyz
        use libgrpp
        implicit none

        ! output matrices
        real(8) :: arep_matrix(:, :)
        real(8) :: so_matrices(:, :, :)

        ! local variables
        integer :: basis_dim
        integer :: iatom1, iatom2
        integer :: iblock1, iblock2
        integer :: ifun1, ifun2
        integer :: nprim
        integer :: ncontr
        integer :: z1, z2
        integer :: L_A, L_B
        integer :: r, s, t
        integer :: cart_A
        integer :: ncart1, ncart2
        integer :: icart1, icart2
        integer :: ishell, jshell
        integer(4) :: nprim_A, nprim_B
        integer :: cart_list_1(100, 3)
        integer :: cart_list_2(100, 3)
        real(8) :: origin_A(3), origin_B(3)
        integer, parameter :: MAX_BUF = 10000
        real(8) :: buf_a(MAX_BUF)
        real(8) :: buf_x(MAX_BUF), buf_y(MAX_BUF), buf_z(MAX_BUF)
        integer :: ic, zc
        integer :: ioffs, joffs
        integer :: i, j
        integer :: iarep
        integer :: L

        integer(4) :: num_oc_shells
        integer(4) :: oc_shells_L(ECP_MAX_OC_SHELLS)
        integer(4) :: oc_shells_J(ECP_MAX_OC_SHELLS)
        integer(4) :: ocpot_num_prim(ECP_MAX_OC_SHELLS)
        integer(4) :: ocpot_powers(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
        real(8) :: ocpot_coeffs(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
        real(8) :: ocpot_alpha(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

        ! outercore shells
        integer(4) :: oc_shells_num_primitives(ECP_MAX_OC_SHELLS)
        real(8) :: oc_shells_coeffs(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
        real(8) :: oc_shells_alpha(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

        print *
        print *, 'evaluation of grpp integrals'
        print *

        arep_matrix = 0.0d0
        so_matrices = 0.0d0

        ioffs = 0
        ishell = 0
        do iatom1 = 1, natoms
            z1 = charges(iatom1)
            do iblock1 = 1, num_blocks(z1)
                do ifun1 = 1, block_num_contr(z1, iblock1)
                    call generate_cartesians(iblock1 - 1, cart_list_1, ncart1)
                    ishell = ishell + 1

                    joffs = 0
                    jshell = 0
                    do iatom2 = 1, natoms
                        z2 = charges(iatom2)
                        do iblock2 = 1, num_blocks(z2)
                            do ifun2 = 1, block_num_contr(z2, iblock2)
                                call generate_cartesians(iblock2 - 1, cart_list_2, ncart2)
                                jshell = jshell + 1

                                L_A = iblock1 - 1
                                L_B = iblock2 - 1
                                nprim_A = block_num_prim(z1, iblock1)
                                nprim_B = block_num_prim(z2, iblock2)

                                print '(a,i3,a,i3)', ' grpp: ishell=', ishell, '    jshell=', jshell

                                ! sum over atoms with pseudopotentials
                                do ic = 1, natoms
                                    zc = charges(ic)
                                    if (n_arep(zc) == 0) then
                                        cycle
                                    end if

                                    ! radially-local part (type-1 integrals)

                                    call libgrpp_type1_integrals( &
                                            coord(iatom1, :), L_A, nprim_A, &
                                            coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                            coord(iatom2, :), L_B, nprim_B, &
                                            coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                            coord(ic, :), ecp_num_prim(zc, n_arep(zc)), ecp_powers(zc, n_arep(zc), :), &
                                            arep(zc, n_arep(zc), :), ecp_alpha(zc, n_arep(zc), :), buf_a &
                                            )
                                    call update_matrix_part(arep_matrix, ioffs, joffs, buf_a, ncart1, ncart2)

                                    ! semilocal averaged part (type-2 integrals)

                                    do iarep = 1, n_arep(zc) - 1
                                        L = iarep - 1
                                        call libgrpp_type2_integrals( &
                                                coord(iatom1, :), L_A, nprim_A, &
                                                coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                                coord(iatom2, :), L_B, nprim_B, &
                                                coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                                coord(ic, :), L, ecp_num_prim(zc, iarep), ecp_powers(zc, iarep, :), &
                                                arep(zc, iarep, :), ecp_alpha(zc, iarep, :), buf_a &
                                                )
                                        call update_matrix_part(arep_matrix, ioffs, joffs, buf_a, ncart1, ncart2)
                                    end do

                                    ! semilocal spin-orbit part

                                    do iarep = 2, n_arep(zc)
                                        L = iarep - 1
                                        call libgrpp_spin_orbit_integrals( &
                                                coord(iatom1, :), L_A, nprim_A, &
                                                coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                                coord(iatom2, :), L_B, nprim_B, &
                                                coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                                coord(ic, :), L, ecp_num_prim(zc, iarep), ecp_powers(zc, iarep, :), &
                                                esop(zc, iarep, :), ecp_alpha(zc, iarep, :), buf_x, buf_y, buf_z &
                                                )
                                        buf_x = buf_x * 2.0 / (2.0 * L + 1)
                                        buf_y = buf_y * 2.0 / (2.0 * L + 1)
                                        buf_z = buf_z * 2.0 / (2.0 * L + 1)
                                        call update_matrix_part(so_matrices(1,:,:), ioffs, joffs, buf_x, ncart1, ncart2)
                                        call update_matrix_part(so_matrices(2,:,:), ioffs, joffs, buf_y, ncart1, ncart2)
                                        call update_matrix_part(so_matrices(3,:,:), ioffs, joffs, buf_z, ncart1, ncart2)
                                    end do

                                    ! non-local part

                                    !        integer(4) :: num_oc_shells
                                    !        integer(4) :: oc_shells_L(ECP_MAX_OC_SHELLS)
                                    !        integer(4) :: oc_shells_J(ECP_MAX_OC_SHELLS)
                                    !        integer(4) :: ecp_num_primitives(ECP_MAX_OC_SHELLS)
                                    !        integer(4) :: ecp_powers(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
                                    !        real(8) :: ecp_coeffs(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
                                    !        real(8) :: ecp_alpha(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

                                    num_oc_shells = 0
                                    do L = 1, ECP_MAXL
                                        do i = 1, n_oc_shells(zc, L)
                                            num_oc_shells = num_oc_shells + 1

                                            oc_shells_L(num_oc_shells) = L - 1
                                            if (L == 1) then
                                                oc_shells_J(num_oc_shells) = 1
                                            elseif (mod(i, 2) == 1) then
                                                oc_shells_J(num_oc_shells) = 2 * (L - 1) - 1
                                            else
                                                oc_shells_J(num_oc_shells) = 2 * (L - 1) + 1
                                            end if

                                            ocpot_num_prim(num_oc_shells) = ecp_num_prim(zc, L)
                                            ocpot_powers(num_oc_shells, :) = ecp_powers(zc, L, :)
                                            ocpot_coeffs(num_oc_shells, :) = ocpot(zc, L, i, :)
                                            ocpot_alpha(num_oc_shells, :) = ecp_alpha(zc, L, :)

                                            oc_shells_num_primitives(num_oc_shells) = ocbas_num_prim(zc, L)
                                            oc_shells_coeffs(num_oc_shells, :) = ocbas_coeffs(zc, L, :, i)
                                            oc_shells_alpha(num_oc_shells, :) = ocbas_alpha(zc, L, :)

                                        end do
                                    end do

                                    call libgrpp_outercore_potential_integrals( &
                                        coord(iatom1, :), L_A, nprim_A, &
                                        coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                        coord(iatom2, :), L_B, nprim_B, &
                                        coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                        coord(ic, :), num_oc_shells, oc_shells_L, oc_shells_J, &
                                        ocpot_num_prim, ocpot_powers, ocpot_coeffs, ocpot_alpha, &
                                        oc_shells_num_primitives, oc_shells_coeffs, oc_shells_alpha, &
                                        buf_a, buf_x, buf_y, buf_z &
                                    )

                                    call update_matrix_part(arep_matrix, ioffs, joffs, buf_a, ncart1, ncart2)
                                    call update_matrix_part(so_matrices(1,:,:), ioffs, joffs, buf_x, ncart1, ncart2)
                                    call update_matrix_part(so_matrices(2,:,:), ioffs, joffs, buf_y, ncart1, ncart2)
                                    call update_matrix_part(so_matrices(3,:,:), ioffs, joffs, buf_z, ncart1, ncart2)

                                end do

                                joffs = joffs + ncart2

                            end do
                        end do
                    end do

                    ioffs = ioffs + ncart1

                end do
            end do
        end do

    end subroutine calculate_ecp_integrals


    subroutine calculate_rpp_integrals_gradients(arep_grad, so_grad)

        use ecp
        use basis
        use xyz
        use libgrpp
        implicit none

        ! output matrices
        real(8) :: arep_grad(:, :, :, :)
        real(8) :: so_grad(:, :, :, :, :)

        ! local variables
        integer :: basis_dim
        integer :: ipoint
        integer :: iatom1, iatom2
        integer :: iblock1, iblock2
        integer :: ifun1, ifun2
        integer :: nprim
        integer :: ncontr
        integer :: z1, z2
        integer :: L_A, L_B
        integer :: r, s, t
        integer :: cart_A
        integer :: ncart1, ncart2
        integer :: icart1, icart2
        integer :: ishell, jshell
        integer(4) :: nprim_A, nprim_B
        integer :: cart_list_1(100, 3)
        integer :: cart_list_2(100, 3)
        real(8) :: origin_A(3), origin_B(3)
        integer, parameter :: MAX_BUF = 10000
        real(8) :: buf_arep_dx(MAX_BUF), buf_arep_dy(MAX_BUF), buf_arep_dz(MAX_BUF)
        real(8) :: buf_so(3,3,MAX_BUF) ! 1st index -- nuc coord, 2nd index -- SO_{x,y,z}
        integer :: ic, zc
        integer :: ioffs, joffs
        integer :: i, j
        integer :: iarep
        integer :: L

        integer(4) :: num_oc_shells
        integer(4) :: oc_shells_L(ECP_MAX_OC_SHELLS)
        integer(4) :: oc_shells_J(ECP_MAX_OC_SHELLS)
        integer(4) :: ocpot_num_prim(ECP_MAX_OC_SHELLS)
        integer(4) :: ocpot_powers(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
        real(8) :: ocpot_coeffs(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
        real(8) :: ocpot_alpha(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

        ! outercore shells
        integer(4) :: oc_shells_num_primitives(ECP_MAX_OC_SHELLS)
        real(8) :: oc_shells_coeffs(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
        real(8) :: oc_shells_alpha(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

        print *
        print *, 'evaluation of gradients of grpp integrals'
        print *

        arep_grad = 0.0d0
        so_grad = 0.0d0

        ioffs = 0
        ishell = 0
        do iatom1 = 1, natoms
            z1 = charges(iatom1)
            do iblock1 = 1, num_blocks(z1)
                do ifun1 = 1, block_num_contr(z1, iblock1)
                    call generate_cartesians(iblock1 - 1, cart_list_1, ncart1)
                    ishell = ishell + 1

                    joffs = 0
                    jshell = 0
                    do iatom2 = 1, natoms
                        z2 = charges(iatom2)
                        do iblock2 = 1, num_blocks(z2)
                            do ifun2 = 1, block_num_contr(z2, iblock2)
                                call generate_cartesians(iblock2 - 1, cart_list_2, ncart2)
                                jshell = jshell + 1

                                L_A = iblock1 - 1
                                L_B = iblock2 - 1
                                nprim_A = block_num_prim(z1, iblock1)
                                nprim_B = block_num_prim(z2, iblock2)

                                print '(a,i3,a,i3)', ' grpp grad: ishell=', ishell, '    jshell=', jshell

                                ! differentiate over 3d points
                                do ipoint = 1, natoms

                                ! sum over atoms with pseudopotentials
                                do ic = 1, natoms
                                    zc = charges(ic)
                                    if (n_arep(zc) == 0) then
                                        cycle
                                    end if

                                    ! radially-local part (type-1 integrals)

                                    call libgrpp_type1_integrals_gradient( &
                                            coord(iatom1, :), L_A, nprim_A, &
                                            coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                            coord(iatom2, :), L_B, nprim_B, &
                                            coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                            coord(ic, :), ecp_num_prim(zc, n_arep(zc)), ecp_powers(zc, n_arep(zc), :), &
                                            arep(zc, n_arep(zc), :), ecp_alpha(zc, n_arep(zc), :), coord(ipoint, :), &
                                            buf_arep_dx, buf_arep_dy, buf_arep_dz &
                                    )
                                    call update_matrix_part(arep_grad(ipoint,1,:,:), ioffs, joffs, buf_arep_dx, ncart1, ncart2)
                                    call update_matrix_part(arep_grad(ipoint,2,:,:), ioffs, joffs, buf_arep_dy, ncart1, ncart2)
                                    call update_matrix_part(arep_grad(ipoint,3,:,:), ioffs, joffs, buf_arep_dz, ncart1, ncart2)

                                    ! semilocal averaged part (type-2 integrals)

                                    do iarep = 1, n_arep(zc) - 1
                                        L = iarep - 1
                                        call libgrpp_type2_integrals_gradient( &
                                                coord(iatom1, :), L_A, nprim_A, &
                                                coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                                coord(iatom2, :), L_B, nprim_B, &
                                                coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                                coord(ic, :), L, ecp_num_prim(zc, iarep), ecp_powers(zc, iarep, :), &
                                                arep(zc, iarep, :), ecp_alpha(zc, iarep, :), coord(ipoint, :), &
                                                buf_arep_dx, buf_arep_dy, buf_arep_dz &
                                        )
                                        call update_matrix_part(arep_grad(ipoint,1,:,:), ioffs, joffs, buf_arep_dx, ncart1, ncart2)
                                        call update_matrix_part(arep_grad(ipoint,2,:,:), ioffs, joffs, buf_arep_dy, ncart1, ncart2)
                                        call update_matrix_part(arep_grad(ipoint,3,:,:), ioffs, joffs, buf_arep_dz, ncart1, ncart2)
                                    end do

                                    ! semilocal spin-orbit part

                                    do iarep = 2, n_arep(zc)
                                        L = iarep - 1

                                        call libgrpp_spin_orbit_integrals_gradient( &
                                                coord(iatom1, :), L_A, nprim_A, &
                                                coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                                coord(iatom2, :), L_B, nprim_B, &
                                                coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                                coord(ic, :), L, ecp_num_prim(zc, iarep), ecp_powers(zc, iarep, :), &
                                                esop(zc, iarep, :), ecp_alpha(zc, iarep, :), coord(ipoint, :), &
                                                buf_so(1,1,:), buf_so(2,1,:), buf_so(3,1,:), &  ! grad_sox_x, grad_sox_y, grad_sox_z
                                                buf_so(1,2,:), buf_so(2,2,:), buf_so(3,2,:), &  ! grad_soy_x, grad_soy_y, grad_soy_z
                                                buf_so(1,3,:), buf_so(2,3,:), buf_so(3,3,:)  &  ! grad_soz_x, grad_soz_y, grad_soz_z
                                        )

                                        ! update output arrays
                                        buf_so = buf_so * 2.0 / (2.0 * L + 1)
                                        do i = 1, 3     ! loop over coordinates wrt differentiation point
                                            do j = 1, 3 ! loop over {SO_x, SO_y, SO_z} components
                                                call update_matrix_part(so_grad(ipoint,i,j,:,:), ioffs, joffs, &
                                                    buf_so(i,j,:), ncart1, ncart2)
                                            end do
                                        end do
                                    end do

                                    ! non-local part

                                    !        integer(4) :: num_oc_shells
                                    !        integer(4) :: oc_shells_L(ECP_MAX_OC_SHELLS)
                                    !        integer(4) :: oc_shells_J(ECP_MAX_OC_SHELLS)
                                    !        integer(4) :: ecp_num_primitives(ECP_MAX_OC_SHELLS)
                                    !        integer(4) :: ecp_powers(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
                                    !        real(8) :: ecp_coeffs(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)
                                    !        real(8) :: ecp_alpha(ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

                                    !num_oc_shells = 0
                                    !do L = 1, ECP_MAXL
                                    !    do i = 1, n_oc_shells(zc, L)
                                    !        num_oc_shells = num_oc_shells + 1

                                    !        oc_shells_L(num_oc_shells) = L - 1
                                    !        if (L == 1) then
                                    !            oc_shells_J(num_oc_shells) = 1
                                    !        elseif (mod(i, 2) == 1) then
                                    !            oc_shells_J(num_oc_shells) = 2 * (L - 1) - 1
                                    !        else
                                    !            oc_shells_J(num_oc_shells) = 2 * (L - 1) + 1
                                    !        end if

                                    !        ocpot_num_prim(num_oc_shells) = ecp_num_prim(zc, L)
                                    !        ocpot_powers(num_oc_shells, :) = ecp_powers(zc, L, :)
                                    !        ocpot_coeffs(num_oc_shells, :) = ocpot(zc, L, i, :)
                                    !        ocpot_alpha(num_oc_shells, :) = ecp_alpha(zc, L, :)

                                    !        oc_shells_num_primitives(num_oc_shells) = ocbas_num_prim(zc, L)
                                    !        oc_shells_coeffs(num_oc_shells, :) = ocbas_coeffs(zc, L, :, i)
                                    !        oc_shells_alpha(num_oc_shells, :) = ocbas_alpha(zc, L, :)

                                    !    end do
                                    !end do

                                    !call libgrpp_outercore_potential_integrals( &
                                    !        coord(iatom1, :), L_A, nprim_A, &
                                    !        coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                    !        coord(iatom2, :), L_B, nprim_B, &
                                    !        coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                    !        coord(ic, :), num_oc_shells, oc_shells_L, oc_shells_J, &
                                    !        ocpot_num_prim, ocpot_powers, ocpot_coeffs, ocpot_alpha, &
                                    !        oc_shells_num_primitives, oc_shells_coeffs, oc_shells_alpha, &
                                    !        buf_a, buf_x, buf_y, buf_z &
                                    !        )

                                    !call update_matrix_part(arep_matrix, ioffs, joffs, buf_a, ncart1, ncart2)
                                    !call update_matrix_part(so_matrices(1,:,:), ioffs, joffs, buf_x, ncart1, ncart2)
                                    !call update_matrix_part(so_matrices(2,:,:), ioffs, joffs, buf_y, ncart1, ncart2)
                                    !call update_matrix_part(so_matrices(3,:,:), ioffs, joffs, buf_z, ncart1, ncart2)

                                end do

                                end do
                                ! end of differentiation wrt 3d coordinates of atoms

                                joffs = joffs + ncart2

                            end do
                        end do
                    end do

                    ioffs = ioffs + ncart1

                end do
            end do
        end do

    end subroutine calculate_rpp_integrals_gradients


    subroutine calculate_overlap_integrals(overlap_matrix)

        use ecp
        use basis
        use xyz
        use libgrpp
        implicit none

        ! output matrices
        real(8) :: overlap_matrix(:, :)

        ! local variables
        integer :: basis_dim
        integer :: iatom1, iatom2
        integer :: iblock1, iblock2
        integer :: ifun1, ifun2
        integer :: nprim
        integer :: ncontr
        integer :: z1, z2
        integer :: L_A, L_B
        integer :: r, s, t
        integer :: cart_A
        integer :: ncart1, ncart2
        integer :: icart1, icart2
        integer :: ishell, jshell
        integer(4) :: nprim_A, nprim_B
        integer :: cart_list_1(100, 3)
        integer :: cart_list_2(100, 3)
        real(8) :: origin_A(3), origin_B(3)
        integer, parameter :: MAX_BUF = 10000
        real(8) :: buf(MAX_BUF)
        integer :: ic, zc
        integer :: mass_number
        integer :: ioffs, joffs
        integer :: i, j
        integer(4) :: err_code
        real(8) :: r_rms
        real(8) :: fermi_c
        real(8) :: fermi_a
        real(8) :: nuclear_model_params(10)
        real(8), parameter :: FERMI_UNITS_TO_ATOMIC = 1.0 / 52917.7249

        print *
        print *, 'evaluation of overlap integrals'
        print *

        overlap_matrix = 0.0d0

        ioffs = 0
        ishell = 0
        do iatom1 = 1, natoms
            z1 = charges(iatom1)
            do iblock1 = 1, num_blocks(z1)
                do ifun1 = 1, block_num_contr(z1, iblock1)
                    call generate_cartesians(iblock1 - 1, cart_list_1, ncart1)
                    ishell = ishell + 1

                    joffs = 0
                    jshell = 0
                    do iatom2 = 1, natoms
                        z2 = charges(iatom2)
                        do iblock2 = 1, num_blocks(z2)
                            do ifun2 = 1, block_num_contr(z2, iblock2)
                                call generate_cartesians(iblock2 - 1, cart_list_2, ncart2)
                                jshell = jshell + 1

                                L_A = iblock1 - 1
                                L_B = iblock2 - 1
                                nprim_A = block_num_prim(z1, iblock1)
                                nprim_B = block_num_prim(z2, iblock2)

                                print '(a,i3,a,i3)', ' overlap: ishell=', ishell, '    jshell=', jshell

                                ! calculate integrals for the shell pair
                                call libgrpp_overlap_integrals( &
                                        coord(iatom1, :), L_A, nprim_A, &
                                        coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                        coord(iatom2, :), L_B, nprim_B, &
                                        coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                        buf &
                                        )

                                call update_matrix_part(overlap_matrix, ioffs, joffs, buf, ncart1, ncart2)

                                joffs = joffs + ncart2

                            end do
                        end do
                    end do

                    ioffs = ioffs + ncart1

                end do
            end do
        end do

    end subroutine calculate_overlap_integrals


    subroutine calculate_nuclear_attraction_integrals(nuc_attr_matrix, nuclear_model)

        use ecp
        use basis
        use xyz
        use libgrpp
        implicit none

        ! output matrices
        real(8) :: nuc_attr_matrix(:, :)
        integer :: nuclear_model

        ! local variables
        integer :: basis_dim
        integer :: iatom1, iatom2
        integer :: iblock1, iblock2
        integer :: ifun1, ifun2
        integer :: nprim
        integer :: ncontr
        integer :: z1, z2
        integer :: L_A, L_B
        integer :: r, s, t
        integer :: cart_A
        integer :: ncart1, ncart2
        integer :: icart1, icart2
        integer :: ishell, jshell
        integer(4) :: nprim_A, nprim_B
        integer :: cart_list_1(100, 3)
        integer :: cart_list_2(100, 3)
        real(8) :: origin_A(3), origin_B(3)
        integer, parameter :: MAX_BUF = 10000
        real(8) :: buf(MAX_BUF)
        integer :: ic, zc
        integer :: mass_number
        integer :: ioffs, joffs
        integer :: i, j
        integer(4) :: err_code
        real(8) :: r_rms
        real(8) :: fermi_c
        real(8) :: fermi_a
        real(8) :: nuclear_model_params(10)
        real(8), parameter :: FERMI_UNITS_TO_ATOMIC = 1.0 / 52917.7249

        print *
        print *, 'evaluation of nuclear attraction integrals'
        print *

        nuc_attr_matrix = 0.0d0

        ioffs = 0
        ishell = 0
        do iatom1 = 1, natoms
            z1 = charges(iatom1)
            do iblock1 = 1, num_blocks(z1)
                do ifun1 = 1, block_num_contr(z1, iblock1)
                    call generate_cartesians(iblock1 - 1, cart_list_1, ncart1)
                    ishell = ishell + 1

                    joffs = 0
                    jshell = 0
                    do iatom2 = 1, natoms
                        z2 = charges(iatom2)
                        do iblock2 = 1, num_blocks(z2)
                            do ifun2 = 1, block_num_contr(z2, iblock2)
                                call generate_cartesians(iblock2 - 1, cart_list_2, ncart2)
                                jshell = jshell + 1

                                L_A = iblock1 - 1
                                L_B = iblock2 - 1
                                nprim_A = block_num_prim(z1, iblock1)
                                nprim_B = block_num_prim(z2, iblock2)

                                print '(a,i3,a,i3)', ' nuc attr: ishell=', ishell, '    jshell=', jshell

                                ! loop over atoms
                                do ic = 1, natoms
                                    zc = charges(ic)

                                    ! define nuclear model
                                    mass_number = mass_numbers(zc + 1)
                                    call libgrpp_estimate_nuclear_rms_radius_johnson_1985(mass_number, r_rms)

                                    if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE) then
                                        ! do nothing
                                    elseif (nuclear_model == LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL) then
                                        nuclear_model_params(1) = r_rms * FERMI_UNITS_TO_ATOMIC
                                    elseif (nuclear_model == LIBGRPP_NUCLEAR_MODEL_GAUSSIAN) then
                                        nuclear_model_params(1) = r_rms * FERMI_UNITS_TO_ATOMIC
                                    else
                                        call libgrpp_estimate_fermi_model_parameters(R_rms, fermi_c, fermi_a, err_code)
                                        nuclear_model_params(1) = fermi_c * FERMI_UNITS_TO_ATOMIC
                                        nuclear_model_params(2) = fermi_a * FERMI_UNITS_TO_ATOMIC
                                    end if

                                    ! calculate integrals for the shell pair
                                    call libgrpp_nuclear_attraction_integrals( &
                                            coord(iatom1, :), L_A, nprim_A, &
                                            coeffs(z1, iblock1, ifun1, :), exponents(z1, iblock1, :), &
                                            coord(iatom2, :), L_B, nprim_B, &
                                            coeffs(z2, iblock2, ifun2, :), exponents(z2, iblock2, :), &
                                            coord(ic, :), zc, nuclear_model, nuclear_model_params, buf &
                                            )

                                    call update_matrix_part(nuc_attr_matrix, ioffs, joffs, buf, ncart1, ncart2)
                                end do

                                joffs = joffs + ncart2

                            end do
                        end do
                    end do

                    ioffs = ioffs + ncart1

                end do
            end do
        end do

    end subroutine calculate_nuclear_attraction_integrals


    subroutine update_matrix_part(A, offset_1, offset_2, B, size_1, size_2)

        implicit none
        real(8), intent(inout) :: A(:,:)
        integer, intent(in) :: offset_1, offset_2
        real(8), intent(in) :: B(:)
        integer, intent(in) :: size_1, size_2

        integer :: i, j

        do i = 1, size_1
            do j = 1, size_2
                A(offset_1 + i, offset_2 + j) = &
                        A(offset_1 + i, offset_2 + j) + B((i - 1) * size_2 + j)
            end do
        end do

    end subroutine update_matrix_part

    subroutine print_matrix(file_name, A, size_1, size_2)

        implicit none
        character(len=*), intent(in) :: file_name
        real(8), intent(in) :: A(:,:)
        integer, intent(in) :: size_1, size_2

        integer :: i, j
        integer :: out
        real(8), parameter :: THRESH = 1d-10

        out = 12
        open(unit=out, file=file_name)

        write (out,'(i6)') size_1

        do i = 1, size_1
            do j = 1, i
                if (abs(A(i,j)) > THRESH) then
                    write(out,'(2i6,e30.16)') i, j, A(i,j)
                end if
            end do
        end do

        close(unit=out)

    end subroutine print_matrix

end module evalints

