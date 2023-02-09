!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2023 Alexander Oleynichenko
!

module libgrpp

    integer(4), parameter :: LIBGRPP_CART_ORDER_DIRAC = 0
    integer(4), parameter :: LIBGRPP_CART_ORDER_TURBOMOLE = 1

    integer(4), parameter :: LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE = 0
    integer(4), parameter :: LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL = 1
    integer(4), parameter :: LIBGRPP_NUCLEAR_MODEL_GAUSSIAN = 2
    integer(4), parameter :: LIBGRPP_NUCLEAR_MODEL_FERMI = 3
    integer(4), parameter :: LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE = 4

    interface

        subroutine libgrpp_set_default_parameters()
            ! no arguments
        end subroutine libgrpp_set_default_parameters


        subroutine libgrpp_set_radial_tolerance(tolerance)
            real(8), intent(in) :: tolerance
        end subroutine libgrpp_set_radial_tolerance


        subroutine libgrpp_set_angular_screening_tolerance(tolerance)
            real(8), intent(in) :: tolerance
        end subroutine libgrpp_set_angular_screening_tolerance


        subroutine libgrpp_set_modified_bessel_tolerance(tolerance)
            real(8), intent(in) :: tolerance
        end subroutine libgrpp_set_modified_bessel_tolerance


        subroutine libgrpp_set_cartesian_order(order)
            integer(4), intent(in) :: order
        end subroutine libgrpp_set_cartesian_order


        subroutine libgrpp_type1_integrals(&
                origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
                origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
                ecp_origin, ecp_num_primitives, ecp_powers, ecp_coeffs, ecp_alpha, &
                matrix &
                )
            ! shell centered on atom A
            real(8), dimension(*), intent(in) :: origin_A
            integer(4), intent(in) :: L_A
            integer(4), intent(in) :: num_primitives_A
            real(8), intent(in) :: coeffs_A(*)
            real(8), intent(in) :: alpha_A(*)
            ! shell centered on atom B
            real(8), dimension(*), intent(in) :: origin_B
            integer(4), intent(in) :: L_B
            integer(4), intent(in) :: num_primitives_B
            real(8), intent(in) :: coeffs_B(*)
            real(8), intent(in) :: alpha_B(*)
            ! effective core potential expansion
            real(8), dimension(*), intent(in) :: ecp_origin
            integer(4), dimension(*), intent(in) :: ecp_num_primitives
            integer(4), dimension(*), intent(in) :: ecp_powers
            real(8), dimension(*), intent(in) :: ecp_coeffs
            real(8), dimension(*), intent(in) :: ecp_alpha
            ! output: matrix with ECP integrals
            real(8), dimension(*), intent(out) :: matrix
        end subroutine libgrpp_type1_integrals


        subroutine libgrpp_type2_integrals(&
                origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
                origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
                ecp_origin, ecp_ang_momentum, ecp_num_primitives, ecp_powers, ecp_coeffs, ecp_alpha, &
                matrix &
                )
            ! shell centered on atom A
            real(8), dimension(*), intent(in) :: origin_A
            integer(4), intent(in) :: L_A
            integer(4), intent(in) :: num_primitives_A
            real(8), intent(in) :: coeffs_A(*)
            real(8), intent(in) :: alpha_A(*)
            ! shell centered on atom B
            real(8), dimension(*), intent(in) :: origin_B
            integer(4), intent(in) :: L_B
            integer(4), intent(in) :: num_primitives_B
            real(8), intent(in) :: coeffs_B(*)
            real(8), intent(in) :: alpha_B(*)
            ! effective core potential expansion
            real(8), dimension(*), intent(in) :: ecp_origin
            integer(4), intent(in) :: ecp_ang_momentum
            integer(4), dimension(*), intent(in) :: ecp_num_primitives
            integer(4), dimension(*), intent(in) :: ecp_powers
            real(8), dimension(*), intent(in) :: ecp_coeffs
            real(8), dimension(*), intent(in) :: ecp_alpha
            ! output: matrix with ECP integrals
            real(8), dimension(*), intent(out) :: matrix
        end subroutine libgrpp_type2_integrals


        subroutine libgrpp_spin_orbit_integrals(&
                origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
                origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
                ecp_origin, ecp_ang_momentum, ecp_num_primitives, ecp_powers, ecp_coeffs, ecp_alpha, &
                so_x_matrix, so_y_matrix, so_z_matrix &
                )
            ! shell centered on atom A
            real(8), dimension(*), intent(in) :: origin_A
            integer(4), intent(in) :: L_A
            integer(4), intent(in) :: num_primitives_A
            real(8), intent(in) :: coeffs_A(*)
            real(8), intent(in) :: alpha_A(*)
            ! shell centered on atom B
            real(8), dimension(*), intent(in) :: origin_B
            integer(4), intent(in) :: L_B
            integer(4), intent(in) :: num_primitives_B
            real(8), intent(in) :: coeffs_B(*)
            real(8), intent(in) :: alpha_B(*)
            ! effective core potential expansion
            real(8), dimension(*), intent(in) :: ecp_origin
            integer(4), intent(in) :: ecp_ang_momentum
            integer(4), dimension(*), intent(in) :: ecp_num_primitives
            integer(4), dimension(*), intent(in) :: ecp_powers
            real(8), dimension(*), intent(in) :: ecp_coeffs
            real(8), dimension(*), intent(in) :: ecp_alpha
            ! output: matrices with ECP integrals
            real(8), dimension(*), intent(out) :: so_x_matrix
            real(8), dimension(*), intent(out) :: so_y_matrix
            real(8), dimension(*), intent(out) :: so_z_matrix
        end subroutine libgrpp_spin_orbit_integrals

    end interface

contains

    subroutine libgrpp_outercore_potential_integrals(&
            origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
            origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
            ecp_origin, num_oc_shells, &
            oc_shells_L, oc_shells_J, ecp_num_primitives, ecp_powers, ecp_coeffs, ecp_alpha, &
            oc_shells_num_primitives, oc_shells_coeffs, oc_shells_alpha, &
            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix &
            )

        implicit none

        ! shell centered on atom A
        real(8), intent(in) :: origin_A(*)
        integer(4), intent(in) :: L_A
        integer(4), intent(in) :: num_primitives_A
        real(8), intent(in) :: coeffs_A(*)
        real(8), intent(in) :: alpha_A(*)
        ! shell centered on atom B
        real(8), intent(in) :: origin_B(*)
        integer(4), intent(in) :: L_B
        integer(4), intent(in) :: num_primitives_B
        real(8), intent(in) :: coeffs_B(*)
        real(8), intent(in) :: alpha_B(*)
        ! effective core potential expansion
        real(8), intent(in) :: ecp_origin(*)
        integer(4) :: num_oc_shells
        integer(4), intent(in) :: oc_shells_L(:)
        integer(4), intent(in) :: oc_shells_J(:)
        integer(4), intent(in) :: ecp_num_primitives(:)
        integer(4), intent(in) :: ecp_powers(:, :)
        real(8), intent(in) :: ecp_coeffs(:, :)
        real(8), intent(in) :: ecp_alpha(:, :)
        ! outercore shells
        integer(4) :: oc_shells_num_primitives(:)
        real(8) :: oc_shells_coeffs(:, :)
        real(8) :: oc_shells_alpha(:, :)
        ! output: matrices with ECP integrals
        real(8), intent(out) :: arep_matrix(*)
        real(8), intent(out) :: so_x_matrix(*)
        real(8), intent(out) :: so_y_matrix(*)
        real(8), intent(out) :: so_z_matrix(*)

        ! local variables
        integer :: ncart1
        integer :: ncart2
        integer :: i, j

        ncart1 = (L_A + 1) * (L_A + 2) / 2
        ncart2 = (L_B + 1) * (L_B + 2) / 2

        arep_matrix(1:ncart1 * ncart2) = 0.0d0
        so_x_matrix(1:ncart1 * ncart2) = 0.0d0
        so_y_matrix(1:ncart1 * ncart2) = 0.0d0
        so_z_matrix(1:ncart1 * ncart2) = 0.0d0

        ! the first non-local term:
        ! \sum_{nlj} U*|nlj><nlj| + |nlj><nlj|*U
        do i = 1, num_oc_shells
            call libgrpp_outercore_potential_integrals_part_1(&
                    origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
                    origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
                    ecp_origin, oc_shells_L(i), oc_shells_J(i), &
                    ecp_num_primitives(i), ecp_powers(i, :), ecp_coeffs(i, :), ecp_alpha(i, :), &
                    oc_shells_num_primitives(i), oc_shells_coeffs(i, :), oc_shells_alpha(i, :), &
                    arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix &
                    )
        end do

        ! the second non-local term:
        ! \sum_{nlj,n'lj} |nlj><nlj| U |n'lj><n'lj|
        do i = 1, num_oc_shells
            do j = 1, num_oc_shells

                call libgrpp_outercore_potential_integrals_part_2(&
                        origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
                        origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
                        ecp_origin, &
                        oc_shells_L(i), oc_shells_J(i), &
                        ecp_num_primitives(i), ecp_powers(i, :), ecp_coeffs(i, :), ecp_alpha(i, :), &
                        oc_shells_num_primitives(i), oc_shells_coeffs(i, :), oc_shells_alpha(i, :), &
                        oc_shells_L(j), oc_shells_J(j), &
                        ecp_num_primitives(j), ecp_powers(j, :), ecp_coeffs(j, :), ecp_alpha(j, :), &
                        oc_shells_num_primitives(j), oc_shells_coeffs(j, :), oc_shells_alpha(j, :), &
                        arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix &
                        )

            end do
        end do

    end subroutine libgrpp_outercore_potential_integrals

end module libgrpp

