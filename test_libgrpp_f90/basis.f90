!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2023 Alexander Oleynichenko
!

module basis

    use xyz
    implicit none

    integer, parameter :: MAXL = 10
    integer, parameter :: MAX_CNTRCT_LEN = 30

    ! basis set data
    integer :: num_blocks(N_ELEMENTS)
    integer :: block_num_prim(N_ELEMENTS, MAXL)
    integer :: block_num_contr(N_ELEMENTS, MAXL)
    real(8) :: exponents(N_ELEMENTS, MAXL, MAX_CNTRCT_LEN)
    real(8) :: coeffs(N_ELEMENTS, MAXL, MAX_CNTRCT_LEN, MAX_CNTRCT_LEN)

contains

    subroutine read_basis(basis_file, charge)

        implicit none
        character(len = *), intent(in) :: basis_file
        integer, intent(in) :: charge

        ! local variables
        integer :: inp
        integer :: i
        integer :: z
        integer :: nblocks
        character :: star
        integer :: iblock
        integer :: nprim
        integer :: ncontr
        integer :: iprim
        integer :: icontr
        real(8) :: norm_factor

        z = charge

        inp = 10
        open(unit = inp, file = basis_file, status = 'old', err = 13)

        do while (.TRUE.)
            11 continue
            read (inp, *, end = 12, err = 11) star, z, nblocks
            if (star == '*' .and. z == charge) then
                num_blocks(z) = nblocks
                do iblock = 1, nblocks
                    read (inp, *, err = 14) nprim, ncontr
                    block_num_prim(z, iblock) = nprim
                    block_num_contr(z, iblock) = ncontr
                    do iprim = 1, nprim
                        read (inp, *, err = 14) exponents(z, iblock, iprim), coeffs(z, iblock, 1:ncontr, iprim)
                    end do
                end do
            end if
        end do

        12 continue
        close(unit = inp)
        return

        13 continue
        print '(a,i3,a)', 'Input file with basis for element ', charge, ' is not found'
        stop
        14 continue
        print *, 'Syntax error in the basis file'
        stop

    end subroutine read_basis


    subroutine normalize_basis(charge)

        implicit none
        integer, intent(in) :: charge
        integer :: iblock
        integer :: iprim
        integer :: icontr
        real(8) :: norm_factor

        print *, 'normalization of basis functions for element Z = ', charge

        do iblock = 1, num_blocks(charge)
            do icontr = 1, block_num_contr(charge,iblock)

                call radial_gto_norm_factor(iblock - 1, block_num_prim(charge,iblock), &
                        coeffs(charge, iblock, icontr, :), &
                        exponents(charge, iblock, :), norm_factor)

                !print *, 'coeffs = ', coeffs(charge, iblock, icontr, :2)
                !print *, 'norm_factor = ', norm_factor

                coeffs(charge, iblock, icontr, :) = coeffs(charge, iblock, icontr, :) * norm_factor

            end do
        end do

    end subroutine normalize_basis


    subroutine print_basis(charge)

        implicit none
        integer, intent(in) :: charge

        character, parameter :: symbols_table(MAXL) = (/ 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M' /)
        integer :: iblock
        integer :: iprim

        print *
        print '(a,i3,a)', ' basis set for element Z = ', charge, ':'
        print *, '------------------------------'
        do iblock = 1, num_blocks(charge)
            print *, symbols_table(iblock)
            do iprim = 1, block_num_prim(charge,iblock)
                print '(f20.10,a,20f16.10)', exponents(charge, iblock, iprim), '  ', &
                        coeffs(charge, iblock, 1:block_num_contr(charge, iblock), iprim)
            end do
        end do
        print *

    end subroutine print_basis

    subroutine cart_basis_size(natoms, charges, size)

        implicit none

        integer, intent(in) :: natoms
        integer, intent(in) :: charges(*)
        integer, intent(out) :: size

        integer :: iatom
        integer :: z
        integer :: iblock
        integer :: L

        size = 0
        do iatom = 1, natoms
            z = charges(iatom)
            do iblock = 1, num_blocks(z)
                L = iblock - 1
                size = size + block_num_contr(z, iblock) * (L + 1) * (L + 2) / 2
            end do
        end do

    end subroutine cart_basis_size

    subroutine generate_cartesians(L, cart_list, ncart)

        implicit none
        integer, intent(in) :: L
        integer, dimension(:, :), intent(out) :: cart_list
        integer, intent(out) :: ncart

        integer :: r, s, t

        ncart = 0
        do r = L, 0, -1
            do s = L, 0, -1
                do t = L, 0, -1
                    if (r + s + t /= L) then
                        cycle
                    end if
                    ncart = ncart + 1
                    cart_list(ncart, 1) = r
                    cart_list(ncart, 2) = s
                    cart_list(ncart, 3) = t
                end do
            end do
        end do

    end subroutine generate_cartesians


    subroutine get_cartesian_fun_label(nlm, label)

        implicit none
        integer, intent(in) :: nlm(3)
        character(len = 16), intent(out) :: label

        integer :: i, j

        ! default: S function
        do i = 1, 16
            label(i:i) = ' '
        end do
        label(1:1) = 'S'

        i = 1
        do j = 1, nlm(1)
            label(i:i) = 'X'
            i = i + 1
        end do
        do j = 1, nlm(2)
            label(i:i) = 'Y'
            i = i + 1
        end do
        do j = 1, nlm(3)
            label(i:i) = 'Z'
            i = i + 1
        end do

    end subroutine get_cartesian_fun_label


    subroutine get_atom_centered_basis_set_labels(natoms, charges, labels)

        implicit none

        integer, intent(in) :: natoms
        integer, intent(in) :: charges(*)
        character(len = 32), intent(out) :: labels(*)

        integer :: size
        integer :: iatom
        integer :: z
        integer :: iblock
        integer :: L
        integer :: ncart
        integer :: cart_list(100, 3)
        integer :: icart
        integer :: ifun
        integer :: count
        character(len = 8) :: element_symbol
        character(len = 16) :: cart_label

        count = 0
        do iatom = 1, natoms
            z = charges(iatom)
            call get_element_symbol(z, element_symbol)
            do iblock = 1, num_blocks(z)

                call generate_cartesians(iblock - 1, cart_list, ncart)

                do ifun = 1, block_num_contr(z, iblock)
                    do icart = 1, ncart
                        count = count + 1
                        call get_cartesian_fun_label(cart_list(icart, :), cart_label)
                        write (labels(count), '(i3,a,a2,a,a)') count, ' ', element_symbol, ' ', cart_label
                    end do
                end do

            end do
        end do

    end subroutine get_atom_centered_basis_set_labels

end module basis