!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2023 Alexander Oleynichenko
!

module ecp

    use xyz
    implicit none

    integer, parameter :: ECP_MAXL = 10
    integer, parameter :: ECP_MAX_CNTRCT_LEN = 50
    integer, parameter :: ECP_MAX_OC_SHELLS = 20

    ! effective core potential data
    integer :: n_arep(N_ELEMENTS)
    integer :: n_esop(N_ELEMENTS)
    integer :: n_oc_shells(N_ELEMENTS, ECP_MAXL)

    integer :: ecp_num_prim(N_ELEMENTS, ECP_MAXL)
    real(8) :: ecp_alpha(N_ELEMENTS, ECP_MAXL, ECP_MAX_CNTRCT_LEN)
    integer :: ecp_powers(N_ELEMENTS, ECP_MAXL, ECP_MAX_CNTRCT_LEN)
    real(8) :: arep(N_ELEMENTS, ECP_MAXL, ECP_MAX_CNTRCT_LEN)
    real(8) :: esop(N_ELEMENTS, ECP_MAXL, ECP_MAX_CNTRCT_LEN)
    real(8) :: ocpot(N_ELEMENTS, ECP_MAXL, ECP_MAX_OC_SHELLS, ECP_MAX_CNTRCT_LEN)

    integer :: ocbas_num_prim(N_ELEMENTS, ECP_MAXL)
    integer :: ocbas_num_contr(N_ELEMENTS, ECP_MAXL)
    real(8) :: ocbas_alpha(N_ELEMENTS, ECP_MAXL, ECP_MAX_CNTRCT_LEN)
    real(8) :: ocbas_coeffs(N_ELEMENTS, ECP_MAXL, ECP_MAX_CNTRCT_LEN, ECP_MAX_CNTRCT_LEN)

contains

    subroutine read_rpp(basis_file, charge)

        implicit none

        ! arguments
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
        integer :: max_oc
        integer :: n_arep_tmp
        integer :: n_oc_pot
        integer :: ioc
        integer :: irep

        n_arep(charge) = 0
        n_esop(charge) = 0
        n_oc_shells(charge, :) = 0

        inp = 10
        open(unit = inp, file = basis_file, status = 'old', err = 13)

        do while (.TRUE.)
            11 continue
            read (inp, *, end = 12, err = 11) star, z, max_oc, n_arep_tmp
            if (star == '*' .and. z == charge) then

                n_arep(z) = n_arep_tmp
                n_esop(z) = n_arep(z)

                do ioc = 1, max_oc
                    read (inp, *, err = 14) nprim, ncontr
                    ocbas_num_prim(z, ioc) = nprim
                    ocbas_num_contr(z, ioc) = ncontr
                    n_oc_shells(z, ioc) = ncontr
                    do iprim = 1, nprim
                        read (inp, *, err = 14) ocbas_alpha(z, ioc, iprim), ocbas_coeffs(z, ioc, iprim, 1:ncontr)
                    end do
                end do

                do irep = 1, n_arep(z)
                    read (inp, *, err = 14) nprim, n_oc_pot
                    ecp_num_prim(z, irep) = nprim
                    do iprim = 1, nprim
                        read (inp, *, err = 14) ecp_powers(z, irep, iprim), ecp_alpha(z, irep, iprim), &
                                arep(z, irep, iprim), esop(z, irep, iprim), &
                                (ocpot(z, irep, ioc, iprim), ioc = 1, n_oc_pot)
                    end do
                end do

                exit

            end if
        end do

        12 continue
        close(unit = inp)
        return

        13 continue
        print '(a,i3,a)', 'Input file with ECP for element ', charge, ' is not found'
        stop
        14 continue
        print *, 'Syntax error in the ECP file'
        stop

    end subroutine read_rpp

    subroutine print_ecp(charge)

        implicit none

        ! arguments
        integer, intent(in) :: charge

        character, parameter :: symbols_table(ECP_MAXL) = (/ 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M' /)
        integer :: iblock
        integer :: iprim
        integer :: irep
        integer :: ioc

        print *
        print '(a,i3,a)', ' generalized relativistic pseudopotential for element Z = ', charge, ':'
        print *, '-------------------------------------------------------------'

        print *, 'outercore basis set:'
        do irep = 1, n_arep(charge)
            if (n_oc_shells(charge, irep) == 0) then
                cycle
            end if
            print *, symbols_table(irep)
            do iprim = 1, ocbas_num_prim(charge, irep)
                print '(f20.10,a,20f16.10)', ocbas_alpha(charge, irep, iprim), '  ', &
                        ocbas_coeffs(charge, irep, iprim, 1:ocbas_num_contr(charge, irep))
            end do

        end do
        print *

        print *, 'semilocal rpp:'
        do irep = 1, n_arep(charge)
            print '(2x,a,31x,a,16x,a)', symbols_table(irep), 'arep', 'esop'
            do iprim = 1, ecp_num_prim(charge, irep)
                print '(i3,3f20.12)', ecp_powers(charge, irep, iprim), ecp_alpha(charge, irep, iprim), &
                        arep(charge, irep, iprim), esop(charge, irep, iprim)
            end do
        end do
        print *

        print *, 'outercore potentials:'
        do irep = 1, n_arep(charge)
            if (n_oc_shells(charge, irep) == 0) then
                cycle
            end if

            ! header line: LJ symbols
            if (irep == 1) then
                print '(34x,10(a,16x))', ('S1/2', ioc = 1, n_oc_shells(charge, irep))
            else
                write (*, '(35x)', advance = 'no')
                do ioc = 1, n_oc_shells(charge, irep), 2
                    write (*, '(2(a,i1,a,16x))', advance = 'no') &
                            symbols_table(irep), 2 * (irep - 1) - 1, '/2', &
                            symbols_table(irep), 2 * (irep - 1) + 1, '/2'
                end do
                print *
            end if

            do iprim = 1, ecp_num_prim(charge, irep)
                print '(i3,20f20.12)', ecp_powers(charge, irep, iprim), ecp_alpha(charge, irep, iprim), &
                        (ocpot(charge, irep, ioc, iprim), ioc = 1, n_oc_shells(charge, irep))
            end do
        end do
        print *

    end subroutine print_ecp

end module ecp