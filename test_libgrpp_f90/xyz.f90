!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2022 Alexander Oleynichenko
!

module xyz

    integer, parameter :: N_ELEMENTS = 131
    character(len = 8), parameter, dimension(N_ELEMENTS) :: element_symbols = &
            (/ 'Gh  ', 'H   ', 'He  ', 'Li  ', 'Be  ', 'B   ', 'C   ', 'N   ', &
               'O   ', 'F   ', 'Ne  ', 'Na  ', 'Mg  ', 'Al  ', 'Si  ', 'P   ', &
               'S   ', 'Cl  ', 'Ar  ', 'K   ', 'Ca  ', 'Sc  ', 'Ti  ', 'V   ', &
               'Cr  ', 'Mn  ', 'Fe  ', 'Co  ', 'Ni  ', 'Cu  ', 'Zn  ', 'Ga  ', &
               'Ge  ', 'As  ', 'Se  ', 'Br  ', 'Kr  ', 'Rb  ', 'Sr  ', 'Y   ', &
               'Zr  ', 'Nb  ', 'Mo  ', 'Tc  ', 'Ru  ', 'Rh  ', 'Pd  ', 'Ag  ', &
               'Cd  ', 'In  ', 'Sn  ', 'Sb  ', 'Te  ', 'I   ', 'Xe  ', 'Cs  ', &
               'Ba  ', 'La  ', 'Ce  ', 'Pr  ', 'Nd  ', 'Pm  ', 'Sm  ', 'Eu  ', &
               'Gd  ', 'Tb  ', 'Dy  ', 'Ho  ', 'Er  ', 'Tm  ', 'Yb  ', 'Lu  ', &
               'Hf  ', 'Ta  ', 'W   ', 'Re  ', 'Os  ', 'Ir  ', 'Pt  ', 'Au  ', &
               'Hg  ', 'Tl  ', 'Pb  ', 'Bi  ', 'Po  ', 'At  ', 'Rn  ', 'Fr  ', &
               'Ra  ', 'Ac  ', 'Th  ', 'Pa  ', 'U   ', 'Np  ', 'Pu  ', 'Am  ', &
               'Cm  ', 'Bk  ', 'Cf  ', 'Es  ', 'Fm  ', 'Md  ', 'Nb  ', 'Lr  ', &
               'Rf  ', 'Db  ', 'Sg  ', 'Bh  ', 'Hs  ', 'Mt  ', 'Ds  ', 'Rg  ', &
               'Cn  ', 'Nh  ', 'Fl  ', 'Mc  ', 'Lv  ', 'Ts  ', 'Og  ', 'E119', &
               'E120', 'E121', 'E122', 'E123', 'E124', 'E125', 'E126', 'E127', &
               'E128', 'E129', 'E130' /)

    ! mass numbers are given for the most abundant isotopes
    ! (or for the most long-lived for radioactive elements)

    integer, parameter, dimension(N_ELEMENTS) :: mass_numbers = &
            (/ 0  , 1  , 4  , 7  , 9  , 11 , 12 , 14 , &
               16 , 19 , 20 , 23 , 24 , 27 , 28 , 31 , &
               32 , 35 , 40 , 39 , 40 , 45 , 48 , 51 , &
               52 , 55 , 56 , 59 , 58 , 63 , 64 , 69 , &
               74 , 75 , 80 , 79 , 84 , 85 , 88 , 89 , &
               90 , 93 , 98 , 98 , 102, 100, 106, 107, &
               114, 115, 120, 121, 130, 127, 132, 133, &
               138, 139, 140, 141, 142, 145, 152, 153, &
               158, 159, 164, 165, 166, 169, 174, 175, &
               180, 181, 184, 187, 192, 193, 195, 197, &
               202, 205, 208, 209, 209, 210, 222, 223, &
               226, 227, 232, 231, 238, 237, 244, 243, &
               247, 247, 251, 252, 257, 258, 259, 266, &
               267, 268, 269, 270, 277, 278, 281, 282, &
               285, 286, 289, 290, 293, 294, 294, 303, &
               304, 305, 306, 307, 308, 309, 310, 311, &
               312, 313, 314 /)

    integer :: natoms
    integer, dimension(:), allocatable :: charges
    real(8), dimension(:, :), allocatable :: coord

contains
    subroutine read_molecule()

        implicit none

        ! local variables
        integer :: inp
        integer :: i
        character(len = 8) :: label
        real(8) :: x, y, z
        integer :: charge

        inp = 10
        open(unit = inp, file = 'molecule.xyz', status = 'old', err = 13)

        ! read number of atoms and allocate arrays
        read (inp, *) natoms
        allocate(charges(natoms))
        allocate(coord(natoms, 3))

        ! skip comment line
        read (inp, *)

        ! read x, y, z coordinates of atoms
        do i = 1, natoms
            read (inp, *, err = 14) label, coord(i, 1:3)
            call get_charge(label, charges(i))
        end do

        close(unit = inp)
        return

        13 continue
        print *, 'Input XYZ file with geometry is not found'
        stop
        14 continue
        print *, 'Syntax error in the XYZ file'
        stop

    end subroutine

    subroutine get_element_symbol(charge, symbol)

        implicit none
        character(len = *), intent(out) :: symbol
        integer, intent(in) :: charge

        symbol = trim(element_symbols(charge + 1))

    end subroutine get_element_symbol

    subroutine get_charge(label, charge)

        implicit none
        character(len = *), intent(in) :: label
        integer, intent(out) :: charge

        integer :: i

        do i = 1, N_ELEMENTS
            if (trim(element_symbols(i)) == trim(label)) then
                charge = i - 1
                return
            end if
        end do

        print *, 'Unknown element: ', trim(label)
        stop

    end subroutine get_charge

    subroutine print_molecule

        implicit none

        integer :: i

        print *
        print *, 'molecular geometry:'
        print *, '-------------------'
        do i = 1, natoms
            print '(i4,3f16.10)', charges(i), coord(i, 1:3)
        end do
        print *

    end subroutine print_molecule

end module xyz
