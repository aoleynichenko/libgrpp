!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2022 Alexander Oleynichenko
!

program new_molgep

    use xyz
    use basis
    use ecp
    use libgrpp
    use evalints
    implicit none

    ! local variables
    integer :: iatom
    integer :: z
    integer :: basis_dim
    integer :: i
    character(len=32), allocatable :: basis_labels(:)
    real(8), allocatable :: nuc_attr_matrix(:,:)
    real(8), allocatable :: arep_matrix(:,:)
    real(8), allocatable :: esop_matrices(:,:,:)

    print *
    print *, '    -----------------------------------------------------'
    print *, '        the fortran front-end to the libgrpp library     '
    print *, '    -----------------------------------------------------'
    print *, '    a. oleynichenko                           27 jun 2022'
    print *, '    -----------------------------------------------------'
    print *

    call read_molecule
    call print_molecule

    ! read basis set (and, optionally, ECP) for each atom
    do iatom = 1, natoms
        z = charges(iatom)
        call read_basis('basis.inp', z)
        call normalize_basis(z)
        call read_rpp('ecp.inp', z)
    end do

    ! print basis sets and ECPs
    do z = 1, N_ELEMENTS
        if (num_blocks(z) /= 0) then
            call print_basis(z)
        end if
        if (n_arep(z) /= 0) then
            call print_ecp(z)
        end if
    end do

    ! contracted basis functions should be normalized to 1
    do iatom = 1, natoms
        z = charges(iatom)
        call normalize_basis(z)
    end do

    ! calculate basis set size
    call cart_basis_size(natoms, charges, basis_dim)

    ! print atom-centered basis set
    allocate(basis_labels(basis_dim))
    call get_atom_centered_basis_set_labels(natoms, charges, basis_labels)
    print *
    print *, 'list of cartesian basis functions:'
    print *, '----------------------------------'
    do i = 1, basis_dim
        print *, basis_labels(i)
    end do
    print *

    ! calculate matrix elements of the ECP operator
    allocate(nuc_attr_matrix(basis_dim,basis_dim))
    allocate(arep_matrix(basis_dim,basis_dim))
    allocate(esop_matrices(3,basis_dim,basis_dim))

    nuc_attr_matrix = 0.0
    arep_matrix = 0.0
    esop_matrices = 0.0

    call calculate_ecp_integrals(arep_matrix, esop_matrices)
    !call calculate_nuclear_attraction_integrals(nuc_attr_matrix, LIBGRPP_NUCLEAR_MODEL_FERMI)

    ! flush nuclear attraction, AREP, ESOP matrices to files
    call print_matrix('libgrpp_f90_coulomb.txt', nuc_attr_matrix, basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_arep.txt', arep_matrix, basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_so_x.txt', esop_matrices(1,:,:), basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_so_y.txt', esop_matrices(2,:,:), basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_so_z.txt', esop_matrices(3,:,:), basis_dim, basis_dim)

    ! cleanup
    deallocate(nuc_attr_matrix)
    deallocate(arep_matrix)
    deallocate(esop_matrices)
    deallocate(basis_labels)

end program new_molgep


