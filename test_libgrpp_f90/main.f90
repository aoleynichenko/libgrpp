!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2023 Alexander Oleynichenko
!

program test_libgrpp_f90

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
    integer :: i, j
    character(len=32), allocatable :: basis_labels(:)
    character(len=50) :: file_name
    character, parameter :: xyz_letters(3) = (/'x', 'y', 'z'/)
    real(8), allocatable :: nuc_attr_matrix(:,:)
    real(8), allocatable :: overlap_matrix(:,:)
    real(8), allocatable :: arep_matrix(:,:)
    real(8), allocatable :: esop_matrices(:,:,:)
    real(8), allocatable :: arep_grad(:,:,:,:)
    real(8), allocatable :: spin_orbit_grad(:,:,:,:,:)

    print *
    print *, '    -----------------------------------------------------'
    print *, '        the fortran front-end to the libgrpp library     '
    print *, '    -----------------------------------------------------'
    print *, '    a. oleynichenko                            9 feb 2023'
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

    ! calculate matrix elements of the GRPP operator
    allocate(arep_matrix(basis_dim,basis_dim))
    allocate(esop_matrices(3,basis_dim,basis_dim))
    arep_matrix = 0.0
    esop_matrices = 0.0

    call calculate_ecp_integrals(arep_matrix, esop_matrices)

    ! calculate gradients of GRPP matrix elements
    ! wrt coordinates of nuclei
    allocate(arep_grad(natoms,3,basis_dim,basis_dim))
    allocate(spin_orbit_grad(natoms,3,3,basis_dim,basis_dim))
    call calculate_rpp_integrals_gradients(arep_grad,spin_orbit_grad)

    ! overlap matrix
    allocate(overlap_matrix(basis_dim,basis_dim))
    overlap_matrix = 0.0

    call calculate_overlap_integrals(overlap_matrix)

    ! nuclear attraction integrals
    allocate(nuc_attr_matrix(basis_dim,basis_dim))
    nuc_attr_matrix = 0.0

    call calculate_nuclear_attraction_integrals(nuc_attr_matrix, LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE)

    ! flush matrices of AO integrals to files
    call print_matrix('libgrpp_f90_arep.txt', arep_matrix, basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_so_x.txt', esop_matrices(1,:,:), basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_so_y.txt', esop_matrices(2,:,:), basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_so_z.txt', esop_matrices(3,:,:), basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_overlap.txt', overlap_matrix, basis_dim, basis_dim)
    call print_matrix('libgrpp_f90_nucattr.txt', nuc_attr_matrix, basis_dim, basis_dim)

    ! gradients of AREP integrals - wrt nuclear coordinates
    do iatom = 1, natoms
        do i = 1, 3
            write (file_name,"('libgrpp_f90_arep_grad_',i0,a,'.txt')") iatom - 1, xyz_letters(i)
            call print_matrix(trim(file_name), arep_grad(iatom,i,:,:), basis_dim, basis_dim)
        end do
    end do

    ! gradients of SO integrals - wrt nuclear coordinates
    do iatom = 1, natoms
        do i = 1, 3
            do j = 1, 3
                write (file_name,"('libgrpp_f90_so_',a,'_grad_',i0,a,'.txt')") xyz_letters(j), iatom - 1, xyz_letters(i)
                call print_matrix(trim(file_name), spin_orbit_grad(iatom,i,j,:,:), basis_dim, basis_dim)
            end do
        end do
    end do

    ! cleanup
    deallocate(overlap_matrix)
    deallocate(nuc_attr_matrix)
    deallocate(arep_matrix)
    deallocate(esop_matrices)
    deallocate(basis_labels)

end program test_libgrpp_f90


