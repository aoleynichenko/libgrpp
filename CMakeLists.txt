#
#  libgrpp - a library for the evaluation of integrals over
#            generalized relativistic pseudopotentials.
#
#  Copyright (C) 2021-2023 Alexander Oleynichenko
#

cmake_minimum_required(VERSION 3.19)

include(CMakeDependentOption)
include(CMakePackageConfigHelpers)
 
set(VERSION_MAJOR 2023)
set(VERSION_MINOR 12)
set(VERSION_PATCH 14)

project(libgrpp
        DESCRIPTION "A library for the evaluation of molecular integrals of the generalized relativistic pseudopotential operator (GRPP) over Gaussian functions."
        HOMEPAGE_URL https://github.com/aoleynichenko/libgrpp
        VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
        LANGUAGES Fortran C)

set(CMAKE_C_STANDARD 11)

# imply -O3 -DNDEBUG, can be changed with cmake -DCMAKE_BUILD_TYPE=DEBUG,etc....

set(CMAKE_BUILD_TYPE Release)

set(libgrpp_APIVERSION ${libgrpp_VERSION_MAJOR}.${libgrpp_VERSION_MINOR})

add_subdirectory(libgrpp)

#add_compile_options(-O3)
#add_compile_options(-pg)


add_executable(test_libgrpp_c.x
        test_libgrpp_c/main.c
        test_libgrpp_c/molecule.c
        test_libgrpp_c/basis.c
        test_libgrpp_c/rpp.c
        test_libgrpp_c/elements.c
        test_libgrpp_c/shell_list.c
        test_libgrpp_c/eval_integrals.c
        test_libgrpp_c/eval_integrals_grad.c
        test_libgrpp_c/print_matrix.c
        test_libgrpp_c/abs_time.c
        )

add_executable(test_libgrpp_f90.x
        test_libgrpp_f90/main.f90
        test_libgrpp_f90/xyz.f90
        test_libgrpp_f90/basis.f90
        test_libgrpp_f90/rpp.f90
        test_libgrpp_f90/evalints.f90
        test_libgrpp_f90/libgrpp.f90
        )

target_link_libraries(test_libgrpp_c.x grpp -lm)# -pg)
target_link_libraries(test_libgrpp_f90.x grpp -lm)


# testing binary
add_executable(libgrpp_compare_matrices.x test/compare_matrices.c)

# enable testing functionality
enable_testing()

# define tests

add_test(NAME Cs     WORKING_DIRECTORY ../test/Cs     COMMAND bash run_test.sh)
add_test(NAME K      WORKING_DIRECTORY ../test/K      COMMAND bash run_test.sh)
#add_test(NAME Tl_qz  WORKING_DIRECTORY ../test/Tl_qz  COMMAND bash run_test.sh)
add_test(NAME H2O    WORKING_DIRECTORY ../test/H2O    COMMAND bash run_test.sh)
add_test(NAME ThO    WORKING_DIRECTORY ../test/ThO    COMMAND bash run_test.sh)
#add_test(NAME UO2    WORKING_DIRECTORY ../test/UO2    COMMAND bash run_test.sh)

foreach (t Cs  K H2O ThO)
    set_property(TEST ${t} PROPERTY ENVIRONMENT "PATH=${CMAKE_BINARY_DIR}:$ENV{PATH}")
endforeach ()

include(GNUInstallDirs)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/grppConfigVersion.cmake"
  VERSION "${libgrpp_VERSION}"
  COMPATIBILITY SameMajorVersion)

configure_file("${PROJECT_SOURCE_DIR}/cmake/grppConfig.cmake.in"
               "${PROJECT_BINARY_DIR}/grppConfig.cmake" @ONLY)
configure_file(cmake/libgrpp.pc.in libgrpp.pc @ONLY)

install(FILES "${PROJECT_BINARY_DIR}/grppConfig.cmake"
              "${PROJECT_BINARY_DIR}/grppConfigVersion.cmake"
	      DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

install(FILES "${PROJECT_BINARY_DIR}/libgrpp.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

install(DIRECTORY libgrpp
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        PATTERN "*.h")

