#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GSL::gsl" for configuration ""
set_property(TARGET GSL::gsl APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(GSL::gsl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libgsl.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS GSL::gsl )
list(APPEND _IMPORT_CHECK_FILES_FOR_GSL::gsl "${_IMPORT_PREFIX}/lib/libgsl.a" )

# Import target "GSL::gslcblas" for configuration ""
set_property(TARGET GSL::gslcblas APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(GSL::gslcblas PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libgslcblas.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS GSL::gslcblas )
list(APPEND _IMPORT_CHECK_FILES_FOR_GSL::gslcblas "${_IMPORT_PREFIX}/lib/libgslcblas.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
