# set minimum cmake version
cmake_minimum_required(VERSION 3.18)

# project name and language
project(fmm2d LANGUAGES C Fortran)

# verbose makefile
set(CMAKE_VERBOSE_MAKEFILE ON)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Grab Python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
# Grab OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# build static lib libfmm2d.a
# source files for libfmm2d.a
file(GLOB_RECURSE source_list "src/*.f" "src/*.f90")
#message(${source_list})
# add fmm2d static lib
add_library(fmm2d STATIC ${source_list})
#compiler options
target_compile_options(fmm2d PRIVATE -fallow-argument-mismatch -fPIC -O3 -march=native -funroll-loops -std=legacy -w)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

# extensions variables
# hfmm2d_fortran
set(hfmm_module_name "hfmm2d_fortran")
set(hfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/helmholtz/hfmm2dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/helmholtz/hfmm2dwrap_vec.f
                          ${CMAKE_SOURCE_DIR}/src/helmholtz/helmkernels2d.f)
set(f2py_helm_module_c "${hfmm_module_name}module.c")
# lfmm2d_fortran
set(lfmm_module_name "lfmm2d_fortran")
set(lfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/laplace/lfmm2dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/lfmm2dwrap_vec.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/lapkernels2d.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/rfmm2dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/rfmm2dwrap_vec.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/rlapkernels2d.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/cfmm2dwrap.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/cfmm2dwrap_vec.f
                          ${CMAKE_SOURCE_DIR}/src/laplace/cauchykernels2d.f)
set(f2py_lap_module_c "${lfmm_module_name}module.c")
# bhfmm2d_fortran
set(bhfmm_module_name "bhfmm2d_fortran")
set(bhfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/biharmonic/bhfmm2dwrap.f
                           ${CMAKE_SOURCE_DIR}/src/biharmonic/bhkernels2d.f)
set(f2py_bh_module_c "${bhfmm_module_name}module.c")
# stfmm2d_fortran
set(stfmm_module_name "stfmm2d_fortran")
set(stfmm_fortran_src_file ${CMAKE_SOURCE_DIR}/src/stokes/stfmm2d.f
                           ${CMAKE_SOURCE_DIR}/src/stokes/stokkernels2d.f)
set(f2py_st_module_c "${stfmm_module_name}module.c")

# Generate extensions' sources
# hfmm2d_fortran
add_custom_target(
  hfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_helm_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_helm_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${hfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "hfmm2d_fortran"
                   ${hfmm_fortran_src_file}
                   --lower
  DEPENDS ${hfmm_fortran_src_file} # Fortran source
)
# lfmm2d_fortran
add_custom_target(
  lfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_lap_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_lap_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${lfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "lfmm2d_fortran"
                   ${lfmm_fortran_src_file}
                   --lower
  DEPENDS ${lfmm_fortran_src_file} # Fortran source
)
# bhfmm2d_fortran
add_custom_target(
  bhfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_bh_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_bh_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${bhfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "bhfmm2d_fortran"
                   ${bhfmm_fortran_src_file}
                   --lower
  DEPENDS ${bhfmm_fortran_src_file} # Fortran source
)
# stfmm2d_fortran
add_custom_target(
  stfmm_genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_st_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_st_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${stfmm_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "stfmm2d_fortran"
                   ${stfmm_fortran_src_file}
                   --lower
  DEPENDS ${stfmm_fortran_src_file} # Fortran source
)

# Set up extensions targets
# hfmm2d_fortran
Python_add_library(hfmm2d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_helm_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${hfmm_fortran_src_file}" # Fortran source(s)
)
# lfmm2d_fortran
Python_add_library(lfmm2d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_lap_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${lfmm_fortran_src_file}" # Fortran source(s)
)
# bhfmm2d_fortran
Python_add_library(bhfmm2d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_bh_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${bhfmm_fortran_src_file}" # Fortran source(s)
)
# stfmm2d_fortran
Python_add_library(stfmm2d_fortran MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_st_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${stfmm_fortran_src_file}" # Fortran source(s)
)

# Dependencies for extensions
# hfmm2d_fortran
target_link_libraries(hfmm2d_fortran PRIVATE Python::NumPy fmm2d)
add_dependencies(hfmm2d_fortran hfmm_genpyf)
target_include_directories(hfmm2d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")
# lfmm2d_fortran
target_link_libraries(lfmm2d_fortran PRIVATE Python::NumPy fmm2d)
add_dependencies(lfmm2d_fortran lfmm_genpyf)
target_include_directories(lfmm2d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")
# bhfmm2d_fortran
target_link_libraries(bhfmm2d_fortran PRIVATE Python::NumPy fmm2d)
add_dependencies(bhfmm2d_fortran bhfmm_genpyf)
target_include_directories(bhfmm2d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")
# stfmm2d_fortran
target_link_libraries(stfmm2d_fortran PRIVATE Python::NumPy fmm2d)
add_dependencies(stfmm2d_fortran stfmm_genpyf)
target_include_directories(stfmm2d_fortran PRIVATE "${F2PY_INCLUDE_DIR}")

add_subdirectory(python)
