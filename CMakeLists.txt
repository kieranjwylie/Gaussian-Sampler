cmake_minimum_required(VERSION 3.15)
project(gaussian_sampler LANGUAGES Fortran)

# Add all Fortran source files in src/
file(GLOB SRC_FILES src/*.f90 src/*.f)

# Create the executable
add_executable(sampler ${SRC_FILES})

# (Optional) Set Fortran standard
set_target_properties(sampler PROPERTIES
    Fortran_STANDARD 2008
    Fortran_STANDARD_REQUIRED YES
)