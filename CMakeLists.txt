cmake_minimum_required (VERSION 4.0)

# Build option for static linking (GNU only)
# to use: cmake -S . -B build -G "Ninja" -DBUILD_STATIC=ON
option(BUILD_STATIC "Link statically to avoid runtime DLLs (GNU only)" OFF)


project (kmc_tian
         DESCRIPTION "kMC_tian project"
         LANGUAGES Fortran)
# enable_language(Fortran)

# Produce a json file with compiler options 
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Disable Fortran preprocessing
set(CMAKE_Fortran_PREPROCESS OFF)

# Set output directories for binaries
if(BUILD_STATIC)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/debug_static)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/release_static)
else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/debug)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/release)
endif()

# Compiler options for gfortran
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -fimplicit-none -fdec-format-defaults -std=f2008")
    set(CMAKE_Fortran_FLAGS "-Wall -fmessage-length=0 ${dialect}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g3 -Og -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG")
endif()

# # Compiler options for Intel compiler
# if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#     set(dialect "-stand f08 -free -implicitnone")
#     set(bounds "-check bounds ")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -check all ${bounds} ${dialect}")
#     set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -no-save-temps -nofpp  ${dialect}")
# endif()

# Compiler flags for Intel (classic and oneAPI/LLVM front end)
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(dialect "-stand f08 -free -implicitnone")
    set(CMAKE_Fortran_FLAGS "${dialect}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -check all -check bounds")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -no-save-temps -nofpp")
endif()

# # Get information from GIT if available
# find_package(Git)
# if (GIT_FOUND)
#     execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty --broken
#                     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#                     OUTPUT_VARIABLE GIT_VERSION_STRING
#                     OUTPUT_STRIP_TRAILING_WHITESPACE)
#     execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
#                     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#                     OUTPUT_VARIABLE GIT_COMMIT_HASH
#                     OUTPUT_STRIP_TRAILING_WHITESPACE)
# else()
#     set(GIT_VERSION_STRING "N/A")
#     set(GIT_COMMIT_HASH    "N/A")
# endif()
# Generate version_info.f90 from template
# configure_file(
#     "${CMAKE_SOURCE_DIR}/version_info.f90.in"
#     "${CMAKE_BINARY_DIR}/version_info.f90"
#     @ONLY)



# Get information from GIT if available
find_package(Git)
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty --broken
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_VERSION_STRING
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_COMMIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    
    # Find the current branch ref file for dependency tracking
    execute_process(COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_HEAD_REF
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    
    if(GIT_HEAD_REF)
        string(REPLACE "refs/" "" GIT_HEAD_REF "${GIT_HEAD_REF}")
        set(GIT_REF_FILE "${CMAKE_SOURCE_DIR}/.git/${GIT_HEAD_REF}")
    endif()
else()
    set(GIT_VERSION_STRING "N/A")
    set(GIT_COMMIT_HASH    "N/A")
endif()

# Initial generation of version_info.f90 from template (for configure time)
configure_file(
    "${CMAKE_SOURCE_DIR}/version_info.f90.in"
    "${CMAKE_BINARY_DIR}/version_info.f90"
    @ONLY)

# Make version_info.f90 depend on Git state files so it regenerates on commits/changes
if(GIT_FOUND)
    set(GIT_DEPENDENCIES 
        "${CMAKE_SOURCE_DIR}/.git/HEAD"
        "${CMAKE_SOURCE_DIR}/.git/index")
    
    if(EXISTS "${GIT_REF_FILE}")
        list(APPEND GIT_DEPENDENCIES "${GIT_REF_FILE}")
    endif()
    
    add_custom_command(
        OUTPUT "${CMAKE_BINARY_DIR}/version_info.f90"
        COMMAND ${CMAKE_COMMAND} -E copy 
            "${CMAKE_SOURCE_DIR}/version_info.f90.in" 
            "${CMAKE_BINARY_DIR}/version_info.f90.in.tmp"
        COMMAND ${CMAKE_COMMAND} 
            -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DINPUT_FILE=${CMAKE_BINARY_DIR}/version_info.f90.in.tmp
            -DOUTPUT_FILE=${CMAKE_BINARY_DIR}/version_info.f90
            -P "${CMAKE_SOURCE_DIR}/cmake/UpdateVersion.cmake"
        DEPENDS 
            "${CMAKE_SOURCE_DIR}/version_info.f90.in"
            ${GIT_DEPENDENCIES}
        COMMENT "Updating version information from Git"
    )
endif()





# Gather sources but exclude any version_info.f90 in the source tree
file(GLOB_RECURSE sources RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/*.f90)
list(FILTER sources EXCLUDE REGEX "version_info\\.f90$")

# Add executable and explicitly include generated version_info.f90 file (from build dir)
add_executable(kmc_tian ${sources} ${CMAKE_BINARY_DIR}/version_info.f90)

# Static linking option for GNU compiler
if(BUILD_STATIC AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_link_options(kmc_tian PRIVATE
        "-static-libgcc"
        "-static-libgfortran"
        # Uncomment the next line for fully static linking (may fail if static system libs unavailable)
        # "-static"
    )
endif()

# Make generated modules visible to sources
target_include_directories(kmc_tian PUBLIC "${CMAKE_BINARY_DIR}")
