cmake_minimum_required (VERSION 4.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)	# produce a json file with compiler options

project (kmc_tian
		 DESCRIPTION "kMC_tian project"
		 LANGUAGES Fortran)
# enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/release)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

	add_compile_options(
	
	# options for all builds
    -Wall -ffree-form -fimplicit-none  -fdec-format-defaults
 
	# options for Debug builds
	$<$<CONFIG:DEBUG>:-g3>          	# Add full debug symbols
	$<$<CONFIG:DEBUG>:-fcheck=all>		
	$<$<CONFIG:DEBUG>:-fcheck=bounds>  	# check bounds
    $<$<CONFIG:DEBUG>:-Og>        		# Don't optimizeor debugging
	
	# Options for Release builds
    $<$<CONFIG:RELEASE>:-O2>    		# Optimize for performance
    $<$<CONFIG:RELEASE>:-DNDEBUG>      	# Define NDEBUG to disable asserts
	
	)
 endif()
 
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -check all ${bounds} ${dialect}")
    set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -no-save-temps -nofpp  ${dialect}")
endif()

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
else()
    set(GIT_VERSION_STRING "N/A")
    set(GIT_COMMIT_HASH    "N/A")
endif()

# Add the version string to the version.f90.in file and make version.f90
configure_file(
	"${CMAKE_SOURCE_DIR}/version_info.f90.in"
    "${CMAKE_BINARY_DIR}/version_info.f90"
    @ONLY)

#
# Compile.
#
file(GLOB_RECURSE sources  *.f90)
add_executable(kmc_tian	${sources} ${CMAKE_CURRENT_BINARY_DIR/version.f90} )

# Add the binary directory to include paths for Fortran
target_include_directories(kmc_tian PUBLIC "${CMAKE_BINARY_DIR}")


