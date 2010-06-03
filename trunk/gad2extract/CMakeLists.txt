# ============================================================================
# CMakeListst.txt file for GAD2EXTRACT PROJECT library                                        
# ============================================================================
cmake_minimum_required(VERSION 2.3)

# project name
project (GAD2EXTRACT)
#########################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
set (CMAKE_BUILD_TYPE Debug)
###################################
#SET(Boost_USE_STATIC_LIBS ON)
#SET(Boost_USE_MULTITHREAD ON)
SET( Boost_ADDITIONAL_VERSIONS  "1.43.0" "1.38")

find_package ( Boost  COMPONENTS program_options )
#thread )
################ TEST BOOST DIR 
MESSAGE (STATUS "BOOST FOUND IN " ${Boost_INCLUDE_DIRS})
MESSAGE (STATUS "BOOST LIB FOUND IN " ${Boost_LIBRARY_DIRS})
MESSAGE (STATUS "BOOST I WILL INCLUDE LIB  " ${Boost_LIBRARIES})
# contains the full path to the top level directory of your build tree 

##############################################
#find_package (GSL)
#find_package (GTS)
#MESSAGE (STATUS "GTS FOUND IN " ${GTS_INCLUDE_DIR})
#MESSAGE (STATUS "Here GTS has: " ${CMAKE_GTS_CXX_FLAGS})
##############################################
#include ("${CMAKE_CURRENT_SOURCE_DIR}/Modules/GenerateDocumentation.cmake")
#GENERATE_DOCUMENTATION("Doxyfile.in")


##############################################

MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
# Find all the sources for the LIB
FILE(GLOB LIBSOURCES lib/*.cpp)


# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH lib)


IF(NOT CMAKE_BUILD_TYPE OR 
		CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR 
		CMAKE_BUILD_TYPE STREQUAL "Debug" )
	SET(CMAKE_BUILD_TYPE  "DEBUG" CACHE STRING
		"Choose the type of build, options are: None=Debug"
		FORCE)
	SET(CMAKE_CXX_FLAGS "-O3  -g3 " )
	#	SET(CMAKE_VERBOSE_MAKEFILE yes)
ELSE(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE  "RELEASE" CACHE STRING
	                 "Choose the type of build, options are: None=Debug"
	                 FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE OR 
	CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR
	CMAKE_BUILD_TYPE STREQUAL "Debug")


if( CMAKE_BUILD_TYPE STREQUAL "RELEASE")
if(CMAKE_COMPILER_IS_GNUCXX)
	#	SET(CMAKE_CXX_FLAGS "-finline-functions -funroll-all-loops -Wall -msse3  -O3 -ffast-math -mfpmath=sse -Wno-non-virtual-dtor" )
SET(CMAKE_CXX_FLAGS "-O3  ")
else(CMAKE_COMPILER_IS_GNUCXX)
        SET(CMAKE_CXX_FLAGS "-O3 -xW -Zp16 -ip -ansi_alias -ansi  -mp -vec-report3 -ipo " )	
endif(CMAKE_COMPILER_IS_GNUCXX)

endif( CMAKE_BUILD_TYPE STREQUAL "RELEASE")

#################
# Fixing missing glib2
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
#${CMAKE_GTS_CXX_FLAGS}")

#################

MESSAGE(STATUS "BUILD TYPE: " ${CMAKE_BUILD_TYPE} )

# create library "armparallel"
add_library (armparallel SHARED ${LIBSOURCES})
#
# ------------------
# PROGRAMS          
# ------------------

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${GAD2EXTRACT_SOURCE_DIR}/src ${GAD2EXTRACT_SOURCE_DIR}/lib ${Boost_INCLUDE_DIRS} )

# Make sure the linker can find the Hello library once it is built.
link_directories (${GAD2EXTRACT_BINARY_DIR}/lib  ${Boost_LIBRARY_DIRS} )

# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.

# Find all LIB's sources files
FILE(GLOB execpp_sources src/*.cpp)

# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

  # Link the executable to the Hello library.
  target_link_libraries (${exe} armparallel ${Boost_LIBRARIES}  )

ENDFOREACH(exe_cpp ${execpp_sources})

# ----------------------------------------------------------
# install target                                            

set (CMAKE_INSTALL_PREFIX $ENV{LOCAL})

# lib
INSTALL(FILES  ${GAD2EXTRACT_BINARY_DIR}/lib/libarmparallel.so DESTINATION lib)

################################################################
