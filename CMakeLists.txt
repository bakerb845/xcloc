cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(xcloc VERSION 2.0.0 LANGUAGES CXX)
enable_testing()

##########################################################################################
#                    Define the compiler standard and requisites                         #
##########################################################################################
include(FindOpenMP)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_C_FLAGS   " ${OpenMP_C_FLAGS}")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (XCLOC_WRAP_PYTHON)
   find_package(PythonInterp 3)
   find_package(PythonLibs 3 REQUIRED)
   message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   find_package(pybind11 REQUIRED)
endif()

find_package(GTest REQUIRED)

##########################################################################################
#                                       Include Directories                              #
##########################################################################################

include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   ${IPP_INCLUDE_DIR}
   ${MKL_INCLUDE_DIR}
)
set(PRIVATE_INCLUDE_DEPENDS ${IPP_INCLUDE_DIR} ${MKL_INCLUDE_DIR})
message("Include directories:" ${PRIVATE_INCLUDE_DEPENDS})

##########################################################################################
#                                        Sort Out Python Wrappers                        #
##########################################################################################
#set(RTSEIS_WRAP_PYTHON YES)
if (RTSEIS_WRAP_PYTHON)
   add_library(pyxcloc MODULE
               python/pyxcloc.cpp
               )
   target_link_libraries(pyrtseis PRIVATE pybind11::module rtseis)
endif()

##########################################################################################
#                                      Define Library Dependencies                       #
##########################################################################################

#configure_file(${CMAKE_SOURCE_DIR}/include/rtseis/config.h.in
#               ${CMAKE_SOURCE_DIR}/include/rtseis/config.h)
set(LIBSRC 
    src/correlationEngine.cpp
    src/correlationEngineParameters.cpp)

set(LIBALL ${MKL_LIBRARY} ${IPP_LIBRARY})
set(LIBALL_PYTHON xcloc ${LIBALL} ${PYTHON_LIBRARIES})

# cmake -DBUILD_SHARED_LIBS=YES /path/to/source
set(BUILD_SHARED_LIBS YES)
add_library(xcloc SHARED ${LIBSRC})
target_link_libraries(xcloc PRIVATE ${MKL_LIBRARY} ${IPP_LIBRARY})

##########################################################################################
#                                      Unit Tests                                        #
##########################################################################################
add_executable(testSerialLibrary
               testing/main.cpp
               testing/correlationEngine.cpp)
target_link_libraries(testSerialLibrary PRIVATE xcloc ${GTEST_BOTH_LIBRARIES})
target_include_directories(testSerialLibrary PRIVATE ${GTEST_INCLUDE_DIRS})
add_test(NAME testSerialLibrary
         COMMAND testSerialLibrary)

##########################################################################################
#                                      CPACK Packaging                                   #
##########################################################################################
set(CPACK_PACKAGE_NAME "Cross-correlation Location (xcloc)")
set(CPACK_PACKAGE_VENDOR "Ben Baker")
set(CPACK_PACKAGE_CONTACT "bakerb845@gmail.com")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A library for locating events with multilateration")
#set(CPACK_PACKAGE_VERSION_MAJOR ${XCLOC_VERSION_MAJOR})
#set(CPACK_PACKAGE_VERSION_MINOR ${XCLOC_VERSION_MINOR})
#set(CPACK_PACKAGE_VERSION_PATCH ${XCLOC_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

##########################################################################################
#                                  Configuration Details                                 #
##########################################################################################

