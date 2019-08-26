cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
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
set(FindIPP_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(GTest REQUIRED)
find_package(FindIPP REQUIRED)

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

#if (NOT RTSEIS)
#include(ExternalProject)
#set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)
#ExternalProject_Add(rtseis
#                    GIT_REPOSITORY https://github.com/uofuseismo/rtseis
#                    CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
#                               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
#                               -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
#                               -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
#                               -DIPP_INCLUDE_DIR=${IPP_INCLUDE_DIR}
#                               -DIPP_CORE_LIBRARY=${IPP_CORE_LIBRARY}
#                               -DIPP_VM_LIBRARY=${IPP_VM_LIBRARY}
#                               -DIPP_IPPS_LIBRARY=${IPP_IPPS_LIBRARY}
#                               -DMKL_INCLUDE_DIR=${MKL_INCLUDE_DIR}
#                               -DMKL_LIBRARY=${MKL_LIBRARY}
#                               -DGTEST_INCLUDE_DIRS=${GTEST_INCLUDE_DIRS}
#                               -DGTEST_BOTH_LIBRARIES=${GTEST_BOTH_LIBRARIES}
#                    )

#configure_file(${CMAKE_SOURCE_DIR}/include/rtseis/config.h.in
#               ${CMAKE_SOURCE_DIR}/include/rtseis/config.h)
set(LIBSRC 
    src/correlationEngine.cpp
    src/correlationEngineParameters.cpp
    src/diffractionStackMigrationEngine.cpp)

set(LIBALL ${MKL_LIBRARY} ${IPP_LIBRARY})
set(LIBALL_PYTHON xcloc ${LIBALL} ${PYTHON_LIBRARIES})

# cmake -DBUILD_SHARED_LIBS=YES /path/to/source
set(BUILD_SHARED_LIBS YES)
add_library(xcloc SHARED ${LIBSRC})
target_link_libraries(xcloc PRIVATE ${LIBALL})

##########################################################################################
#                                      Unit Tests                                        #
##########################################################################################
add_executable(testSerialLibrary
               testing/main.cpp
               testing/correlationEngine.cpp)
target_link_libraries(testSerialLibrary
                      PRIVATE xcloc ${MKL_LIBRARY} ${IPP_LIBRARY} ${GTEST_BOTH_LIBRARIES})
target_include_directories(testSerialLibrary
                           PRIVATE ${PRIVATE_INCLUDE_DEPENDS} ${GTEST_INCLUDE_DIRS})
add_test(NAME testSerialLibrary
         COMMAND testSerialLibrary)

##########################################################################################
#                                       Installation                                     #
##########################################################################################
include(GNUInstallDirs)
if (XCLOC_WRAP_PYTHON)
   install(TARGETS xcloc pyxcloc
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
   install(TARGETS xcloc
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

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


