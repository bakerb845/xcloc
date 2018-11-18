CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)
PROJECT(xcloc)
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)
ENABLE_TESTING()

#SET(CMAKE_C_COMPILER ${C_COMPILER}) #"icc")
#SET(CMAKE_CXX_COMPILER ${CXX_COMPILER}) #"icpc")

SET(XCLOC_VERSION_MAJOR 1)
SET(XCLOC_VERSION_MINOR 1)
SET(XCLOC_VERSION_PATCH 1)
SET(XCLOC_VERSION ${XCLOC_VERSION_MAJOR}.${XCLOC_VERSION_MINOR}.${XCLOC_VERSION_PATCH})
MESSAGE("Configuring version: " ${XCLOC_VERSION})
# Add modules for finding requisites
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/modules)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
SET(CMAKE_TESTING_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testing)
SET(CMAKE_BINARY_DIR               ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

INCLUDE(FindOpenMP QUIET)
IF (OPENMP_FOUND)
   MESSAGE("OpenMP Found")
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()

#ADD_SUBDIRECTORY(include)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/xcloc_config.h)
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${IPP_INCLUDE_DIR}
  ${MKL_INCLUDE_DIR}
  #${H5_C_INCLUDE_DIR}
  #${XML2_INCLUDE_DIR}
)
IF (XCLOC_USE_MPI)
   ADD_DEFINITIONS(-DXCLOC_USE_MPI_F)
   INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH}  ${MPI_C_INCLUDE_PATH})
ENDIF()
IF (XCLOC_PROFILE)
   INCLUDE_DIRECTORIES(${ADVISOR_INCLUDE_DIR})
ENDIF()

SET(LIB_ALL #${SACIO_LIBRARY} #${ISCL_LIBRARY}
    ${MKL_LIBRARY} ${IPP_LIBRARY})
    #${XML2_LIBRARY} ${H5_LIBRARY} ${INIPARSER_LIBRARY})
IF (XCLOC_USE_MPI)
    SET(LIB_ALL ${LIB_ALL} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
ENDIF()
SET(LIB_ALL ${LIB_ALL} -lpthread -ldl -lm)


#SET(SRC xcloc.c h5ioGrid.c h5ioUtils.c xdmf.c)
#SET(SRC ${SRC} )
#SET(IO_SRC h5ioGrid.c h5ioUtils.c xdmf.c)
SET(SRC ${IO_SRC}
    firFilter.c fdxc.F90 ipps.F90 constants.F90
    memory.F90 spxc.F90
    dsmxc.F90 utils.F90 sort.c)
IF (XCLOC_USE_MPI)
   SET(SRC ${SRC} utils_mpi.F90 fdxc_mpi.F90 dsmxc_mpi.F90 xcloc_mpi.F90)
ENDIF()
SET(SRC ${SRC} xcloc.F90)
SET(HEADERS include/xcloc_config.h include/xcloc.h include/xcloc_enum.h include/xcloc_finter.h)

# Executables
#ADD_EXECUTABLE(dales dales.c ${SRC})
# Unit tests
#SET(UNIT_TESTS_SRC unit_tests/acousticGreens2D.c)
ADD_LIBRARY(xcloc_shared SHARED ${SRC})
#ADD_LIBRARY(xcloc_static STATIC ${SRCS})
#ADD_EXECUTABLE(interloc interloc.c ${SRC})
SET(TEST_SRC
    testing/test_suite.c
    testing/serialTest.c
    testing/serialLocationTest.c
    testing/acousticGreens2D.c
    testing/xclocTest.c
    testing/iscl_hack.c)
IF (XCLOC_USE_MPI)
   SET(TEST_SRC ${TEST_SRC} testing/parallelTest.c testing/parxclocTest.c)
ENDIF()
ADD_EXECUTABLE(xcf_tests ${TEST_SRC})
#ADD_EXECUTABLE(xcfft unit_tests/check_xcfft.c ${UNIT_TESTS_SRC} ${SRC})
#ADD_EXECUTABLE(xc2d unit_tests/check_2d.c ${UNIT_TESTS_SRC} ${SRC})
#ADD_EXECUTABLE(xcloc unit_tests/xcloc.c ${UNIT_TESTS_SRC} ${SRC})
#ADD_EXECUTABLE(scaling unit_tests/scaling.c ${UNIT_TESTS_SRC} ${SRC})

TARGET_LINK_LIBRARIES(xcloc_shared ${LIB_ALL})
TARGET_LINK_LIBRARIES(xcf_tests xcloc_shared ${LIB_ALL})
#TARGET_LINK_LIBRARIES(interloc ${LIB_ALL})
#TARGET_LINK_LIBRARIES(xcfft ${LIB_ALL})
#TARGET_LINK_LIBRARIES(xc2d ${LIB_ALL})
#TARGET_LINK_LIBRARIES(xcloc ${LIB_ALL})
#TARGET_LINK_LIBRARIES(scaling ${LIB_ALL})

# Set the targets
SET_TARGET_PROPERTIES(xcloc_shared PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#SET_TARGET_PROPERTIES(xcloc_static PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#SET_TARGET_PROPERTIES(interloc  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) 
#SET_TARGET_PROPERTIES(xcfft     PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
#SET_TARGET_PROPERTIES(xc2d      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
#SET_TARGET_PROPERTIES(xcloc     PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
SET_TARGET_PROPERTIES(xcf_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})
#SET_TARGET_PROPERTIES(scaling   PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_TESTING_OUTPUT_DIRECTORY})

##########################################################################################
#                                       Unit Tests                                       #
##########################################################################################

#ADD_TEST(NAME xcfft
#         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
#         COMMAND xcfft)
ADD_TEST(NAME xcf_tests
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND xcf_tests)

##########################################################################################
#                                      CPACK Packaging                                   #
##########################################################################################
SET(CPACK_PACKAGE_NAME "Cross-correlation Location (xcloc)")
SET(CPACK_PACKAGE_VENDOR "Ben Baker")
SET(CPACK_PACKAGE_CONTACT "bakerb845@gmail.com")
SET(CPACK_PACKAGE_LICENSE "MIT")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A library for locating events with multilateration")
SET(CPACK_PACKAGE_VERSION_MAJOR ${XCLOC_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${XCLOC_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${XCLOC_VERSION_PATCH})
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

##########################################################################################
#                                  Configuration Details                                 #
##########################################################################################
SET(XCLOC_LIBRARY xcloc)
SET(INCLUDE_INSTALL_DIR include/)
SET(LIB_INSTALL_DIR lib/)
SET(SYSCONFIG_INSTALL_DIR etc/xcloc/)
INCLUDE(CMakePackageConfigHelpers)
CONFIGURE_PACKAGE_CONFIG_FILE(${CMAKE_MODULE_PATH}/FindXCLOC.cmake.in
                              ${CMAKE_BINARY_DIR}/FindXCLOC.cmake
                              INSTALL_DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/xcloc/cmake
                              PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/xclocConfigVersion.cmake
                                 VERSION ${XCLOC_VERSION} 
                                 COMPATIBILITY SameMajorVersion)



IF (NOT (${CMAKE_LIBRARY_OUTPUT_DIRECTORY} STREQUAL ${CMAKE_INSTALL_PREFIX}/lib))
   INSTALL(TARGETS xcloc_shared
           PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_WRITE GROUP_READ WORLD_EXECUTE WORLD_WRITE WORLD_READ
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include
           LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
   INSTALL(FILES ${HEADERS}
           DESTINATION ${CMAKE_INSTALL_PREFIX}/include/xcloc)
ENDIF()
#IF (NOT (${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} STREQUAL ${CMAKE_INSTALL_PREFIX}/lib))
#   INSTALL(TARGETS xcloc_static
#           ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib))
#ENDIF()
#INSTALL(TARGETS ${TARGET_NAMES}
#        RUNTIME DESTINATION ${CMAKE_BINARY_DIR})
INSTALL(FILES
        ${CMAKE_BINARY_DIR}/FindXCLOC.cmake
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/xclocConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/xcloc)


