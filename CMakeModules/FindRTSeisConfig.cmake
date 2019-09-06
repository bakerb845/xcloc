# Already in cache, be silent
if (RTSEIS_INCLUDE_DIR AND RTSEIS_LIBRARY)
   set (RTSEIS_FIND_QUIETLY TRUE)
   set (RTSEIS_FOUND TRUE)
endif()

find_path(RTSEIS_INCLUDE_DIR
          NAMES rtseis
          HINTS $ENV{RTSEIS_DIR}/include /usr/local/include)
find_library(RTSEIS_LIBRARY
             NAMES rtseis
             PATHS /usr/local/lib 
                   /usr/lib
                   /usr/local/lib64
                   /usr/lib64
                   $ENV{RTSEIS_DIR}/lib)
if (RTSEIS_INCLUDE_DIR AND RTSEIS_LIBRARY)
   set(RTSEIS_FOUND TRUE)
   include(FindPackageHandleStandardArgs)
   find_package_handle_standard_args(RTSEIS DEFAULT_MSG RTSEIS_LIBRARY RTSEIS_INCLUDE_DIR)
   mark_as_advanced(RTSEIS_INCLUDE_DIR RTSEIS_LIBRARY)
endif()
