CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/xcloc_config.h)
SET(headers xcloc_config.h xcloc.h xcloc_enum.h xcloc_envelope.h
    xcloc_migrate.h xcloc_migrateMPI.h
    xcloc_rmsFilter.h xcloc_xcfft.h xcloc_xcfftMPI.h)
