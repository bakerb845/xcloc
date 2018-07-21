#ifndef XCLOC_ENUM_H__
#define XCLOC_ENUM_H__
#include "xcloc_config.h"

enum xclocNumbering_enum
{
    XCLOC_C_NUMBERING = 0,       /*!< C based indexing. */
    XCLOC_FORTRAN_NUMBERING = 1  /*!< Fortran based indexing. */
};

enum xclocPrecision_enum
{
    XCLOC_SINGLE_PRECISION = 0, /*!< Single precision. */
    XCLOC_DOUBLE_PRECISION = 1  /*!< Double precision. */
};

enum xclocMigrateXC_enum
{
    XCLOC_MIGRATE_PHASE_XCS = 0, /*!< Migrate phase correlograms. */
    XCLOC_MIGRATE_XCS = 1        /*!< Migrate cross correlograms. */
};

enum xclocAccuracy_enum
{
    XCLOC_HIGH_ACCURACY = 0,     /*!< High accuracy.  Expect the maximum
                                      allowable precision in exchange
                                      for the lowest performance. */
    XCLOC_MEDIUM_ACCURACY  = 1,  /*!< Expect to lose about 1 to 2
                                      bits of precision in exchange for
                                      better performacne than high
                                      accuracy mode. */
    XCLOC_EXTENDED_ACCURACY = 2  /*!< Expect to lose about half
                                      the number of bits worth of
                                      precision in exchange for the
                                      highest performance. */
};

enum xclocXCFilterType_enum
{
    XCLOC_SPXC_DONOT_FILTER = 0,    /*!< Do not filter the correlograms. */
    XCLOC_SPXC_ENVELOPE_FILTER = 1, /*!< Compute envelope of correlograms. */
    XCLOC_SPXC_RMS_FILTER = 2 /*!< Compute RMS Filter of correlograms. */
};

enum xclocXCVerbose_enum
{
    XCLOC_PRINT_NONE =-1,     /*!< Print nothing. */
    XCLOC_PRINT_ERRORS = 0,   /*!< Print errors. */
    XCLOC_PRINT_WARNINGS = 1, /*!< Print errors and warnings. */
    XCLOC_PRINT_INFO = 2,     /*!< Print errors, warnings, and general information. */
    XCLOC_PRINT_DEBUG = 3     /*!< Print everything. */
};



#endif
