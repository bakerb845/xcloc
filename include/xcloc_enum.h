#ifndef XCLOC_ENUM_H__
#define XCLOC_ENUM_H__

enum xclocPrecision_enum
{
    XCLOC_SINGLE_PRECISION = 0, /*!< Single precision. */
    XCLOC_DOUBLE_PRECISION = 1  /*!< Double precision. */
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

#endif
