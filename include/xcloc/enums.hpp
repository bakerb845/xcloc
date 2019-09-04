#ifndef XCLOC_ENUMS_HPP
#define XCLOC_ENUMS_HPP
namespace XCLoc
{
/*!
 * @brief Defines the floating point precision.
 */
enum class Precision
{
    DOUBLE, /*!< Uses double (64 bit) precision calculations. */
    FLOAT   /*!< Uses float (32 bit) precision calcuations. */
};
/*!
 * @brief Defines the floating point accuracy used in some of the floating
 *        point computations in MKL.
 */
enum class MKLFloatingPointAccuracy
{
    HIGH_ACCURACY, /*!< The default high floating point accuracy. */ 
    LOW_ACCURACY,  /*!< This offers better perfformance than HIGH_ACCURACY.
                        For double precision this has about 12 digits of
                        accuracy.  For single precision this has about 4 digits
                        of accuracy. */
    EP_ACCURACY    /*!< This offers the best performance but for double 
                        precision this may have less than 8 digits of accuracy
                        and for float precision this may have less than four
                        digits of accuracy.  This mode is not recommended. */
};
/*!
 * @brief Defines the correlogram post-processing strategy.
 */
enum class CorrelogramFilteringType
{
    NO_FILTERING,           /*!< Do not modify raw correlograms. */
    FIR_ENVELOPE_FILTERING  /*!< Compute FIR-based envelope of correlograms.
                                 This helps to reduce side-lobe artifacts
                                 in the migration caused by imperfect
                                 velocity models. */
};

}
#endif
