#ifndef XCLOC_ACOUSTICGREENS2D_HPP
#define XCLOC_ACOUSTICGREENS2D_HPP
#include <vector>
#include <memory>
class RickerWavelet;
/*!
 * @brief Computes the analytic Green's functions for a line-source in a 
 *        whole space.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class AcousticGreens2D
{
public:
    /*!
     * @brief Constructor.
     */
    AcousticGreens2D();
    /*!
     * @brief Destructor.
     */
    ~AcousticGreens2D();
    /*!
     * @brief Sets the source time function.
     * @param[in] stf  A Ricker wavelet source time function.
     * @throws std::invalid_argument if the source time function is invalid.
     */
    void setRickerSourceTimeFunction(const RickerWavelet &stf); 

    /*!
     * @brief Computes the Greens functions.
     */
    void compute();
private:
    class AcousticGreens2DImpl;
    std::unique_ptr<AcousticGreens2DImpl> pImpl;
};
#endif
