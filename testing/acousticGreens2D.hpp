#ifndef XCLOC_ACOUSTICGREENS2D_HPP
#define XCLOC_ACOUSTICGREENS2D_HPP
#include <vector>
#include <memory>
class ISourceTimeFunction;
class RickerWavelet;
/*!
 * @brief Computes the analytic Green's functions for a line-source in a 
 *        whole space.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class AcousticGreens2D
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    AcousticGreens2D();
    /*!
     * @brief Copy constructor.
     * @param[in] greens  The 2D Green's function from which to initialize this.
     */
    AcousticGreens2D(const AcousticGreens2D &greens);
    /*!
     * @brief Move constructor.
     * @param[in,out] greens  The 2D Green's functions from whose memory will be
     *                        moved to this class.  On exit, greens' behavior
     *                        is undefined.
     */
    AcousticGreens2D(AcousticGreens2D &&greens) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] greens  The 2D acoustic Green's function to copy.
     * @result A deep copy of the greens.
     */
    AcousticGreens2D& operator=(const AcousticGreens2D &greens);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] greens  The 2D acoustic Green's functions whose memory
     *                        is to be moved to this.  On exit, greens' 
     *                        behavior is undefined. 
     * @result The memory from greens moved to this.
     */
    AcousticGreens2D& operator=(AcousticGreens2D &&greens) noexcept;
    /*! @} */
    
    /*!
     * @brief Destructor.
     */
    ~AcousticGreens2D();
    /*!
     * @brief Sets the source time function.
     * @param[in] stf  A Ricker wavelet source time function.
     * @throws std::invalid_argument if the source time function is invalid.
     */
    void setSourceTimeFunction(const ISourceTimeFunction &stf); 
    /*!
     * @brief Determines whether or not the source time function was set.
     * @result True indicates that the source time function was set.
     */
    bool haveSourceTimeFunction() const noexcept;

    /*!
     * @brief Sets the whose spaces' acoustic velocity.
     * @param[in] vel   The seismic velocity in m/s.
     * @throws std::invalid_argument if vel is not positive.
     */ 
    void setVelocity(double vel);
    /*!
     * @brief Sets the whose spaces' density.
     * @param[in] density  The density in kg/m**3
     * @throws std::invalid_argument if density is not positive.
     */
    void setDensity(double density);
    /*!
     * @brief Sets the quality factor.
     * @param[in] q   The quality factor.
     * @throws std::invalid_argument if q is not positive.
     */
    void setQualityFactor(double q);

    /*!
     * @brief Sets the receiver positions.
     * @param[in] positions  The (x,y) positions (meters) of the receivers.
     * @throws std::invalid_argument if positions.size() is 0.
     */
    void setReceiverPositions(const std::vector<std::pair<double, double>> &positions);
    /*!
     * @brief Sets the source position.
     * @param[in] position  The (x,y) position (meters) of the source.
     */
    void setSourcePosition(const std::pair<double, double> &position) noexcept;
                
    /*!
     * @brief Computes the Greens functions.
     * @throws std::runtime_error if the velocicity, density, source time 
     *         function, receiver positions, or source locations were not set.
     */
    void compute();
    /*!
     * @brief Gets the Green's function for the irec'th receiver.
     * @param[in] irec
     * @throws std::invalid_argument if irec is less than 0 or greater than or
     *         equal to the number of receivers.
     * @throws std::runtime_error if the Greens functions were not computed.
     */
    std::vector<double> getGreensFunction(int irec) const;
private:
    class AcousticGreens2DImpl;
    std::unique_ptr<AcousticGreens2DImpl> pImpl;
};
#endif
