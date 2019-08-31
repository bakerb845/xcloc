#ifndef XCLOC_CORRELOGRAM_POSTPROCESSOR_HPP
#define XCLOC_CORRELOGRAM_POSTPROCESSOR_HPP
#include <memory>
namespace XCLoc
{
/// Forward declarations
class CorrelogramPostProcessorParameters;
template<class T> class CorrelationEngine;
template<class T>
class CorrelogramPostProcessor
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    CorrelogramPostProcessor();
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~CorrelogramPostProcessor();
    /*!
     * @brief Releases memory and resets all variables.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the post-processor.
     * @param[in] parameters  The post-processing parameters.
     * @param[in] engine      The the correlogram engine that will contain the
     *                        raw correlograms.
     * @throws std::invalid_argument if the parameters are not valid.
     */
    void initialize(const CorrelogramPostProcessorParameters &parameters,
                    std::shared_ptr<const XCLoc::CorrelationEngine<T>> engine);
    /*!
     * @brief Flag indicating that the post-processor is initialized.
     */
    bool isInitialized() const noexcept;
private:
    class CorrelogramPostProcessorImpl;
    std::unique_ptr<CorrelogramPostProcessorImpl> pImpl; 
};
}
#endif
