#ifndef XCLOC_DIFFRACTIONSTACKMIGRATIONENGINE_HPP
#define XCLOC_DIFFRACTIONSTACKMIGRATIONENGINE_HPP 1
#include <memory>
#include "xcloc/enums.hpp"

namespace XCLoc
{
template<class T> class CorrelationEngine;
/*!
 * @class DiffractionStackMigrationEngine "diffractionStackMigration.hpp" "xcloc/diffractionStackMigration.hpp"
 * @brief This computes the diffraction stack migration of the correlograms.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T=double>
class DiffractionStackMigrationEngine
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    DiffractionStackMigrationEngine();
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~DiffractionStackMigrationEngine();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Sets a shared pointer to the correlogram engine.
     * @param[in] correlograms  A pointer to the correlogram engine.
     */
    void setCorrelationEngine(std::shared_ptr<const CorrelationEngine<T>> &correlograms);

    /*!
     * @brief Computes the diffraction stack migration of the correlograms.
     */
    void compute();
    /*!
     * @brief Indicates whether or not the correlogram engine has been set.
     * @result True indicates that the correlogram engine has been set.
     */
    bool haveCorrelationEngine() const noexcept;
private:
    class DSMImpl;
    std::unique_ptr<DSMImpl> pImpl; 
};
}
#endif
