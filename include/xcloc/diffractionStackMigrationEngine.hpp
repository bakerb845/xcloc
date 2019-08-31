#ifndef XCLOC_DIFFRACTIONSTACKMIGRATIONENGINE_HPP
#define XCLOC_DIFFRACTIONSTACKMIGRATIONENGINE_HPP 1
#include <memory>
#include "xcloc/enums.hpp"

namespace XCLoc
{
template<class T> class Correlograms;
template<class T> class CorrelationEngine;
template<class T> class TravelTimeTable;
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

    /*! @name Travel time tables.
     * @{
     */
    void setSamplingRate(const double df);
    /*!
     * @brief Adds a travel time.
     * @param[in] waveformID  The waveform ID.
     * @param[in] phase       The phase label for the table, e.g., P or S.
     * @param[in] table       The travel time table. 
     */
    void addTravelTimeTable(int waveformID, const std::string &phase,
                            const TravelTimeTable<T> &table);
    /* @} */

    /*! @name Correlograms
     * @{
     */
    /*!
     * @brief Sets a shared pointer to the correlogram engine.
     * @param[in] correlograms  A pointer to the correlogram engine.
     */
    void setCorrelationEngine(std::shared_ptr<const CorrelationEngine<T>> correlograms);
    /*!
     * @brief Sets the correlograms that will be migrated.
     * @param[in] df            The sampling rate of the correlograms in Hz.
     * @param[in] correlograms  Pointer to the correlograms.
     * @note This uses dependency injection.  Hence, in a previous step the
     *       user would compute the (processed) correlograms and this class
     *       will immediately obtain the correlograms to migrate. 
     */
    void setCorrelograms(const double df,
                         std::shared_ptr<const Correlograms<T>> correlograms);
    /*!
     * @brief Determines if a pointer to the correlograms were set.
     * @result True indicates that the correlograms were set.
     */
    bool haveCorrelograms() const noexcept;
    /*! @} */

    /*!
     * @brief Create the travel time tables to be used in the migration. 
     * @
     */
    void createMigrationTables();

    /*!
     * @brief Computes the diffraction stack migration of the correlograms.
     */
    void compute();
    /*!
     * @brief Indicates whether or not the correlogram engine has been set.
     * @result True indicates that the correlogram engine has been set.
     */
    bool haveCorrelationEngine() const noexcept;
    /*! @} */
private:
    class DSMImpl;
    std::unique_ptr<DSMImpl> pImpl; 
};
}
#endif
