#ifndef XCLOC_CORRELOGRAMDIFFRACTIONSTACKMIGRATION_HPP
#define XCLOC_CORRELOGRAMDIFFRACTIONSTACKMIGRATION_HPP 1
#include <memory>
#include "xcloc/enums.hpp"

namespace XCLoc
{
class TravelTimeTableName;
template<class T> class Correlograms;
template<class T> class TravelTimeTable;
/*!
 * @class CorrelogramDiffractionStackMigration "correlogramDiffractionStackMigration.hpp" "xcloc/correlogramDiffractionStack.hpp"
 * @brief This computes the diffraction stack migration of the correlograms.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T=double>
class CorrelogramDiffractionStackMigration
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    CorrelogramDiffractionStackMigration();
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~CorrelogramDiffractionStackMigration();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Travel time tables.
     * @{
     */
    /*!
     * @brief Adds a travel time.
     * @param[in] tableName   The travel time table name.
     * @param[in] table       The travel time table. 
     * @throws std::invalid_argument if table name is invalid.
     */
    void addTravelTimeTable(const TravelTimeTableName &tableName,
                            const TravelTimeTable<T> &table);
    /* @} */

    /*! @name Correlograms
     * @{
     */
    /*!
     * @brief Sets the sampling rate of the correlograms.
     * @param[in] samplingRate  The sampling rate of the correlograms in Hz.
     * @throws std::invalid_argument if the sampling rate is not positive.
     */
    void setCorrelogramSamplingRate(double samplingRate);
    /*!
     * @brief Gets the sampling rate.
     * @result The sampling rate of the correlograms in Hz.
     * @throws std::runtime_error if \c setCorrelogramSamplingRate() was not
     *         called.
     * @sa \c haveCorrelogramSamplingRate()
     */
    double getCorrelogramSamplingRate() const;
    /*!
     * @brief Determines whether or not the correlogram sampling rate was set.
     * @result True indicates that the correlogram sampling rate was set.
     */
    bool haveCorrelogramSamplingRate() const noexcept;
    /*!
     * @brief Sets a shared pointer to the correlograms.
     * @param[in] correlograms  A pointer to the correlograms.
     */
    void setCorrelograms(std::shared_ptr<const Correlograms<T>> correlograms);
    /*!
     * @brief Determines if a pointer to the correlograms were set.
     * @result True indicates that the pointer to the correlograms were set.
     */
    bool haveCorrelogramPointer() const noexcept;
    /*! @} */

    /*!
     * @brief Create the travel time tables to be used in the migration. 
     * @throws std::runtime_error if the correlograms pointer and sampling rate
     *         has not been set or the travel time tables have not all be set.
     * @note This should be called only once after the correlogram information
     *       and travel time tables have been set.
     */
    void createMigrationTables();
    /*!
     * @brief Determines whether or not the travel time tables for the migration
     *        have been created.
     * @result True indicates that the travel time tables used in the migration
     *         have been created.
     */
    bool haveMigrationTables() const noexcept;

    /*!
     * @brief Computes the diffraction stack migration of the correlograms.
     * @throws std::runtime_error if the correlograms have not been computed in
     *         a previous step or the travel time models have not yet been set.
     */
    void compute();
private:
    class DSMImpl;
    std::unique_ptr<DSMImpl> pImpl; 
};
}
#endif
