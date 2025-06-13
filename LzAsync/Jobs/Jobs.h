#ifndef __JOBS_LIBS__
#define __JOBS_LIBS__

// ======================================================================================
/**
 * @file Jobs.h
 * @brief Unified job management system for task scheduling and tracking
 * @details
 * Provides core job abstraction used throughout the dispatch queue system,
 * encapsulating executable units with complete lifecycle tracking and
 * hierarchical relationships.
 */

#include <chrono>
#include <string>
#include <functional>

using std::string;

// ======================================================================================

/**
 * @brief Fundamental executable job unit
 * @details
 * Encapsulates work items with:
 * - Unique identification and hierarchical relationships
 * - Precise timing instrumentation
 * - Execution status tracking
 * - Callback function management
 * 
 * Forms the atomic unit of work in the dispatch queue system.
 */
class Jobs 
{
public:
    using time_point = std::chrono::time_point<std::chrono::steady_clock>;     ///< High-resolution timestamp type for precise execution tracking

    /**
     * @brief Constructs an executable job unit
     * @param pnameJob Descriptive identifier for the job
     * @param pFunc Callable unit of work to execute
     * @param pParentJobId Hierarchical parent identifier (-1 for root jobs)
     * @details
     * Automatically assigns unique sequential ID and maintains
     * parent-child relationships for job dependency tracking.
     */
    Jobs(const string& pnameJob, std::function<void(void)> pFunc, 
         int64_t pParentJobId = -1);
    

    /**
     * @brief Constructs a termination marker job
     * @param pTime Timestamp of termination event
     * @param pT0 System reference timestamp
     * @param pNameJob Description of termination cause
     * @param pError Optional error information
     */
    Jobs(time_point pTime, time_point pT0, const string& pNameJob = "Quitting", const string& pError = "unamed");

    
    /**
     * @brief Records job initiation metrics
     * @param pT0 System reference timestamp
     * @param pStart Activation timestamp
     */
    void JobStart(const time_point& pT0, time_point pStart);
    

    /**
     * @brief Records job completion metrics
     * @param pT0 System reference timestamp
     * @param pStop Completion timestamp
     * @param pSucces Execution status
     * @param pError Optional failure description
     */
    void JobDone(const time_point& pT0, time_point pStop, bool pSucces, const std::string& pError = "");


    /**
     * @brief Retrieves job's descriptive identifier
     * @return Human-readable job name
     */
    std::string GetJobName() const {return mJobName;};

    /**
     * @brief Gets total jobs created in system
     * @return Global job counter value
     */
    int64_t GetNbJob() const {return smNbJob;};

    /**
     * @brief Gets unique job identifier
     * @return This job's unique ID
     */
    int64_t GetJobId() const {return mJobId;};

    /**
     * @brief Gets hierarchical parent identifier
     * @return Parent job ID or -1 if root job
     */
    int64_t GetJobIdParent() const {return mParentJobId;};


    /**
     * @brief Gets relative start time
     * @return Milliseconds since reference epoch
     */
    int64_t GetStartmsec() const {return mStart_msec;};

    /**
     * @brief Gets relative completion time
     * @return Milliseconds since reference epoch
     */
    int64_t GetStopmsec() const {return mStop_msec;};


    /**
     * @brief Checks execution status
     * @return True if job completed successfully
     */
    bool IsSuccess() const {return mSucces;};

    /**
     * @brief Retrieves failure information
     * @return Description of failure condition
     */
    std::string GetError() const {return mError;};

    /**
     * @brief Gets executable unit
     * @return Callable job function
     */
    std::function<void(void)> GetFunc(){return mFunc;};

    /**
     * @brief Resets global job counter
     *
     */
    static void ResetJobCounter() { smNbJob = 0; }

private:

    std::string mJobName;               ///< Human-readable job identifier
    static int64_t smNbJob;             ///< Global job counter/ID generator
    int64_t mJobId;                     ///< Unique instance identifier
    int64_t mParentJobId;               ///< Parent job reference (-1 = none)

    int64_t mStart_msec;                ///< Relative start time (ms since epoch)
    int64_t mStop_msec;                 ///< Relative completion time (ms since epoch)

    bool mSucces;                       ///< Execution status flag
    std::string mError;                 ///< Failure description if !mSucces

    std::function<void(void)> mFunc;    ///< Encapsulated executable unit

};


#endif // __JOBS_LIBS__
