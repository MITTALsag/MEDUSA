#ifndef __mDataTHREAD_H__
#define __mDataTHREAD_H__


// ======================================================================================
/**
 * @file DataThread.h
 * @brief Management system for thread-related data operations
 * @details
 * Provides comprehensive tracking of thread activities including:
 * - Job execution logs
 * - Mutex state transitions
 * - Tag lifecycle events
 * 
 * The system is designed for performance analysis and debugging of threaded applications.
 */

#include "../Jobs/Jobs.h"
#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <string>
#include <stdexcept>
#include <LzServices/LzLog.h>




// ======================================================================================
/**
 * @brief Mutex operation states
 * @details Tracks the lifecycle of mutex operations:
 * - NEED: Thread requests mutex access
 * - HAVE: Thread successfully acquires mutex
 * - FREE: Thread releases mutex
 */
enum class ThreadMutexState { NEED, HAVE, FREE };

/**
 * @brief Tag operation markers
 * @details Indicates tag boundary events:
 * - START: Beginning of tagged section
 * - STOP: End of tagged section
 */
enum class TagState { START, STOP };


// ======================================================================================
using std::string;
using std::array;
using std::map;
using std::vector;
using std::tuple;
using std::get;


// ======================================================================================
/**
 * @brief Job execution recorder
 * @details Maintains chronological record of job executions
 * with timing and status information for performance analysis.
 */

class JobsManager 
{
public:

    /**
     * @brief Adds a job record to the tracking system
     * @param pJobs The job object containing execution details to be recorded
     * @note The job is appended to the internal collection and persists until cleared
     */
    void AddLog(const Jobs& pJobs) { mJobs.push_back(pJobs); }

    /**
     * @brief Retrieves all tracked job records
     * @return const reference to the vector containing all job entries
     * @warning The returned reference becomes invalid if the manager is modified
     */
    const vector<Jobs>& GetJobs() const { return mJobs; }

private:
    vector<Jobs> mJobs; ///< Thread-safe collection of job execution records
};


// ======================================================================================
/**
 * @brief Mutex operation tracker
 * @details Records precise timing of mutex state transitions
 * with nanosecond resolution for contention analysis.
 * 
 * Stores data as [request_time, acquire_time, release_time] tuples
 * indexed by mutex name.
 */

class MutexContainer 
{
public:
    using MutexLogs = array<int64_t, 3>;  ///< Timestamp array tracking [request_time, acquisition_time, release_time]


    /**
     * @brief Retrieves complete mutex timing data
     * @return const reference to the mutex log database
     * @details 
     * Returns a read-only view of all collected mutex timing information,
     * organized by mutex name with associated timestamp triplets
     */
    const map<string, vector<MutexLogs>>& getData() const { return mData; }

    /**
     * @brief Records a mutex state transition timestamp
     * @param pTime The precise timestamp of the event (millisecond since DispatchQ::sT0)
     * @param pName The identifier of the mutex being tracked
     * @param pHow The transition type (NEED/HAVE/FREE)
     * @throws std::runtime_error if attempting to record a state transition without existing log
     * @throws std::invalid_argument if provided an invalid ThreadMutexState value
     * @note 
     * - Each mutex log contains three timestamps corresponding to:
     *   1. NEED (request time)
     *   2. HAVE (acquisition time) 
     *   3. FREE (release time)
     * - State transitions must follow NEED->HAVE->FREE sequence
     */
    void SetTimes(const int64_t pTime, const string& pName, const ThreadMutexState& pHow);

private:
    map<string, vector<MutexLogs>> mData; ///< Database of mutex operations indexed by name, storing timestamp sequences
};

// ======================================================================================
/**
 * @brief Tagged section monitor
 * @details Captures timing of marked code sections
 * for performance benchmarking and profiling.
 * 
 * Records [start_time, end_time] pairs associated with
 * each named tag for duration analysis.
 */
class TagContainer 
{
public:
    using Tag = tuple<string, array<int64_t, 2>>;  ///< Tag record containing name and [start_time, end_time] pair


    /**
     * @brief Retrieves all collected tag timing data
     * @return const reference to vector of tag records
     * @details
     * Returns chronological collection of all tag markers with their
     * associated timing information in [start, end] pairs
     */
    const vector<Tag>& getData() const { return mData; }

    /**
     * @brief Records a tag boundary event timestamp
     * @param pTime Nanosecond timestamp of the event
     * @param pName Identifier of the tagged section
     * @param pTagState Type of boundary (START/STOP)
     * @throws std::runtime_error if STOP received without matching START
     * @throws std::invalid_argument if invalid TagState provided
     * @note
     * - Tags must follow proper START-STOP pairing
     * - Timestamps should be monotonic for accurate duration calculation
     */
    void SetTimes(const int64_t pTime, const string& pName, const TagState& pTagState);

private:
    vector<Tag> mData; ///< Collection of tagged section records with timing data
};

// ======================================================================================
/**
 * @brief Thread instrumentation facade
 * @details Unified interface for comprehensive thread monitoring,
 * combining job tracking, mutex analysis, and section profiling.
 * 
 * Provides thread-safe aggregation of performance metrics
 * across all monitored thread activities.
 */
class DataThread 
{
public:
// ======================================================================================
    // Jobs

    /**
     * @brief Records a job execution entry
     * @param pLog Job record containing execution details
     */
    void AddLog(const Jobs& pLog) { mJobsManager.AddLog(pLog); }

    /**
     * @brief Retrieves all job execution records
     * @return const reference to collection of job records
     */
    const vector<Jobs>& GetJobs() const { return mJobsManager.GetJobs(); }


// ======================================================================================
    // MutexLogs

    /**
     * @brief Retrieves complete mutex timing database
     * @return const reference to mutex operation records
     * @details
     * Returns map of mutex names to their operation timelines,
     * with each entry containing [request, acquire, release] timestamps
     */
    const map<string, vector<MutexContainer::MutexLogs>>& getMutexLogs() const { return mMutexLogs.getData(); }

    /**
     * @brief Records a mutex state transition
     * @param pTime Event timestamp in miliseconds since DispatchQ::sT0
     * @param pName Mutex identifier
     * @param pHow Transition type (NEED/HAVE/FREE)
     * @see MutexContainer::SetTimes for detailed behavior
     */
    void setMutexTime(const int64_t pTime, const string& pName, const ThreadMutexState& pHow) {mMutexLogs.SetTimes(pTime, pName, pHow);}

// ======================================================================================
    // Tags

    /**
     * @brief Retrieves all tagged section records
     * @return const reference to tag database
     * @details
     * Returns collection of all tagged sections with their
     * [start, end] timing pairs
     */
    const vector<TagContainer::Tag>& getTags() const { return mTags.getData(); }


    /**
     * @brief Records a tag boundary event
     * @param pTime Event timestamp in milliseconds since DispatchQ::sT0
     * @param pName Tag identifier
     * @param pTagState Boundary type (START/STOP)
     * @see TagContainer::SetTimes for detailed behavior
     */
    void setTagsTime(const int64_t pTime, const string& pName, const TagState& pTagState) {mTags.SetTimes(pTime, pName, pTagState);}


private:
    JobsManager mJobsManager;       ///< Job execution recorder
    MutexContainer mMutexLogs;      ///< Mutex operation tracker
    TagContainer mTags;             ///< Section profiling system
};

#endif // __mDataTHREAD_H__