#pragma once

//======= Includes =======
#include "LzLib_Async.h"
#include "Jobs/Jobs.h"
#include "DataThread/DataThread.h"
#include "GraphCycleThreads/GraphCycleThreads.h"

#include <LzServices/List.h>

#include <thread>
#include <functional>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <chrono>
#include <string>
#include <iostream>
#include <map>
#include <unistd.h>
#include "json.hpp"

//======= DLL exports =======
namespace LzAsync
{
    class DLL_EXPORT_LzAsync DispatchQ;
}

//======= Namespace and Declarations ========
namespace LzAsync
{
using LzServices::List;
using std::string;
using std::vector;
using std::mutex;
using std::queue;
using std::pair;
using std::condition_variable;
using std::atomic_bool;
using std::map ;
using std::cout ;
using std::endl ;
using nlohmann::json_abi_v3_12_0::ordered_json;



// ======================================================================================
/**
 * @brief The mode of logging 
 * @details Tracks the lifecycle of mutex operations:
 * - OFF: No Logs
 * - END: The most precise logs, only at the end of the execution will not work if the program crashes.
 * - REALTIME: Logs are exported in real-time, can detect basics deadlocks (circular waiting mutexes), but not as precise as END mode.
 */
enum class LogsMod {OFF,END, REALTIME};

/**
 * @brief The DispatchQ class
 *
 * Log entries are enabled by default (@param mEnableLogs). 
 * Some marginal perf gain can be expected by disabling them at the expense of a lesser traceability.
 *
 * @warning This class deals with multiple worker threads but it
 * is not ITSELF thread-safe i.e., it is NOT SAFE to have multiple
 * thread accessing a shared DispatchQ object. As this would be an
 * exotic use-case, an API lock was not implemented and is missing
 * in the implementation to ensure thread-safety.
 *
 * The only exception being the "dispatch" method which can be invoked
 * by any thread, including the worker threads themselves.
 *
 */
class DispatchQ
{
public:
    // Thread types
    using func = std::function<void(void)>; // prec : job
    using func_and_name = std::pair<const string, func>; // prec : named_job


    /**
     * @brief DispatchQ
     * @param pName
     * @param pThread_cnt
     */
    DispatchQ( const string & pName="<unnamed>", size_t pThread_cnt=1 );

    /**
     * @brief ~DispatchQ
     * Destructor: wait for all threads to finish and free the resources
     * @warning This method is blocking and will wait for all threads to finish
     */
    virtual ~DispatchQ();
	
    /**
     * @brief getName
     * @return Name of the DispatchQ
     * @warning This method is not thread-safe. It should be called only from the main thread.
     */
    const string & getName() const { return mName; }

    /**
     * @brief startProfiling
     * @details modify the mEnableLogs variable to true and set the sT0 variable to the current time.
     * Also clears all the previous logs.
     * @warning This method is not thread-safe. It should be called only from the main thread.
     */

    //////////////////////// A MOFIfER LEs comms ///////////////////////////
    void startProfiling(const string& pToDir, const LogsMod pLogsMod = LogsMod::END);

    /**
     * @brief stopProfiling
     * @details modify the mEnableLogs variable to false.
     * @warning This method is not thread-safe. It should be called only from the main thread.
     */
    void stopProfiling() { mLogsMod = LogsMod::OFF; }  //Refaire plus tard

    /**
     * @brief Get the Nb Errors object
     * 
     * @return int 
     */
    int GetNbErrors() { return mNb_Errors; }

    /**
     * @brief Get the Nb Success object
     * 
     * @return int 
     */
    int GetNbSuccess() { return mNb_Success; }


    /**
     * @brief getLogEntries
     * @return vector<string> containing the log entries
     * @details Returns the log entries of the DispatchQ.
     * The log entries are stored in a vector of strings.
     * The log entries are cleared after this call.
     */
    vector<std::thread::id> GetThreadIds() const;


    /**
     * @brief waitForIdle
     * @param pRestartPool
     * Blocking wait for proper termination of all jobs in queue
     * Waits until the Q is empty AND no more workers are doing any work
     *
     * This makes it possible for workers to push jobs to the queue in order
     * to parallelize recursive processing.
     */
    void waitForIdle( bool pRestartPool ) { finish(mIdleQuit, pRestartPool); }

    /**
     * @brief completeJobs
     * @param pRestartPool
     * Blocking wait for proper termination of all jobs in queue
     * FINAL STATE: all jobs completed, no more jobs in queue, thread pool still running
     */
    void completeJobs( bool pRestartPool ) { finish(mEmptyQuit, pRestartPool); }

    /**
     * @brief interruptJobs
     * @param pRestartPool
     * @param pExportLogs To export logs to file if mLogsMod is not OFF (default: true) (and false just in the destructor)
     * Blocking wait for proper termination of all jobs in queue
     * FINAL STATE: jobs partially completed, no more jobs in queue, thread pool still running
     */
    void interruptJobs( bool pRestartPool, bool pExportLogs = true ) { finish(mEmptyQuit, pRestartPool, pExportLogs); }
    

    /**
     * @brief MutexNeed
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the mutex
     * @see StatusChangeMutex
     */
    friend void MutexNeed(DispatchQ& pDQ, const string& pName);

    /**
     * @brief MutexHave
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the mutex
     * @see StatusChangeMutex
     */
    friend void MutexHave(DispatchQ& pDQ, const string& pName);

    /**
     * @brief MutexFree
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the mutex
     * @see StatusChangeMutex
     */
    friend void MutexFree(DispatchQ& pDQ, const string& pName);


    /**
     * @brief TagStart
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the Tag
     * @see ChangeTags
     */
    friend void TagStart(DispatchQ& pDQ, const std::string& pName);

    /**
     * @brief TagStop
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the Tag
     * @see ChangeTags
     */
    friend void TagStop(DispatchQ& pDQ, const std::string& pName);


    /**
     * @brief restartPool
     * @param pNbWorkers
     * @param pToDir
     * @warning Will INTERRUPT all jobs, and terminate all current worker threads before restarting.
     * @details This method is blocking and will wait for all threads to finish
     * if mLogsMod is not OFF, the logs will be exported to mToFile file before restarting
     * The user have to call again startprofiling after this function to re-enable logs
     */
    void restartPool( size_t pNbWorkers);


    /**
     * @brief dispatch: Dispatch and copy
     * @param op
     */
    void dispatch( const func & pJob, const string & pName="<unnamed>" );
    void dispatch( const vector<func> & pJobs, const string & pName="<unnamed>" );
    void dispatch( const List<func> & pJobs, const string & pName="<unnamed>" );
    void dispatch( func_and_name && pNamedJob );
    void dispatch( const List<func_and_name> & pNamedJobs );

    void addLogEntry( const std::thread::id pId, const Jobs & pLog, const bool update ) ;


    /**
     * @brief dispatch: Dispatch and move
     * @param op
     */
    void dispatch( func && pJob, const string & pName="<unnamed>" );

	// Deleted operations
    DispatchQ( const DispatchQ & ) = delete;
    DispatchQ & operator=( const DispatchQ & ) = delete;
    DispatchQ( DispatchQ && ) = delete;
    DispatchQ & operator=( DispatchQ && ) = delete;


    /**
     * @brief ExportLogsJson
     * @param pToDir file name (without extention) to export the logs to
     * @details Export the logs to a json file. The file will be created if it does not exist.
     * The file will be overwritten if it already exists.
     * The logs are exported in the following format:
     * {
        "final_time": int,
        "nb_threads": int,
        "nb_jobs_OK": int,
        "nb_jobs_KO": int,
        "threads": [
            {
                "nomT": "Thread_0",
                "idT": "int64_t",
                "jobs": [
                    {
                        "nomJ": string,
                        "id": int64_t,
                        "id parent": int64_t,
                        "debut": int64_t,
                        "fin": int64_t,
                        "success": bool,
                        "error": string
                    },
                    ...
                ],
                "StateSwitch": [
                    {
                        "Name": string,
                        "Need": int64_t,
                        "Have": int64_t,
                        "Free": int64_t
                    },
                    ...
                ],
                "Tags": [
                    {
                        "Name": string,
                        "Start": int64_t,
                        "Stop": int64_t
                    },
                    ...
                ]
            },
        ...
        ]
    }
    */
    void ExportLogsJson(const std::string & pToDir);

protected:
    /**
     * @brief finish
     * @param pQuitReason
     * @param pRestartPool
     * @param pExportLogs If true, export logs to file if mLogsMod is not OFF (default: true) (false just in the destructor)
     */
    void finish( bool & pQuitReason, bool pRestartPool, bool pExportLogs = true );


    /**
     * @brief statusChangeMutex
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the mutex
     * @param pThreadMutexState state of the thread mutex (NEED, HAVE, FREE)
     * @details This function is used to change the state of the thread mutex.
     * It is used to log when a thread is waiting or handling a mutex.
     * @example std::mutex mtx; // Simple mutex
     * MutexNeed(pDQ, "mtx");
     * mtx.lock();     // Lock the mutex
     * MutexHave(pDQ, "mtx");
     * // ... Critical Section ...
     * MutexFree(pDQ, "mtx");
     * mtx.unlock();   // Déverrouille le mutex
     * @warning The user must ensure that pName is unique for each mutex, 
     * and that the exact same name is consistently used across 
     * all three statusChangeMutex instances associated with a single mutex
     */
    void statusChangeMutex(const string& pName, const ThreadMutexState& pThreadMutexState);


    /**
     * @brief ChangeTags
     * 
     * @param pDQ reference to the DispatchQ object
     * @param pName name of the tag
     * @param pTagState state of the tag (START, STOP)
     * @details This function is used to change the state of a part of code.
     * It is used to log when a thread is starting or stopping a special part of code.
     * @example TagStart(pDQ, "MyTag");
     * // ... Special part of code ...
     * TagStop(pDQ, "MyTag");
     * @warning The user must ensure that pName is unique for each tag,
     * and that the exact same name is consistently used across
     * all two ChangeTags instances associated with a single tag
     */
    void ChangeTags(const std::string& pName,  const TagState& pTagState);



private:
    // JSON management

    /**
     * @brief InitRealTimeJson
     * 
     */
    void InitRealTimeJson();

    /**
     * @brief UpdateHeaderJson
     * 
     */
    void UpdateHeaderJson();

    /**
     * @brief UpdateJsonJobStart
     * 
     * @param pJob 
     */
    void UpdateJsonJobStart(const Jobs& pJob);

    /**
     * @brief UpdateJsonJobStop
     * 
     * @param pJob 
     * @param pIsSuccess 
     */
    void UpdateJsonJobStop(const Jobs& pJob);

    /**
     * @brief UpdateQuittingJsonJob
     * 
     * @param pJob 
     */
    void UpdateQuittingJsonJob(const Jobs& pJob);

    /**
     * @brief UpdateJsonMutex
     * 
     * @param pThreadMutexState 
     * @param pName 
     * @param pTime 
     */
    void UpdateJsonMutex(const ThreadMutexState& pThreadMutexState,const string pName, const int64_t pTime) ;


    void UpdateJsonDeadLock();

    /**
     * @brief UpdateJsonTag
     * 
     * @param pTagState 
     * @param pName 
     * @param pTime 
     */
    void UpdateJsonTag(const TagState& pTagState,const string pName, const int64_t pTime);

    /**
     * @brief ExportRealTimeJson
     * 
     * @param pToDir is the path to the directory 
     */
    void ExportRealTimeJson(const string& pToDir);


private:
    // Main Thread Function
    
    /**
     * @brief threadFunc is run by the DQ threads.
     * It deals with the jobs: processing and exception management.
     */
    void threadFunc();
    

private:
        
    using time_point = std::chrono::time_point<std::chrono::steady_clock>;

    time_point sT0;                             ///< Time point for the start of the profiling
    
    string mName;                               ///< Name of the DispatchQ
    vector<std::thread> mThreads;               ///< Vector of threads
    mutex mTexMex;                              ///< Mutex for the queue    

    queue<Jobs> mQ;                             ///< Queue of jobs (FIFO) : protected by mTexMex
    std::atomic<int> mNbWorksInProgress{ 0 };   ///< Number of works in progress (atomic)
    condition_variable mCV;                     ///< Condition variable for the queue (used to wake up threads when a job is added to the queue)
    
    //Quit status
    bool mForceQuit{ false };                   ///< Force quit status (true if the user wants to quit the thread pool)
    bool mEmptyQuit{ false };                   ///< Empty quit status (true if the user wants to quit the thread pool when the queue is empty)
    bool mIdleQuit{ false };                    ///< Idle quit status (true if the user wants to quit the thread pool when all threads are idle)

    size_t mNb_Success {0};                     ///< Number of successful jobs
    size_t mNb_Errors {0};                      ///< Number of failed jobs
    size_t mNbQuitting {0};                     ///< Number of quitting threads (atomic)


    // Profiling or not?
    map<std::thread::id, DataThread> mLogs;     ///< Map of logs (key: thread id, value: DataThread object)
    mutex mLogMex;                              ///< Mutex for the logs (protected by mLogMex)            
    map<std::thread::id, Jobs*> mCurrentJobs;   ///< Map of current jobs (key: thread id, value: Jobs object) (For the heredity of Jobs)
    mutex mCurrentJobsMutex;                    ///< Mutex for the current jobs (protected by mCurrentJobsMutex)
    LogsMod mLogsMod {LogsMod::OFF};            ///< Logging mode (OFF, END, REALTIME)
    string mToDir;                             ///< File to export the logs to (if empty, no file is created)

    // Real Time Logs
    ordered_json mRealTimeJson;                 ///< Json to Real Time Logs
    std::thread mRealTimeThread;                ///< Thread for real time logs
    condition_variable mRealTimeCV;
    std::atomic<size_t> mNbUpdate{ 0 };         ///< Number of updates for real time logs
    bool mAllDqWorkerQuit {false};              ///< boolean if all DQ workers quit (inable by default)
    size_t mRunningJobs = 0;                    ///< Number of running jobs in real time logs


    GraphCycleThreads mGraphCycle;             ///< Graph cycle for the threads (used to detect deadlocks)
    bool mCycleDetectionEnabled{ false };      ///< Enable cycle detection (true if cycle detection is enabled)

};


//================================================================================
/**
 * @brief friend functions for DispatchQ class
 * 
 * @param pDQ 
 * @param pName 
 */

// Mutex logs
inline void MutexNeed(DispatchQ& pDQ, const string& pName) { pDQ.statusChangeMutex(pName, ThreadMutexState::NEED); }
inline void MutexHave(DispatchQ& pDQ, const string& pName) { pDQ.statusChangeMutex(pName, ThreadMutexState::HAVE); }
inline void MutexFree(DispatchQ& pDQ, const string& pName) { pDQ.statusChangeMutex(pName, ThreadMutexState::FREE); }

//
inline void TagStart(DispatchQ& pDQ, const std::string& pName) { pDQ.ChangeTags(pName, TagState::START); }
inline void TagStop(DispatchQ& pDQ, const std::string& pName) { pDQ.ChangeTags(pName, TagState::STOP); }



}
