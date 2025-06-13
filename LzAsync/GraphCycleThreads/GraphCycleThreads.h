#ifndef __GRAPH_CYCLE__
#define __GRAPH_CYCLE__

// ======================================================================================
/**
 * @file GraphCycle.h
 * @brief Class for unify all jobs in one structure
 * @details This class is used in Logs class to be the jobs of DispatchQ class. 
It contains all the informations needed to create logs.
 */

#include <chrono>
#include <string>
#include <unordered_set>
#include <thread>
#include <array>
#include <vector>
#include <map>
#include <stdexcept>
#include <mutex>
#include "../DataThread/DataThread.h"
#include "../LzLib_Async.h"
#include <iostream>
#include "../json.hpp"



using std::string;
using std::vector;
using std::map;
using std::mutex;
using nlohmann::json_abi_v3_12_0::ordered_json;
using time_point = std::chrono::time_point<std::chrono::steady_clock>;

// ======================================================================================
/**
 * @brief Thread dependency cycle detector
 * @details Tracks mutex dependencies between threads to detect potential deadlocks
 * by modeling thread-wait relationships as a directed graph and checking for cycles.
 * 
 * Maintains mappings between threads, owned mutexes, and requested mutexes to
 * construct a wait-for graph where cycles indicate deadlock conditions.
 */
class GraphCycleThreads
{
public:
    /**
    * @brief Default constructor
    **/
    GraphCycleThreads() {}  ;

    /**
    * @brief Initializes thread mapping structure
    * @param pThreads Reference to vector of threads to monitor
    * @param sT0 The begining time point
    * @details
    * Creates integer mappings for thread IDs to enable efficient graph operations
    * and establishes initial graph vertex count based on thread pool size.
    */
    void initialize(vector<std::thread>& pThreads, time_point sT0);

    // GraphCycleThreads(vector<std::thread>& pThreads);
    
    /**
     * @brief Adds directed edge to wait-for graph
     * @param pU Source vertex (thread integer ID)
     * @param pV Destination vertex (thread integer ID)
     * @details
     * Represents that thread U is waiting for a resource held by thread V.
     * Edge additions are thread-safe through internal mutex protection.
     */
    void AddEdge(int pU, int pV);

    /**
     * @brief Removes the outgoing edge from a vertex
     * @param pU Vertex (thread integer ID) to clear edges from
     * @details
     * Typically called when a thread releases all resources, eliminating
     * its dependencies in the wait-for graph.
     */
    void EraseEdge(int pU);

    /**
    * @brief Detects cycles in the wait-for graph
    * @details
    * Implements Floyd’s cycle detection algorithm (tortoise and hare) to find
    * potential deadlocks by detecting cycles in the dependency graph.
    * @return True if a cycle (deadlock) is detected, false otherwise.
    */
    bool CycleDetector();

    /**
    * @brief Updates mutex state and checks for cycles
    * @param pName Name of mutex being modified
    * @param pThreadMutexState New state of the mutex (NEED/HAVE/FREE)
    * @details
    * Atomic operation that:
    * 1. Takes mutex lock (mMexCycle) for thread-safe updates
    * 2. Updates both mutex tracking maps (mMapMutexThreadHave/mMapThreadMutexNeed)
    * 3. Modifies wait-for graph edges through called UpdateMutex* methods
    * 4. Maintains consistent state between maps and graph
    */
    bool UpdateCycle(const string& pName, const ThreadMutexState& pThreadMutexState);

    /**
    * @brief Records mutex acquisition and resolves dependencies
    * @param pIdThread ID of thread acquiring mutex
    * @param pMutexName Name of acquired mutex
    * @details
    * Performs three atomic operations under mutex protection:
    * 1. Updates mMapMutexThreadHave (mutex → owner mapping)
    * 2. Removes entry from mMapThreadMutexNeed if present
    * 3. Removes all incoming wait edges in mGraph for this thread-mutex pair
    * 
    * Eliminates the recorded dependency as the thread no longer waits
    * for this resource.
    */
    void UpdateMutexHave(std::thread::id pIdThread, string pMutexName) ;

    /**
    * @brief Records mutex request by thread and updates wait-for graph
    * @param pIdThread ID of thread requesting mutex
    * @param pMutexName Name of requested mutex
    * @details
    * Performs three atomic operations under mutex protection:
    * 1. Updates mMapThreadMutexNeed (thread → required mutex mapping)
    * 2. Checks current owner in mMapMutexThreadHave
    * 3. If mutex is owned by another thread (T_owner), adds edge:
    *    T_requester → T_owner to the wait-for graph (mGraph)
    * 
    * This edge creation represents a potential deadlock dependency.
    * Cycle detection should be triggered after this operation.
    */
    void UpdateMutexNeed(std::thread::id pIdThread, string pMutexName) ;

   /**
    * @brief Records mutex release and updates dependency graph
    * @param pIdThread ID of thread releasing mutex
    * @param pMutexName Name of released mutex
    * @details
    * Performs three atomic operations under mutex protection:
    * 1. Removes ownership record from mMapMutexThreadHave
    * 2. Checks if any thread was waiting for this mutex (via mMapThreadMutexNeed)
    * 3. For each waiting thread, removes corresponding edge from mGraph
    * 
    * Note: Does not automatically grant the mutex to waiting threads,
    * but removes the deadlock-inducing wait edges from the graph.
    */
    void UpdateMutexFree(std::thread::id pIdThread, string pMutexName) ;

    /**
    * @brief Stores detected cycle path
    * @param pStart Starting vertex of the cycle
    * @details
    * Called internally when CycleDetector identifies a cycle.
    * Populates mCyclePath with vertex sequence representing the deadlock loop.
    */
    void SaveCycle(int pStart);


/**
* @brief Serializes detected deadlock information into JSON format
* @param pJson Reference to ordered_json object to receive deadlock data
* @details
* Constructs a detailed JSON structure containing:
* - Timestamp of deadlock detection (milliseconds since program start)
* - Complete cycle information including:
*   * Thread IDs involved in deadlock
*   * Mutexes each thread is waiting for (Need)
*   * Mutexes each thread currently holds (Have)
* 
* JSON Structure Example:
* @code{.json}
*   "DeadLock": {
*     "Time": 123456,
*      "Cycle": [
*       {
*         "ThreadId": "1",
*         "Need": "mutexA",
*         "Have": ["mutexB"]
*       },
*       {
*         "ThreadId": "2",
*         "Need": "mutexB",
*         "Have": ["mutexA"]
*       }
*     ]
*   }
* }
* @endcode
* 
* @note
* - Uses ordered_json to maintain element ordering
* - Timestamp is relative to program start time (T0)
*/
    void UpdateJsonDeadLock(ordered_json &pJson);

private:
    map<string,int>  mMapMutexThreadHave ;              //Map to record mutex ownership by thread
    map<int,string>  mMapThreadMutexNeed ;              //Map to record mutex requests by thread    

    map<std::thread::id,int> mMapThreadInt ;            //Map to convert thread id into int fast
    map<int,std::thread::id> mMapIntThread ;            //Map to convert int into thread id fast
    
    int mNbVerticies;                                   //Number of vertices in the graph (= nb threads) 
    std::vector<int> mGraph;                            //Adjacency list representation of the wait-for graph   

    std::vector<std::thread::id> mCyclePath ;           //Path of vertices forming the detected cycle
    mutex mMexCycle;

    using time_point = std::chrono::time_point<std::chrono::steady_clock>;
    time_point T0;                                      ///< Time point for the start of the profiling
    //ordered_json mRealTimeJson;                          //Json to store the cycle information if a cycle is detected
};
// ======================================================================================
#endif 