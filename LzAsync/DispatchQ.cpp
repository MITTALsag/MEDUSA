#include "DispatchQ.h"
#include <LzServices/LzLog.h>
#include <regex>

#include "json.hpp"



namespace LzAsync
{

//================================================================================
DispatchQ::DispatchQ( const string & pName/*="<unnamed>"*/, size_t pNbWorkers/*=1*/ )
 : mName(pName), 
 mThreads(pNbWorkers)
{
    // Check
    if( pNbWorkers == 0 )
        LzLogException("", "Cannot create a dispatch queue with 0 workers!")

	// Log
    LzLogM("", "DispatchQ::DispatchQ: Creating dispatch queue: '"<<pName<<"'. Thread count= "<<pNbWorkers<<".")

    // Create thread pool
    for( std::thread & iTh : mThreads )
        iTh = std::thread( &DispatchQ::threadFunc, this );

    // Log
    LzLogM("", "DispatchQ::DispatchQ: Pool '"<<mName<<"' created with "<<mThreads.size()<<" worker(s).")

}

//================================================================================
DispatchQ::~DispatchQ()
{
    // Log
    LzLogN("", "DispatchQ::~DispatchQ: '"<<mName<<"'...")

    // Blocking: terminate jobs and do NOT restart pool
    interruptJobs(false, false);
    
    // FINAL STATE: jobs partially completed, no more jobs in queue, thread pool NOT running
}

//================================================================================
void DispatchQ::startProfiling(const string& pToDir, const LogsMod pLogsMod /*=LogsMod::END*/)
{
    // this is thread safe, since it is called only from the main thread

    // Reset profiling data
    sT0 = std::chrono::steady_clock::now();
    mNb_Errors = 0;
    mNb_Success = 0;
    mNbQuitting = 0;
    mRunningJobs = 0;
    mCurrentJobs.clear();
    mLogsMod = pLogsMod;
    mToDir = pToDir;
    Jobs::ResetJobCounter();

    switch(pLogsMod)
    {
        case LzAsync::LogsMod::END:
        {
            LzLogM("", "DispatchQ::startProfiling: Profiling enabled in END mode, logs will be exported to '"<<pToDir<<"'.")
            mLogs.clear();
            break;
        }
        case LzAsync::LogsMod::REALTIME:
            LzLogM("", "DispatchQ::startProfiling: Profiling enabled in REALTIME mode, logs will be exported to '"<<pToDir<<"'.")
            mRealTimeThread = std::thread(&DispatchQ::ExportRealTimeJson, this, pToDir);
            mRealTimeJson.clear(); // Clear previous real-time logs
            mNbUpdate = 0; // Reset update counter for real-time logs
            mAllDqWorkerQuit = false; // Reset quitting status for real-time logs
            InitRealTimeJson();
            mCycleDetectionEnabled = true; // Enable Cycle Detection
            // Initialize GraphCycle
            mGraphCycle.initialize(mThreads, sT0);
            break;
        case LzAsync::LogsMod::OFF:
            LzLogM("", "DispatchQ::startProfiling: Profiling disabled.")
            break;
    }
}


//================================================================================
void DispatchQ::finish( bool & pQuitReason, bool pRestartPool, bool pExportLogs /*=true*/ )
{
    // Lock mutex
    std::unique_lock<std::mutex> lLock(mTexMex);

    // Set predicate value
    pQuitReason = true;

    // Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again
    lLock.unlock();
    mCV.notify_all();

    // Wait for threads to finish before we exit
    for( size_t i=0 ; i<mThreads.size() ; i++ )
    {
        // Still working?
        if( mThreads[i].joinable() )
        {
            // Log and join
            LzLogM("", "DispatchQ::finish: Waiting for thread "<<i<<" (id= "<<mThreads[i].get_id()<<"), in '"<<mName<<"' to join...")
            mThreads[i].join();
            mNbQuitting++;
        }
        else
        {
            // Log
            LzLogM("", "DispatchQ::finish: Thread "<<i<<" (id= "<<mThreads[i].get_id()<<"), in '"<<mName<<"' is not joinable.")
        }
    }

    // Log
    LzLogM("", "DispatchQ::finish: Joined all "<<mThreads.size()<<" thread(s) in '"<<mName<<"'.")

    // HERE: Accessing data members without locking the mex.
    // -----
    //       This is not a problem since all workers have stopped working AND we assume
    //       that there is only one thread using this DispatchQ object. See note about
    //       thread-*UN*safety of this class in the header.

    // Empty queue
    mQ = std::queue<Jobs>();
    //
    // use std::deque(backend of queue) and q_.clear() ?

    // Set number of works in progress to 0 (this should already be the case)
    mNbWorksInProgress = 0;

    // It ain't over til it's over! The DQ can be reused
    pQuitReason = false;


    // Restart pool thread pool?
    if( pRestartPool )
    {
        // Log
        LzLogM("", "DispatchQ::finish: Restarting '"<<mName<<"' with "<<mThreads.size()<<" thread(s).")

        // Restart
        for( size_t i=0 ; i<mThreads.size() ; i++ )
            mThreads[i] = std::thread( &DispatchQ::threadFunc, this );
    }
    else
    if (pExportLogs)
    {
        LzLogM("", "DispatchQ::finish: Not restarting '"<<mName<<"'.")
        // Log

        if (mLogsMod == LzAsync::LogsMod::END)
        {
            if (mToDir.empty())
            {
                throw std::runtime_error("DispatchQ::finish: No file specified for logs export in END mode.");
            }
            // Export logs to file
            ExportLogsJson(mToDir);
        }
        else if (mLogsMod == LzAsync::LogsMod::REALTIME)
        {

            mAllDqWorkerQuit = true;
            // Notify mRealTimeThread to stop
            mRealTimeCV.notify_one();
            // Wait for mRealTimeThread to finish
            if (mRealTimeThread.joinable())
            {
                mRealTimeThread.join();
                LzLogM("", "DispatchQ::finish: Real time thread joined.")
            }
        }
    }
}

//================================================================================

void DispatchQ::InitRealTimeJson()
{
    mRealTimeJson["final_time"] = nullptr;
    mRealTimeJson["nb_threads"] = nullptr;
    mRealTimeJson["nb_jobs_OK"] = nullptr;
    mRealTimeJson["nb_jobs_KO"] = nullptr;
    mRealTimeJson["nb_running_jobs"] = nullptr;
    mRealTimeJson["nb_quitting"] = nullptr;

    vector<ordered_json> lThreadVector;
    ordered_json lTh;
    for(int iTh = 0; iTh < mThreads.size(); iTh++)
    {
        lTh["nomT"] = "Thread_" + std::to_string(iTh);
        
        std::ostringstream lSS;
        // Conversion of this id into a string
        lSS << mThreads[iTh].get_id(); 
        lTh["idT"] = lSS.str();

        vector<ordered_json> lEmptyVector;
        lTh["jobs"] = lEmptyVector;
        lTh["StateSwitch"] = lEmptyVector;
        lTh["Tags"] = lEmptyVector;

        lThreadVector.push_back(lTh);
    }
    mRealTimeJson["threads"] = lThreadVector;
    mRealTimeJson["DeadLock"]= nullptr;
}

//================================================================================

void DispatchQ::ExportLogsJson (const std::string & pToDir)
{
    LzLogM("", "DispatchQ::ExportLogsJson exporting in file : "<<pToDir+".json")
    std::ofstream lF(pToDir+".json"); /*Opening the file pToDir.json*/
    if(!lF)
    {
        LzLogM("", "DispatchQ::ExportLogsJson cannot open :"<<pToDir+".json")
        return ;
    }
    ordered_json lJson ;
    time_point lNow = LzServices::Chrono_Now();  
    int64_t lMsecNow = std::chrono::duration_cast<std::chrono::milliseconds>(lNow - sT0).count();

    lJson["final_time"] = lMsecNow;
    lJson["nb_threads"] = mLogs.size();
    lJson["nb_jobs_OK"] = mNb_Success;
    lJson["nb_jobs_KO"] = mNb_Errors;
    lJson["nb_running_jobs"] = 0;
    lJson["nb_quitting"] = mNbQuitting;
    vector<ordered_json> lThreadsVector ; //  JSON vector for the threads
    int iThread = 0 ;

    for (const auto & iPair : mLogs)
    {
        ordered_json lTreadJson;

        std::thread::id lthreadId = iPair.first;
        DataThread ldataThread = iPair.second;

        lTreadJson["nomT"] = "Thread_"+std::to_string(iThread++);

        std::ostringstream lSS;
        // Conversion of this id into a string
        lSS << lthreadId; 
        lTreadJson["idT"] = lSS.str();

        vector<ordered_json> lJsonJobsVector ;

        // typeof(iPair) = pair<tread_id, DataThread>
        const vector<Jobs>& ljobs = ldataThread.GetJobs();

        for (const Jobs& iLog : ljobs)  
        {
            ordered_json lJsonJob ;
            lJsonJob["nomJ"] = iLog.GetJobName() ;
            lJsonJob["id"] = iLog.GetJobId();
            lJsonJob["id parent"] = iLog.GetJobIdParent();
            lJsonJob["debut"] = iLog.GetStartmsec();
            lJsonJob["fin"] = iLog.GetStopmsec() ;
            lJsonJob["success"] = iLog.IsSuccess();
            lJsonJob["error"] = iLog.GetError() ;
            lJsonJobsVector.push_back(lJsonJob);
        }   

        vector<ordered_json> lJsonMutexVector;

        const std::map<string, vector<MutexContainer::MutexLogs>>& lMutexLogs = ldataThread.getMutexLogs();

        for (const auto& [iMutexName, iMutexData] : lMutexLogs) 
        {
            for (const auto& iMutexLog : iMutexData) 
            {
                ordered_json lJsonMutex;
                lJsonMutex["NameM"] = iMutexName;
                lJsonMutex["NEED"] = iMutexLog[0];
                lJsonMutex["HAVE"] = iMutexLog[1];
                lJsonMutex["FREE"] = iMutexLog[2];
                lJsonMutexVector.push_back(lJsonMutex);
                
            }
        }

        // Tags part
        ordered_json lTags;
        vector<ordered_json> lJsonTagVector ;
        const vector<TagContainer::Tag>& lStatusTag = ldataThread.getTags();

        for (const auto& iTagData : lStatusTag)
        {
            ordered_json lJsonTag ;

            lJsonTag["NameT"] = get<0>(iTagData); ;
            lJsonTag["START"] = get<1>(iTagData)[0];
            lJsonTag["STOP"] = get<1>(iTagData)[1];

            lJsonTagVector.push_back(lJsonTag);
        }

        lTreadJson["jobs"] = lJsonJobsVector ;
        lTreadJson["StateSwitch"] = lJsonMutexVector; 
        lTreadJson["Tags"] = lJsonTagVector; 
        lThreadsVector.push_back(lTreadJson); 
    }
    lJson["threads"] = lThreadsVector ;
    lJson["DeadLock"]= nullptr;
    lF<<lJson.dump(4)<<endl ;
    lF.close() ;
    LzLogM("", "DispatchQ::ExportLogsJson finish")
    return ;
}

//================================================================================

void DeleteFile(const string& pDir, const string& pPattern) {

    // Regular expression to match files like "pPattern"
    std::regex pattern(pPattern);

    try {
        // Iterate over all files in the directory
        for (const auto& entry : std::filesystem::directory_iterator(pDir)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                
                // Check if the filename matches the pattern
                if (std::regex_match(filename, pattern)) {
                    // Delete the file
                    std::filesystem::remove(entry.path());
                    // LzLogM("", "delete " << entry.path());
                }
            }
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

//================================================================================  

void DispatchQ::ExportRealTimeJson (const std::string & pToDir)
{
    LzLogN("", "DispatchQ::ExportRealTimeJson: Exporting in directory : "<< pToDir)

    // To save file names
    int lcounter = 0;
    constexpr int NB_FILE = 10; // Number of file to save
    std::array<string, NB_FILE> lNameFileSave;
    bool first = true;

    // Verify if the directory exit
    if (!std::filesystem::exists(pToDir)) 
    {
        // Create it
        if (std::filesystem::create_directory(pToDir)) 
        {
            LzLogM("","Directory creation suceed" << pToDir);
        } 
        else 
        {
            throw std::runtime_error("Directory creation error " + pToDir);
        }
    } 
    else 
    {
        LzLogM("","Directory already exist." << pToDir);
        DeleteFile(pToDir, "_LOG_\\d+\\.json");
    }

    while (true)
    {
        // Lock the mutex to access the JSON data
        std::unique_lock<std::mutex> lLock(mLogMex);

        // Wait for a notification
        bool lTimeOut = mRealTimeCV.wait_for(lLock, std::chrono::milliseconds(10000), [&] { return mNbUpdate > 0 || mAllDqWorkerQuit; });
                
        if (mAllDqWorkerQuit && mNbUpdate == 0)
        {
            // LzLogM("", "DispatchQ::ExportRealTimeJson: All DQ workers quit. Exiting thread.");
            break;
        }
        else if (!lTimeOut)
        {
            // LzLogM("", "DispatchQ::ExportRealTimeJson: Wake up with time condition.")
            UpdateHeaderJson();
        }
        else
        {
            // LzLogM("", "DispatchQ::ExportRealTimeJson: Wake up with update condition.")
        }

        // Cpoy the JSON data to a local variable
        ordered_json localJson = mRealTimeJson;
        mNbUpdate = 0;

        // Unlock the mutex before writing to the file
        lLock.unlock(); 


        // Write the JSON data to the file
        
        int64_t lTimeStamp = localJson["final_time"];
        
        string lFileName = pToDir + "/_LOG_" + std::to_string(lTimeStamp) + ".json";
        string lFileBefore;
        if (lcounter < NB_FILE)
        {
            if (!first)
            {
                lFileBefore = lNameFileSave[lcounter];
            }
            lNameFileSave[lcounter] = lFileName;
            lcounter++;
            if (first)
                first = false;
        }
        else
        {
            lFileBefore = lNameFileSave[0];
            lNameFileSave[0] = lFileName;
            lcounter = 1;
        }
        DeleteFile(pToDir, "R(_LOG_\\d+\\.json)");

        //Delete before file : 
        // Check if the file exists
        if (std::filesystem::exists(lFileBefore)) 
        {
            // Attempt to delete the file
            if (std::filesystem::remove(lFileBefore)) 
            {
                // LzLogM("", "File successfully deleted: " << lFileBefore);
            } 
            else 
            {
                throw std::runtime_error("Error deleting file: " + lFileBefore);
            }
        }


        std::ofstream lF(lFileName);
        if (!lF)
        {
            throw std::runtime_error("DispatchQ::ExportRealTimeJson : Cannot open file " + lFileName);
        }

        lF << localJson.dump(4) << std::endl;
        lF.close();

    }

}

//================================================================================

void DispatchQ::UpdateHeaderJson()
{
    // Update the header of the JSON file

    // Calcule final time (relative to execution)
    time_point lNow = LzServices::Chrono_Now();  
    int64_t lMsecNow = std::chrono::duration_cast<std::chrono::milliseconds>(lNow - sT0).count();
    mRealTimeJson["final_time"] = lMsecNow;
    mRealTimeJson["nb_threads"] = mThreads.size();
    mRealTimeJson["nb_jobs_OK"] = mNb_Success ;
    mRealTimeJson["nb_jobs_KO"] = mNb_Errors;
    mRealTimeJson["nb_running_jobs"] = mRunningJobs;
    mRealTimeJson["nb_quitting"] = mNbQuitting;

}

//================================================================================

void DispatchQ::UpdateJsonJobStart(const Jobs& pJob)
{
    std::unique_lock<std::mutex> llock(mLogMex);

    // Increment the number of running jobs
    mRunningJobs++;
    
    // Update Header
    UpdateHeaderJson();

    // Find the thread corresponding to the job
    std::thread::id lThreadId = std::this_thread::get_id();
    std::ostringstream lSS;
    lSS << lThreadId; 
    string lIdThread = lSS.str();

    // Iterate over the threads to find the one that matches
    for (auto& thread : mRealTimeJson["threads"]) 
    {
        if (thread["idT"] == lIdThread) 
        {
            ordered_json lJsonJob;
            lJsonJob["nomJ"] = pJob.GetJobName();
            lJsonJob["id"] = pJob.GetJobId();
            lJsonJob["id parent"] = pJob.GetJobIdParent();
            lJsonJob["debut"] = pJob.GetStartmsec();
            lJsonJob["fin"] = -1;
            lJsonJob["success"] = -1;
            lJsonJob["error"] = -1;
            thread["jobs"].push_back(lJsonJob);
            break;
        }
    }

    mNbUpdate++;
    llock.unlock();
    mRealTimeCV.notify_one();

}

//================================================================================


void DispatchQ::UpdateJsonJobStop(const Jobs& pJob) 
{
    std::unique_lock<std::mutex> llock(mLogMex);
    // Decrement the number of running jobs
    mRunningJobs--;

    // Update the number of successful or failed jobs
    if (pJob.IsSuccess()) 
        mNb_Success++ ;
    else 
        mNb_Errors++ ;
    
    // Update Header
    UpdateHeaderJson();

    // Find the thread corresponding to the job
    std::thread::id lThreadId = std::this_thread::get_id();
    std::ostringstream lSS;
    lSS << lThreadId; 
    string lIdThread = lSS.str();

    // Iterate over the threads to find the one that matches
    for (auto& thread : mRealTimeJson["threads"])
    {
        if (thread["idT"] == lIdThread) 
        {
            ordered_json& lJsonJob = thread["jobs"].back();
            lJsonJob["fin"] = pJob.GetStopmsec();
            lJsonJob["success"] = pJob.IsSuccess();
            lJsonJob["error"] = pJob.GetError();
            break;
        }
    }

    mNbUpdate++;
    llock.unlock();
    mRealTimeCV.notify_one();
}

//================================================================================
void DispatchQ::UpdateQuittingJsonJob(const Jobs& pJob) 
{
    std::unique_lock<std::mutex> llock(mLogMex);
    // Decrement the number of running jobs
    mNbQuitting++;

    // Update Header
    UpdateHeaderJson();

    // Find the thread corresponding to the job
    std::thread::id lThreadId = std::this_thread::get_id();
    std::ostringstream lSS;
    lSS << lThreadId; 
    string lIdThread = lSS.str();

    // Iterate over the threads to find the one that matches
    for (auto& thread : mRealTimeJson["threads"])
    {
        if (thread["idT"] == lIdThread) 
        {
            ordered_json lJsonJob;
            lJsonJob["nomJ"] = pJob.GetJobName();
            lJsonJob["id"] = pJob.GetJobId();                   //-1
            lJsonJob["id parent"] = pJob.GetJobIdParent();      //-1    
            lJsonJob["debut"] = pJob.GetStartmsec();
            lJsonJob["fin"] = pJob.GetStopmsec(); 
            lJsonJob["success"] = pJob.IsSuccess();
            lJsonJob["error"] = pJob.GetError();
            thread["jobs"].push_back(lJsonJob);
            break;
        }
    }

    mNbUpdate++;
    llock.unlock();
    mRealTimeCV.notify_one();
}

//================================================================================

void DispatchQ::UpdateJsonMutex(const ThreadMutexState& pThreadMutexState,const string pName, const int64_t pTime)
{
    std::unique_lock<std::mutex> llock(mLogMex);
    mNbUpdate++;

    // Update Header
    UpdateHeaderJson();

    // Find the thread corresponding to the job
    std::thread::id lThreadId = std::this_thread::get_id();
    std::ostringstream lSS;
    lSS << lThreadId; 
    string lIdThread = lSS.str();

    // Iterate over the threads to find the one that matches
    for (auto& iThread : mRealTimeJson["threads"])
    {
        if (iThread["idT"] == lIdThread) 
        {
            // Update the Mutex state
            switch(pThreadMutexState)
            {
                case ThreadMutexState::NEED:
                {
                    // If the Mutex is in NEED state, we create a new Mutex
                    ordered_json lJsonMutex;
                    lJsonMutex["NameM"] = pName;
                    lJsonMutex["NEED"] = pTime ;
                    lJsonMutex["HAVE"] = -1 ;
                    lJsonMutex["FREE"] = -1 ;
                    iThread["StateSwitch"].push_back(lJsonMutex) ;
                    break;
                }
                case ThreadMutexState::HAVE:
                case ThreadMutexState::FREE:
                {
                    // We do a reverse loop to find the last Mutex with the same name
                    // and update its state
                    int iIndexVector =  iThread["StateSwitch"].size() - 1 ;
                    ordered_json lJsonLastMutex = iThread["StateSwitch"][iIndexVector] ;

                    while (lJsonLastMutex["NameM"] != pName && iIndexVector > 0) 
                        lJsonLastMutex = iThread["StateSwitch"][--iIndexVector];
                    
                    // If we reach the beginning of the vector and the name is not found
                    // we throw an error
                    if (lJsonLastMutex["NameM"] != pName)
                        throw std::runtime_error("DispatchQ::UpdateJson : Mutex "+pName+" not found in JSON");
                    else
                    if (pThreadMutexState == ThreadMutexState::HAVE)
                        iThread["StateSwitch"][iIndexVector]["HAVE"] = pTime ;
                    else
                    if (pThreadMutexState == ThreadMutexState::FREE)
                        iThread["StateSwitch"][iIndexVector]["FREE"] = pTime ;
                }
            }
        }
    }
    llock.unlock();

    mRealTimeCV.notify_all();
}

//================================================================================

void DispatchQ::UpdateJsonTag(const TagState& pTagState,const string pName, const int64_t pTime)
{
    std::unique_lock<std::mutex> llock(mLogMex);
    mNbUpdate++;

    // Update Header
    UpdateHeaderJson();


    // Find the thread corresponding to the job
    std::thread::id lThreadId = std::this_thread::get_id();
    std::ostringstream lSS;
    lSS << lThreadId; 
    string lIdThread = lSS.str();


    // Iterate over the threads to find the one that matches
    for (auto& iThread : mRealTimeJson["threads"])
    {
        if (iThread["idT"] == lIdThread) 
        {
            switch(pTagState)
            {
                // If the tag is in START state, we create a new tag
                case TagState::START:
                {
                    // If the Tag is in START state, we create a new Tag
                    ordered_json lJsonTag;
                    lJsonTag["NameT"] = pName;
                    lJsonTag["START"] = pTime ;
                    lJsonTag["STOP"] = -1 ;
                    iThread["Tags"].push_back(lJsonTag) ;
                    break;
                }
                case TagState::STOP:
                {
                    // We do a reverse loop to find the last Tag with the same name
                    // and update its state
                    int iIndexVector =  iThread["Tags"].size() - 1 ;
                    ordered_json lJsonLastTag = iThread["Tags"][iIndexVector] ;

                    while (lJsonLastTag["NameT"] != pName && iIndexVector > 0) 
                        lJsonLastTag = iThread["Tags"][--iIndexVector];
                    
                    // If we reach the beginning of the vector and the name is not found
                    if (lJsonLastTag["NameT"] != pName)
                        throw std::runtime_error("DispatchQ::UpdateJson : Tags "+pName+" not found in JSON");
                    
                    iThread["Tags"][iIndexVector]["STOP"] = pTime;
                }
            }
        }
    }

    llock.unlock();
    mRealTimeCV.notify_all();
}

//================================================================================
void DispatchQ::restartPool( size_t pNbWorkers)
{
    // Check
    if( pNbWorkers == 0 )
        LzLogException("", "DispatchQ::restartPool("<<pNbWorkers<<"): Cannot restart dispatch queue with 0 workers!")

    // Log
    LzLogM("", "DispatchQ::restartPool("<<pNbWorkers<<"): Restarting '"<<mName<<"'...")

    // Blocking: terminate jobs and do NOT restart pool
    interruptJobs( false ); // Will export logs if startProfiling was called before.
    // FINAL STATE: jobs partially completed, no more jobs in queue, thread pool NOT running

    // HERE: Accessing data members without locking the mex.
    // -----
    //       This is not a problem since all workers have stopped working AND we assume
    //       that there is only one thread using this DispatchQ object. See note about
    //       thread-*UN*safety of this class in the header.

    
    // Stop profiling (the cleanup is done when we restart profiling)
    stopProfiling();

	// ResizeiThread
    mThreads.clear();
    mThreads.resize( pNbWorkers );
    for( std::thread & iTh : mThreads )
        iTh = std::thread( &DispatchQ::threadFunc, this );


    // Log
    LzLogM("", "DispatchQ::restartPool("<<pNbWorkers<<"): Pool '"<<mName<<"' restarted with "<<mThreads.size()<<" worker(s).")
}


//================================================================================
vector<std::thread::id> DispatchQ::GetThreadIds() const
{
    vector<std::thread::id> lThreadIds;
    lThreadIds.reserve(mThreads.size());
    for (const std::thread& iTh : mThreads)
    {
        lThreadIds.push_back(iTh.get_id());
    }
    return lThreadIds;
}

//================================================================================
void DispatchQ::dispatch(const std::function<void(void)> & pFunc, const string & pName/*="<unnamed>"*/ )
{
    int64_t lParentId = -1;
    {
        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
        auto it = mCurrentJobs.find(std::this_thread::get_id());
        if (it != mCurrentJobs.end()) 
        {
            lParentId = it->second->GetJobId();
        }
    }

    Jobs lJob(pName, pFunc, lParentId);    
    
    // Lock mutex
    std::unique_lock<std::mutex> lock(mTexMex);
    // Push single job
    mQ.push(lJob);

	// Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
	lock.unlock();
    mCV.notify_all();
}

//================================================================================
void DispatchQ::dispatch( const vector<func> & pFuncs, const string & pName/*="<unnamed>"*/ )
{
    int64_t lParentId = -1;
    {
        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
        auto it = mCurrentJobs.find(std::this_thread::get_id());
        if (it != mCurrentJobs.end()) 
        {
            lParentId = it->second->GetJobId();
        }
    }

    // Lock mutex
    std::unique_lock<std::mutex> lock(mTexMex);

    // Push vector of jobs
    for( size_t i=0 ; i<pFuncs.size() ; i++ )
    {
        Jobs itJobs(pName+" "+std::to_string(i), pFuncs[i], lParentId);
        mQ.push(itJobs);
    }

    // Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
    lock.unlock();
    mCV.notify_all();
}

//================================================================================
void DispatchQ::dispatch( const List<func> & pFuncs, const string & pName/*="<unnamed>"*/ )
{
    int64_t lParentId = -1;
    {
        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
        auto it = mCurrentJobs.find(std::this_thread::get_id());
        if (it != mCurrentJobs.end()) 
        {
            lParentId = it->second->GetJobId();
        }
    }
    
    // Lock mutex
    std::unique_lock<std::mutex> lock(mTexMex);

    // Push vector of jobs
    size_t i = 0;
    BrowseList( iJ, pFuncs )
    {
        Jobs itJobs(pName+" "+std::to_string(i), pFuncs.GetAt(iJ), lParentId);
        mQ.push(itJobs);
        i++;
    }

    // Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
    lock.unlock();
    mCV.notify_all();
}

//================================================================================
void DispatchQ::dispatch( func_and_name && pFN )
{
    // LzLogM("","dispatch(4)")

    int64_t lParentId = -1;
    {
        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
        auto it = mCurrentJobs.find(std::this_thread::get_id());
        if (it != mCurrentJobs.end()) 
        {
            lParentId = it->second->GetJobId();
        }
    }

    Jobs lJob(std::move(pFN).first, std::move(pFN).second, lParentId);

    // Lock mutex
    std::unique_lock<std::mutex> lock(mTexMex);

    // Push and move single job
    mQ.push( lJob ); // à voir si besoin de std::move

    // Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
    lock.unlock();
    mCV.notify_all();
}

//================================================================================
void DispatchQ::dispatch( const List<func_and_name> & pFNs )
{
    int64_t lParentId = -1;
    {
        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
        auto it = mCurrentJobs.find(std::this_thread::get_id());
        if (it != mCurrentJobs.end()) 
        {
            lParentId = it->second->GetJobId();
        }
    }


    // Lock mutex
    std::unique_lock<std::mutex> lock(mTexMex);

    // Push vector of jobs
    size_t i = 0;
    BrowseList( iJ, pFNs )
    {
        Jobs itJobs(pFNs.GetAt(iJ).first, pFNs.GetAt(iJ).second, lParentId) ;
        
        mQ.push( itJobs ); 
        i++;
    }

    // Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
    lock.unlock();
    mCV.notify_all();
}

//================================================================================
void DispatchQ::dispatch( func && pFunc, const string & pName/*="<unnamed>"*/ )
{
    // LzLogM("","dispatch(6)")

    int64_t lParentId = -1;
    {
        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
        auto it = mCurrentJobs.find(std::this_thread::get_id());
        if (it != mCurrentJobs.end()) 
        {
            lParentId = it->second->GetJobId();
        }

    }

    //LzLogM("","Parent id:"<<lParentId)
    Jobs lJob(pName, std::move(pFunc), lParentId);

    // Lock mutex
    std::unique_lock<std::mutex> lock(mTexMex);

    // Push and move single job
    mQ.push(lJob);//mQ.push( { pName, std::move(pFunc) } );

	// Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
	lock.unlock();
    mCV.notify_all();
}

//================================================================================

void DispatchQ::addLogEntry( const std::thread::id pId, const Jobs & pJobs, const bool pUpdate)
{
    // Lock mutex
    std::unique_lock<std::mutex> lLock(mLogMex);

    if (pJobs.IsSuccess() && pUpdate) 
        mNb_Success++;

    else if (!pJobs.IsSuccess() && pUpdate) 
        mNb_Errors++;

    mLogs[pId].AddLog(pJobs);

}

//================================================================================


void DispatchQ::statusChangeMutex(const string& pName, const ThreadMutexState& pThreadMutexState)
{
    try
    {
        time_point lnow = LzServices::Chrono_Now();
        int64_t lMsecNow = std::chrono::duration_cast<std::chrono::milliseconds>(lnow - sT0).count();

        if (mLogsMod==LzAsync::LogsMod::END)
            mLogs[std::this_thread::get_id()].setMutexTime(lMsecNow, pName, pThreadMutexState);
        else
        if (mLogsMod==LzAsync::LogsMod::REALTIME)
        {
           UpdateJsonMutex(pThreadMutexState,pName,lMsecNow);
            if(mGraphCycle.UpdateCycle(pName,pThreadMutexState))
             {
                //std::unique_lock<std::mutex> llock(mLogMex);

                mGraphCycle.UpdateJsonDeadLock(mRealTimeJson);

                std::this_thread::sleep_for( std::chrono::milliseconds(100) );
                //Sleep before updating the header, for visualisation purpose.

                UpdateHeaderJson();

                mNbUpdate++;

                mRealTimeCV.notify_one();

                LzLogM("", "/!\\ Warning: DeadLock detected ! Forced Interruption, logs exported in" + mToDir)
                
                // We don't want to execute exit 1 before the update 
                std::this_thread::sleep_for( std::chrono::milliseconds(50) );
                exit(1);

            }    
            
        }

    }
    catch (const std::exception& e)
    {
        LzLogException("", "Error in statusChangeMutex :\nException : " << e.what())
        throw; // Rethrow the exception to propagate it
    }
    catch (...)
    {
        LzLogException("", "Unknown error in statusChangeMutex")
        throw; // Rethrow the exception to propagate it
    }

}

//================================================================================

void DispatchQ::ChangeTags(const string& pName,  const TagState& pTagState) 
{
    try
    {
        time_point lnow = LzServices::Chrono_Now();
        int64_t lMsecNow = std::chrono::duration_cast<std::chrono::milliseconds>(lnow - sT0).count();

        if (mLogsMod==LzAsync::LogsMod::END)
            mLogs[std::this_thread::get_id()].setTagsTime(lMsecNow, pName, pTagState);
        else
        if (mLogsMod==LzAsync::LogsMod::REALTIME)
           UpdateJsonTag(pTagState,pName,lMsecNow);

    }
    catch (const std::exception& e)
    {
        LzLogException("", "Error in ChangeTags: " << e.what())
        throw; // Rethrow the exception to propagate it
    }
    catch (...)
    {
        LzLogException("", "Unknown error in statusChangeMutex")
        throw; // Rethrow the exception to propagate it
    }

}



//================================================================================
void DispatchQ::threadFunc()
{

    // Get thread id string
    const auto lTID = std::this_thread::get_id();

    // Log
    LzLogM("DQ worker", "Starting")

    // Stick around while I'm still needed
    while( true )
    {   
        // Lock mutex
        std::unique_lock<std::mutex> lLock(mTexMex);

        // Wait until we have data or a quit signal
        mCV.wait(lLock, [this]{ return mQ.size() != 0
                             || mForceQuit
                             || mEmptyQuit /* && mQ.size()==0 (can be logically omitted) */
                             || (mIdleQuit && mNbWorksInProgress==0) /* && mQ.size()==0 (can be logically omitted) */; });
        
        // Event 1: Force quit?
        if( mForceQuit )
        {
            // Log
            LzLogM("DQ worker", "Force quitting")

            // Profiling?
            if(mLogsMod==LzAsync::LogsMod::END)
            {
                time_point lTimePoint = LzServices::Chrono_Now();
                Jobs lJob(lTimePoint, sT0, "Quitting", "Empty quitting");
                addLogEntry(std::this_thread::get_id(), lJob, false);
            }
            else 
            if (mLogsMod==LzAsync::LogsMod::REALTIME)
            {
                time_point lTimePoint = LzServices::Chrono_Now();
                Jobs lJob(lTimePoint, sT0, "Quitting", "Empty quitting");
                UpdateQuittingJsonJob(lJob);
            }

            // Exit thread loop
            break;
        }

        // Event 2: Empty quit?
        //
        // DQ must be empty with no worker currently working on the chain
        //
        if( mEmptyQuit && mQ.size()==0 )
        {
            // Log
            LzLogM("DQ worker", "Empty quitting")



            // Profiling?
            if( mLogsMod==LzAsync::LogsMod::END)
            {
                time_point lTimePoint = LzServices::Chrono_Now();
                Jobs lJob(lTimePoint, sT0, "Quitting", "Empty quitting");
                addLogEntry(std::this_thread::get_id(), lJob, false);
            }
            else
            if (mLogsMod==LzAsync::LogsMod::REALTIME)
            {
                time_point lTimePoint = LzServices::Chrono_Now();
                Jobs lJob(lTimePoint, sT0, "Quitting", "Empty quitting");
                UpdateQuittingJsonJob(lJob);
            }

            // Exit thread loop
            break;
        }

        // Event 3: Idle quit?
        //
        // DQ must be empty AND no other worker should be running (since a running worker
        // might push a new job requiring this worker to still be active).
        //
        if( mIdleQuit && mQ.size()==0 && mNbWorksInProgress==0 )
        {
            // Log
            LzLogM("DQ worker", "Idle quitting")

            // Profiling?
            if( mLogsMod==LzAsync::LogsMod::END )
            {
                time_point lTimePoint = LzServices::Chrono_Now();
                Jobs lJob(lTimePoint, sT0, "Quitting", "Idle quitting");
                addLogEntry(std::this_thread::get_id(), lJob, false);
            }
            else
            if (mLogsMod==LzAsync::LogsMod::REALTIME)
            {
                time_point lTimePoint = LzServices::Chrono_Now();
                Jobs lJob(lTimePoint, sT0, "Quitting", "Idle quitting");
                UpdateQuittingJsonJob(lJob);
            }

            // Exit thread loop
            break;
        }

        // Event 0: Queue is NOT empty
        if( mQ.size() != 0 )
        {
            // Pop job

            Jobs lJob = mQ.front(); 
            mQ.pop();

            // Unlock now that we're done messing with the queue
            lLock.unlock();

            // Update the number of works in progress
            // No need to notify since by increasing this value it cannot reach 0
            mNbWorksInProgress++;

            // Profiling?
            if ( mLogsMod==LzAsync::LogsMod::END || mLogsMod==LzAsync::LogsMod::REALTIME)
            {
                // add mCurrentJobs the actually job
                std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
                mCurrentJobs[std::this_thread::get_id()] = &lJob;
            }

            time_point lStart;
            time_point lStop;

            try
            {
                // Profiling?
                if( mLogsMod==LzAsync::LogsMod::END || mLogsMod==LzAsync::LogsMod::REALTIME)
                {
                    lStart = LzServices::Chrono_Now();
                    lJob.JobStart(sT0, lStart);

                    if (mLogsMod==LzAsync::LogsMod::REALTIME)
                    {
                        UpdateJsonJobStart(lJob);
                    }
                }

                // Do the job
                lJob.GetFunc()();
                
                // Update log book
                if (mLogsMod==LzAsync::LogsMod::END || mLogsMod==LzAsync::LogsMod::REALTIME)
                {
                    lStop = LzServices::Chrono_Now();

                    {
                        // now it is not a current job
                        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
                        mCurrentJobs.erase(std::this_thread::get_id());
                    }

                    lJob.JobDone(sT0, lStop, true);

                    if (mLogsMod==LzAsync::LogsMod::END)
                    {
                        addLogEntry(std::this_thread::get_id(), lJob, true);
                    }
                    else // Real time
                    {
                        UpdateJsonJobStop(lJob);
                    }
                }
  
            }
            catch( const std::exception & pExc )
            {
                // Log
                LzLogE("DQ worker", "DispatchQ '"<<mName<<"' encountered an exception while running job '"<<lJob.GetJobName()<<"'!\nException: '"<< pExc.what()<<"'!")

                if (mLogsMod==LzAsync::LogsMod::END || mLogsMod==LzAsync::LogsMod::REALTIME)
                {
                    lStop = LzServices::Chrono_Now();
                    {
                        // now it is not a current job
                        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
                        mCurrentJobs.erase(std::this_thread::get_id());
                    }

                    lJob.JobDone(sT0, lStop, false, pExc.what());
                    
                    // Update log book
                    if( mLogsMod==LzAsync::LogsMod::END )
                    {
                        addLogEntry(std::this_thread::get_id(), lJob, true);
                    }
                    else    // Real time
                    {
                        UpdateJsonJobStop(lJob);
                    }
                }
                

            }
            catch( ... )
            {
                // Log
                LzLogE("DQ worker", "DispatchQ '"<<mName<<"' encountered an exception while running job '"<<lJob.GetJobName()<<"'! Unknown exception.")

                if (mLogsMod==LzAsync::LogsMod::END || mLogsMod==LzAsync::LogsMod::REALTIME)
                {
                    lStop = LzServices::Chrono_Now();
                    {
                        // now it is not a current job
                        std::unique_lock<std::mutex> lock(mCurrentJobsMutex);
                        mCurrentJobs.erase(std::this_thread::get_id());
                    }
                    
                    lJob.JobDone(sT0, lStop, false, "Unknown exception.");    

                    // Update log book
                    if( mLogsMod==LzAsync::LogsMod::END )
                    {
                        addLogEntry(std::this_thread::get_id(), lJob, true);
                    }
                    else // Real time
                    {
                        UpdateJsonJobStop(lJob);
                    }
                }

            }

            // Short cut (optional): skip lock acquisition if asked to mForceQuit
            if( mForceQuit )
            {

                LzLogM("DQ worker", "Force quitting")

                // Profiling?
                if( mLogsMod==LzAsync::LogsMod::END )
                {
                    time_point lTimePoint = LzServices::Chrono_Now();
                    Jobs lJob(lTimePoint, sT0, "Quitting", "Force quitting (2)");
    
                    addLogEntry(std::this_thread::get_id(), lJob, false);
                }
                else
                if (mLogsMod==LzAsync::LogsMod::REALTIME)
                {
                    time_point lTimePoint = LzServices::Chrono_Now();
                    Jobs lJob(lTimePoint, sT0, "Quitting", "Force quitting (2)");
                    UpdateQuittingJsonJob(lJob);
                }

                // Exit thread loop
                break;
            }

            // The update of the nb of works in progress must be done WHILE HOLDING THE MUTEX to ensure
            // that no thread can start to wait on mCV after checking an outdated value of mNbWorksInProgress
            //
            // In other words, the instruction 'while( predicate )' is NOT ATOMIC since it is equivalent to
            // 'while( !predicate ) wait()'. Its 'atomicity' with respect to the value changes in mNbWorksInProgress
            // must thus be ensured using the mutex. Making mNbWorksInProgress atomic is not enough.
            //
            // Update the number of works in progress
            mNbWorksInProgress--;

            // Need to notify all since the nb of works in progress has decreased
            
            mCV.notify_all();
        }
    }
}

}
