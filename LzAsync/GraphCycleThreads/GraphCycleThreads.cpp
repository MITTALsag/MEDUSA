#include "GraphCycleThreads.h"

void GraphCycleThreads::initialize(vector<std::thread>& pThreads, time_point sT0)
{
    mNbVerticies = pThreads.size();
    mMapThreadInt.clear();
    mMapIntThread.clear();
    mMapThreadMutexNeed.clear();
    mMapMutexThreadHave.clear();
    mGraph.clear();
    mCyclePath.clear();
    T0 =sT0 ;

    int i = 0;
    for (auto& iThread : pThreads)       // use a reference to avoid copying thread objects
    {
        mMapThreadInt[iThread.get_id()] = i;
        mMapIntThread[i]=iThread.get_id();
        mGraph.push_back(-1);
        i++;
    }
}

//================================================================================

bool GraphCycleThreads::UpdateCycle(const string& pName, const ThreadMutexState& pThreadMutexState)
{
    std::unique_lock<std::mutex> lock(mMexCycle);

    switch (pThreadMutexState)
    {
        case ThreadMutexState::NEED:
        {
            UpdateMutexNeed(std::this_thread::get_id(),pName);
            return CycleDetector(); // Check for cycles after updating the mutex need
            break;
        }
        case ThreadMutexState::HAVE:
        {
            UpdateMutexHave(std::this_thread::get_id(),pName);
            return false;
            break;

        }
        case ThreadMutexState::FREE:
        {
            UpdateMutexFree(std::this_thread::get_id(),pName);
            return false;
            break;
        }
    }
    throw(std::runtime_error("GraphCycleThreads::UpdateCycle: Invalid ThreadMutexState"));
    return false;
}

//================================================================================


void GraphCycleThreads::UpdateMutexHave(std::thread::id pIdThread, string pMutexName)
{
    mMapMutexThreadHave[pMutexName] = mMapThreadInt[pIdThread];
    mMapThreadMutexNeed.erase(mMapThreadInt[pIdThread]); 
    for (int i=0; i<mNbVerticies; i++) //Updating the Graph
    {
        if(mMapThreadMutexNeed[i]==pMutexName)
        {
            AddEdge(i, mMapThreadInt[pIdThread]);
        }
    }
}

//================================================================================


void GraphCycleThreads::UpdateMutexNeed(std::thread::id pIdThread, string pMutexName)
{
    mMapThreadMutexNeed[mMapThreadInt[pIdThread]] = pMutexName ; 
    auto lInd = mMapMutexThreadHave.find(pMutexName);
    if (lInd != mMapMutexThreadHave.end()) //Updating the Graph
    {
        AddEdge(mMapThreadInt[pIdThread],lInd->second);
    }
}

//================================================================================


void GraphCycleThreads::UpdateMutexFree(std::thread::id pIdThread, string pMutexName)
{
    mMapMutexThreadHave.erase(pMutexName);; 
    for (int i=0; i<mNbVerticies; i++)      //Updating the Graph
    {
        if(mMapThreadMutexNeed[i]==pMutexName)
        {
            EraseEdge(i);
        }
    }

}

//================================================================================


void GraphCycleThreads::AddEdge(int pU, int pV) 
{
    mGraph[pU] = pV ;
}

//================================================================================


void GraphCycleThreads::EraseEdge(int pU)
{
    mGraph[pU] = -1;
}

//================================================================================

bool GraphCycleThreads::CycleDetector() //Floyd Algorithm that find cycle in O(nb_verticies)
{ 
    std::vector<bool> visited(mNbVerticies, false);
    for (int i = 0; i < mNbVerticies; ++i) 
    {
        if (!visited[i]) 
        {
            int lSlow = i, lFast = i;
            while (lFast!=-1 && mGraph[lFast] != -1 && mGraph[mGraph[lFast]] != -1) 
            {
                visited[lSlow] = true ;
                visited[lFast] = true ;
                lSlow = mGraph[lSlow] ;
                lFast = mGraph[mGraph[lFast]] ;
                if (lSlow == lFast) //Cycle detected
                {
                    SaveCycle(lSlow) ;
                    return true ;
                }
            }
        }
    }
    return false ;
}

//================================================================================

void GraphCycleThreads::SaveCycle(int pStart)
{
    mCyclePath.clear();
    int it = pStart ;
    int lCount =0;
    do
    {
        if (it==-1)
        {
            throw std::runtime_error("GraphCycleThreads::SaveCycle: Invalid start vertex, it=-1");
            break;
        }
        mCyclePath.push_back(mMapIntThread[it]);
        it = mGraph[it];
        lCount++ ;
    } while ( it != pStart && lCount < mNbVerticies);
}

void GraphCycleThreads::UpdateJsonDeadLock(ordered_json &pJson) //
{
    //std::unique_lock<std::mutex> llock(mLogMex);
    time_point lnow = LzServices::Chrono_Now();
    int64_t lMsecNow = std::chrono::duration_cast<std::chrono::milliseconds>(lnow - T0).count();
    ordered_json lJsonCycle;
    ordered_json lJsonThread;
    vector<string> lMutexHave;
    vector<ordered_json> lJsonThreadVector;
    for (const auto& iThread : mCyclePath) 
    {
        lMutexHave.clear();
        std::ostringstream lSS;
        lSS << iThread; // Conversion of this id into a string
        lJsonThread["ThreadId"] = lSS.str();
        lJsonThread["Need"]=  mMapThreadMutexNeed[mMapThreadInt[iThread]];
        for (const auto& [iMutex, lThread]: mMapMutexThreadHave)
        {
            if(mMapIntThread[lThread]==iThread)lMutexHave.push_back(iMutex);
        }
        lJsonThread["Have"] = lMutexHave;
        lJsonThreadVector.push_back(lJsonThread);
    }
    lJsonCycle["Time"]=lMsecNow;
    lJsonCycle["Cycle"]= lJsonThreadVector ;
    pJson["DeadLock"] = lJsonCycle;
}