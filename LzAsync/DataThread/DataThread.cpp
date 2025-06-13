#include "DataThread.h"
#include <iostream>
#include <mutex>
#include <iomanip>




void TagContainer::SetTimes(const int64_t pTime, const string& pName, const TagState& pTagState)
{
    switch (pTagState) 
    {
        case TagState::START: 
        {
            //Create a new Tag log and set Start time
            Tag lNewTag = {pName, {pTime, -1}};
            mData.push_back(lNewTag); 
            break;
        }
        case TagState::STOP: 
        {

            int it;
            // Reverse loop to find the last tag with the same name
            for (it = mData.size() - 1; it >= 0; --it) 
            {
                Tag& tag = mData[it];
                if (get<0>(tag) == pName) 
                {
                    // Set Stop time
                    get<1>(tag)[1] = pTime; 
                    break;
                }
            }
            if (it < 0) 
            {
                // If no tag found, throw an error
                throw std::runtime_error("In TagContainer::SetTimes : Tag logs for " + pName + " are empty");
            }
            break;
        }
        default:
            throw std::invalid_argument("In TagContainer::SetTimes : Invalid TagState value");
    }
}

//======================================================================================
void MutexContainer::SetTimes(const int64_t pTime, const string& pName, const ThreadMutexState& pHow)
{
    switch (pHow) 
    {
        case ThreadMutexState::NEED: 
        {
            
            //Create a new mutex log and set Need time
            MutexLogs lNewMutexLog = {pTime, -1, -1};
            mData[pName].push_back(lNewMutexLog); 
            break;
        }
        case ThreadMutexState::HAVE: 
        {
            // displayData(mData, "Have");


            // Check if the mutex log exists (must be created in NEED state before)
            auto it = mData.find(pName);
            if (it == mData.end()) 
            {
                throw std::runtime_error("In MutexContainer::SetTimes in HAVE state: Mutex logs for " + pName + " don't exist");
            }

            vector<MutexLogs>& lMutexLogs = it->second;
            if (lMutexLogs.empty()) 
            {
                throw std::runtime_error("In MutexContainer::SetTimes : Mutex logs for " + pName + " are empty");
            }
            // Set Have time (is the last one)
            lMutexLogs.back()[1] = pTime; 
            break;
        }
        case ThreadMutexState::FREE: 
        {
            // Check if the mutex log exists (must be created in NEED and Have states before)
            vector<MutexLogs>& lMutexLogs = mData.at(pName);
            if (lMutexLogs.empty()) 
            {
                throw std::runtime_error("In MutexContainer::SetTimes : Mutex logs for " + pName + " are empty");
            }
            // Set Free time (is the last one)
            lMutexLogs.back()[2] = pTime; 
            break;
        }
        default:
            throw std::invalid_argument("In MutexContainer::SetTimes : Invalid ThreadMutexState value");
    }
}