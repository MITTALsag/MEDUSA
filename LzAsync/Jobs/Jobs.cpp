#include "Jobs.h"



int64_t Jobs::smNbJob = 0;  ///< Number of jobs, default : 0


//================================================================================
Jobs::Jobs(const std::string& pnameJob, std::function<void(void)> pFunc, 
              int64_t pParentJobId)
    : mJobName(pnameJob), 
      mFunc(pFunc),
      mParentJobId(pParentJobId),
      mJobId(smNbJob++)
{}

//================================================================================
Jobs::Jobs(time_point pTime, time_point pT0, const string& pNameJob, const string& pError /*=""*/)
    : mJobName(pNameJob),
      mError(pError),
      mSucces(true),
      mJobId(-1),
      mParentJobId(-1)
{
  mStart_msec = std::chrono::duration_cast<std::chrono::milliseconds>(pTime - pT0).count();
  mStop_msec = std::chrono::duration_cast<std::chrono::milliseconds>(pTime - pT0).count();
}

//================================================================================
void Jobs::JobDone(const time_point& pT0, time_point pStop, bool pSucces, const std::string& pError)
{
    // Convert the timestamps to relative milliseconds
    mStop_msec =  std::chrono::duration_cast<std::chrono::milliseconds>(pStop  - pT0).count();
    mSucces = pSucces;
    mError = pError;
}


//================================================================================

void Jobs::JobStart(const time_point& pT0, time_point pStart)
{
  // Convert the timestamps to relative milliseconds
  mStart_msec = std::chrono::duration_cast<std::chrono::milliseconds>(pStart - pT0).count();

}
