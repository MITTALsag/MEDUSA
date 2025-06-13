#include "DispatchQ_UnitTest.h"
#include "DispatchQ_UnitTestData.h"
#include "DispatchQ.h"

#include <chrono>
#include <mutex>


namespace LzAsync
{
using namespace std::chrono_literals;


//================================================================================
void DispatchQ_UnitTest_ALL()
{
    // DispatchQ_UnitTest_1( /*pShouldLock=*/false );
    // DispatchQ_UnitTest_2();
    TestRecursiveDQ__Counter();
}

//================================================================================
void DispatchQ_UnitTest_1( bool pShouldLock )
{
    std::condition_variable CV;
    std::mutex MX;
    std::atomic_bool F{ false };

    // Thread 1
    auto Waiter = [&]
    {
        auto _Pred = [&]{ return F.load(); };
        std::unique_lock<mutex> _Lck(MX);

        // Equivalent implementation of wait( predicate )
#if 0
        CV.wait( _Lck, _Pred );
#else
        while( !_Pred() )
        {
            LzLogM("", "Waiting....")
            std::this_thread::sleep_for( 200ms );
            CV.wait(_Lck);
        }
#endif

        LzLogM("", "Waiter done")
    };


    // Thread 2
    auto Notifier = [&]
    {
        if( pShouldLock )
        {
            // No mutex
            LzLogM("", "F false")
            F = true;
            LzLogM("", "F true")
        }
        else
        {
            // Yes mutex
            std::unique_lock<mutex> L(MX);
            LzLogM("", "F false")
            F = true;
            LzLogM("", "F true")
            L.unlock();
        }

        CV.notify_one();

        LzLogM("", "Notifier done")
    };


    // Test loop
    for( int t=0 ; t<100 ; t++ )
    {
        std::thread T1( Waiter );
        std::thread T2( Notifier );

        T1.join();
        T2.join();
    }
}

// //================================================================================
void DispatchQ_UnitTest_2()
{
    const int NB_TESTS = 1;
    const int NB_JOBS = 20;

    std::atomic<size_t> lNb_Success{ 0 };
    std::atomic<size_t> lNb_Errors{ 0 };

    std::uniform_int_distribution<uint16_t> lShortSleep(0, 50);
    std::uniform_int_distribution<uint16_t> lLongSleep(100, 1000);


    LzAsync::DispatchQ TDQ("TEST DQ", 3);
    TDQ.startProfiling("./json_files/UnitTest2");

    // Run test
    for( int t=0 ; t<NB_TESTS ; t++ )
    {
        LzLogM("", "======================================================= TEST "<<t)

        // Push new jobs
        for( int j=0 ; j<NB_JOBS ; j++ )
        {
            const string lLabel = "Job "+std::to_string(j)+" in test"+std::to_string(t);
            auto job = [&]
                {
                    // Ping
                    LzLogN("DQ worker ", lLabel)

                    // Sleep
                    const uint16_t lMSecs = lShortSleep(LzServices::RandomEngine());
                    std::this_thread::sleep_for(std::chrono::milliseconds(lMSecs));

                    // Boum?
                    if( lMSecs % 20 == 0 )
                    {
                        lNb_Errors++;
                        throw std::runtime_error("*** Fatal error ***");
                    }

                    // Big boum?
                    if( lMSecs % 30 == 0 )
                    {
                        lNb_Errors++;
                        throw 666;
                    }

                    // Success
                    lNb_Success++;
                };

            // Dispatch
            TDQ.dispatch( job, lLabel );
        }

        const uint16_t lMSecs = lLongSleep(LzServices::RandomEngine());

        // Test behavior
        if( lMSecs % 4 == 0 )
        {
            // Log
            LzLogM("", "------------------------------- COMPLETE, "<<lMSecs<<" msec")

            // Wait for completion, and restart pool
            TDQ.completeJobs(true);
            for( int j=0 ; j<50 ; j++ )
            {
                lNb_Success++;
                TDQ.dispatch([]
                    {
                        LzLogM("", "Dummy..."); 
                        std::this_thread::sleep_for( std::chrono::milliseconds(20) ); 
                    }, "Dummy job after resize");
            }
        }
        else
        if( lMSecs % 4 == 1 )
        {
            // Log
            LzLogM("", "------------------------------- INTERRUPT, "<<lMSecs<<" msec")

            // Wait for interruption, and restart pool
            std::this_thread::sleep_for(std::chrono::milliseconds(lMSecs));
            TDQ.interruptJobs(true);
        }
        else
        if( lMSecs % 4 == 2 )
        {
            // Log
            LzLogM("", "------------------------------- NUKE ALL, "<<lMSecs<<" msec")

            // Wait for interruption and termination
            std::this_thread::sleep_for(std::chrono::milliseconds(lMSecs));
            TDQ.interruptJobs(false);

            // Restart pool with 3 workers
            lNb_Success = 0;
            lNb_Errors = 0;
            TDQ.restartPool( 3 );
            TDQ.startProfiling("./json_files/UnitTest2NukeAllAfter");
        }
        else
        {
            // Log
            LzLogM("", "------------------------------- RESIZE, "<<lMSecs<<" msec")

            // Interrupt and restart pool
            std::this_thread::sleep_for(std::chrono::milliseconds(lMSecs));
            // stop profiling in there so call again startProfiling
            TDQ.restartPool( 10 );
            TDQ.startProfiling("./json_files/UnitTest2ResizeAfter");

            //clean variables because we clean logs (and clean job done before)
            lNb_Success = 0;
            lNb_Errors = 0;

            // Do some dummy shit
            for( int j=0 ; j<50 ; j++ )
            {
                lNb_Success++;
                TDQ.dispatch([]
                    {
                        LzLogM("", "Dummy..."); 
                        std::this_thread::sleep_for( std::chrono::milliseconds(20) ); 
                    }, "Dummy job after resize");
            }

            // Wait for completion, no need to restart
            TDQ.completeJobs(false);

        }
    }

    // Shut down DQ
    TDQ.interruptJobs(false);

    // Log final entries
    {
        //Count OK/KO
        const int lCount_Success = TDQ.GetNbSuccess();
        const int lCount_Errors = TDQ.GetNbErrors();

        // Log all
        LzLogN("", "Log entries for '"<<TDQ.getName()<<"'. "<< lCount_Errors <<" error(s).")
        // for( size_t e=0 ; e<lEnts.size() ; e++ )
        //     LzLogM("", e<<" "<<lEnts[e])


        // ASSERT equal
        // if( lEnts.size() != lCount_Success + lCount_Errors )
        //     LzLogException("", "Total entries count error! ("<<lEnts.size()<<" != "<<lCount_Success<<" + "<<lCount_Errors<<").")

        if( lNb_Success != lCount_Success )
            LzLogException("", "Success entries count error! ("<<lNb_Success<<" != "<<lCount_Success<<").")

        if( lNb_Errors != lCount_Errors )
            LzLogException("", "Errors entries count error! ("<<lNb_Errors<<" != "<<lCount_Errors<<").")

        // Clear all
        // TDQ.clearLogEntries();

        // ASSERT empty
        // lEnts = TDQ.getLogEntries();
        // if( lEnts.size() != 0 )
        //     LzLogException("", "Empty entries count error!")
    }
}

//================================================================================
#define TEST_ASYNC_DQ_RECURSION
#ifdef TEST_ASYNC_DQ_RECURSION
static std::shared_ptr<DispatchQ> spTestAsyncDQ;
static std::mutex sTestAsyncIntMex;
#endif
static constexpr bool sTestAsyncLog = false;

//================================================================================
static void RecursiveCount( const std::shared_ptr<vector<int>> pInts, int & pSum )
{
    const size_t N = pInts->size();

    //--------------------------------------
    // Reached a leaf?
    //--------------------------------------
    if( N < 3 )
    {
#ifdef TEST_ASYNC_DQ_RECURSION
        // Acquire mutex
        //logs mutex

        LzAsync::MutexNeed(*spTestAsyncDQ,  "mutex : sTestAsyncIntMex");
        std::unique_lock<std::mutex> lLock( sTestAsyncIntMex );
        LzAsync::MutexHave(*spTestAsyncDQ,  "mutex : sTestAsyncIntMex");

        
#endif
        if( sTestAsyncLog )
            LzLogM("", "Leaf. N= "<<N)

        for( size_t i=0 ; i<N ; i++ )
            pSum += (*pInts)[i];
    
        // std::this_thread::sleep_for( std::chrono::milliseconds(20) );
        LzAsync::MutexFree(*spTestAsyncDQ, "mutex : sTestAsyncIntMex");
        
    }
    else
    {
        // Look for split


        // Compute split
        //const size_t lSplit = N / 2;
        const size_t lSplit = N / 3;

        TagStart(*spTestAsyncDQ, "Balise Test");
        // std::this_thread::sleep_for( std::chrono::milliseconds(20) );


        if( sTestAsyncLog )
            LzLogM("", "Splitting. N= "<<N<<" ==> "<<lSplit)

        // Split data
        std::shared_ptr<vector<int>> lInts_A = std::make_shared<vector<int>>( lSplit );
        for( size_t i=0 ; i<lSplit ; i++ )
            (*lInts_A)[i] = (*pInts)[i];
        //
        std::shared_ptr<vector<int>> lInts_B = std::make_shared<vector<int>>( N - lSplit );
        for( size_t i=lSplit ; i<N ; i++ )
            (*lInts_B)[i - lSplit] = (*pInts)[i];

    TagStop(*spTestAsyncDQ, "Balise Test");

#ifdef TEST_ASYNC_DQ_RECURSION
        // Async
        spTestAsyncDQ->dispatch( [lInts_A, &pSum] { RecursiveCount( lInts_A, pSum ); }, "RecursiveCount : with size = "+std::to_string(lInts_A->size()) );
        spTestAsyncDQ->dispatch( [lInts_B, &pSum] { RecursiveCount( lInts_B, pSum ); }, "RecursiveCount : with size = "+std::to_string(lInts_B->size()) );
#else
        // Sync
        RecursiveCount( lInts_A, pSum );
        RecursiveCount( lInts_B, pSum );
#endif
    }
}

//=================================================================================
void TestRecursiveDQ__Counter()
{
    // Init a vector of ints
    const size_t N = 1000;
    std::shared_ptr<vector<int>> lInts = std::make_shared<vector<int>>( N );
    for( size_t i=0 ; i<N ; i++ )
        (*lInts)[i] = i+1;

    // Compute expected result
    const int SUM = N*(N + 1)/2;

#ifdef TEST_ASYNC_DQ_RECURSION
    //--------------------------------------
    // Set-up async mode
    //--------------------------------------

    // Create DQ
    spTestAsyncDQ = std::make_shared<DispatchQ>( "TestAsync -- DQ", 8 );
#endif

    //--------------------------------------
    // Launch recursive counting
    //--------------------------------------

    {
        LzLogTimeN("", "Counting ints")

        int lSum = 0;

#ifdef TEST_ASYNC_DQ_RECURSION

        spTestAsyncDQ->startProfiling("./json_files/TestRecursiveDQCounter2", LzAsync::LogsMod::REALTIME);
        // Count ints
        spTestAsyncDQ->dispatch( [&]{ RecursiveCount( lInts, lSum ); }, "RecursiveCount : with size = "+std::to_string(lInts->size()) );

        // Wait for all jobs to finish
        spTestAsyncDQ->waitForIdle( /*pRestartPool=*/false );

        // spTestAsyncDQ->ExportLogsJson("./json_files/TestRecursiveDQCounter.json");


        // Check errors
        // vector<string> lDQErrs = spTestAsyncDQ->getLogEntries();
        const size_t lNbErrs = spTestAsyncDQ->GetNbErrors();
        if( lNbErrs != 0 )
        {


            // for(const string & lErr := lDQErrs)
            //     LzLogE("","lErr")


            // Log
            LzLogN("", "Found "<<lNbErrs<<" error(s) during async test:")
            // Fail
            LzLogException("", "Found "<<lNbErrs<<" error(s) during async test!")
        }
#else
        RecursiveCount( lInts, lSum );
#endif

        // Check result
        if( lSum != SUM )
            LzLogException("", "Unexpected sum found! "<<lSum<<" != "<<SUM<<".")
        else
            LzLogM("", "Summation is correct.")
    }
}

//=================================================================================
void TestDeadlock1WithDispatchQ() 
{
    // Creat a 2 thread Dispatch queue
    LzAsync::DispatchQ dispatchQ("DeadlockTestQ", 2);

    // Share Mutexes that will cause the deadlock
    std::mutex mutexA, mutexB;

    // Job 1 : take mutexA then try to take mutexB
    auto job1 = [&]() {
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simulate some work

        // Try to take mutexB but block because job2 have it
        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");


        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexB");
    };

    // Job 2 : take mutexB then try to take mutexA
    auto job2 = [&]() {
        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simulate some work

        // Try to take mutexA but block because job1 have it
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        MutexFree(dispatchQ, "mutexB");
        MutexFree(dispatchQ, "mutexA");
    };


    dispatchQ.startProfiling("./json_files/Deadlock1RealTimeTest", LzAsync::LogsMod::REALTIME);

    dispatchQ.dispatch(job1, "Job1_Deadlock");
    dispatchQ.dispatch(job2, "Job2_Deadlock");

    dispatchQ.interruptJobs(false);
    std::cout << "this line will never be printed" << std::endl;

}

//=================================================================================

void TestDeadloc2kWithDispatchQ() 
{
    // Crée une DispatchQ avec 4 threads (pour exécuter les jobs en parallèle)
    LzAsync::DispatchQ dispatchQ("DeadlockTestQ", 4);

    // Mutex partagés qui causeront le deadlock
    std::mutex mutexA, mutexB, mutexC, mutexD, mutexE;

    // Job 1 : Prend mutexA puis essaie de prendre mutexB
    auto job1 = [&]() {
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        std::cout << "Job 1 a pris mutexA, attend mutexB..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simule un délai

        // Essaie de prendre mutexB (mais bloqué car job2 l'a déjà)
        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");
        std::cout << "Job 1 a pris mutexB (NE SERA JAMAIS AFFICHÉ)" << std::endl;

        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexB");
    };

    // Job 2 : Prend mutexB puis essaie de prendre mutexA
    auto job2 = [&]() {
        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        std::cout << "Job 2 a pris mutexB, attend mutexA..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simule un délai

        // Essaie de prendre mutexA (mais bloqué car job1 l'a déjà)
        MutexNeed(dispatchQ, "mutexC");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexC");

        std::cout << "Job 2 a pris mutexA (NE SERA JAMAIS AFFICHÉ)" << std::endl;

        MutexFree(dispatchQ, "mutexB");
        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexC");


    };

    // Job 1 : Prend mutexA puis essaie de prendre mutexB
    auto job3 = [&]() {
        MutexNeed(dispatchQ, "mutexC");
        std::unique_lock<std::mutex> lockC(mutexC);
        MutexHave(dispatchQ, "mutexC");
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        MutexNeed(dispatchQ, "mutexD");
        std::unique_lock<std::mutex> lockD(mutexD);
        MutexHave(dispatchQ, "mutexD");

        std::cout << "Job 3 a pris mutexC et mutexD, attend mutexE..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simule un délai

        // Essaie de prendre mutexB (mais bloqué car job2 l'a déjà)
        MutexNeed(dispatchQ, "mutexE");
        std::unique_lock<std::mutex> lockE(mutexE);
        MutexHave(dispatchQ, "mutexE");
        std::cout << "Job 3 a pris mutexE (NE SERA JAMAIS AFFICHÉ)" << std::endl;

        MutexFree(dispatchQ, "mutexC");
        MutexFree(dispatchQ, "mutexD");
        MutexFree(dispatchQ, "mutexE");
    };

    // Job 1 : Prend mutexA puis essaie de prendre mutexB
    auto job4 = [&]() {
        MutexNeed(dispatchQ, "mutexE");
        std::unique_lock<std::mutex> lockE(mutexE);
        MutexHave(dispatchQ, "mutexE");

        std::cout << "Job 4 a pris mutexE, attend mutexA..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simule un délai

        // Essaie de prendre mutexB (mais bloqué car job2 l'a déjà)
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");
        std::cout << "Job 4 a pris mutexB (NE SERA JAMAIS AFFICHÉ)" << std::endl;

        MutexFree(dispatchQ, "mutexE");
        MutexFree(dispatchQ, "mutexA");
    };



    dispatchQ.startProfiling("./json_files/Deadlock2RealTimeTest", LzAsync::LogsMod::REALTIME);
    // Envoie les jobs dans la DispatchQ (exécution en parallèle)
    dispatchQ.dispatch(job1, "Job1_Deadlock");
    dispatchQ.dispatch(job2, "Job2_Deadlock");
    dispatchQ.dispatch(job3, "Job3_Deadlock");
    dispatchQ.dispatch(job4, "Job4_Deadlock");
    // Attend que les jobs soient terminés (ne terminera jamais à cause du deadlock)
    dispatchQ.interruptJobs(false);
    std::cout << "Cette ligne ne sera jamais atteinte !" << std::endl;

}

//=================================================================================
void VerySimpleTest()
{
    LzAsync::DispatchQ dispatchQ("VerySimpleTest", 3);

    auto Smalljob = [&]()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    };
    
    auto Bigjob = [&]()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(75));
        dispatchQ.dispatch(Smalljob, "Dispatch 6");
        dispatchQ.dispatch(Smalljob, "Dispatch 7");
    };

    auto ErrorsJob = [&]()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(25));
        throw std::runtime_error("Very Simple Test: This is a test.");
    };

    auto job2 = [&]()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        dispatchQ.dispatch(ErrorsJob, "Dispatch 2");
        dispatchQ.dispatch(Smalljob, "Dispatch 3");
        dispatchQ.dispatch(ErrorsJob, "Dispatch 4");
        dispatchQ.dispatch(Bigjob, "Dispatch 5");

    };

    dispatchQ.startProfiling("./json_files/VerySimpleTest");
    dispatchQ.dispatch(job2, "Dispatch 1");
    dispatchQ.waitForIdle(false);


}

//=================================================================================

void SimpleTest()
{
    // Create a new dispatch queue
    LzAsync::DispatchQ dispatchQ("SimpleTest", 3);

    bool pInstant = true;
    int lcount = 0;

    std::function<void(void)> lJobDispatch;
    lJobDispatch = [&]()
    {
        bool localInstant = pInstant;

        dispatchQ.dispatch([localInstant]() {
            if (!localInstant)
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }, "Job OK " + std::to_string(lcount));

        dispatchQ.dispatch([localInstant]() {
            if (!localInstant)
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            throw std::runtime_error("Simple Test: This is a test.");
        }, "Job KO " + std::to_string(lcount));

        lcount++;
        pInstant = !pInstant;

        if (lcount < 10)
            lJobDispatch();
    };

    dispatchQ.startProfiling("./json_files/SimpleTest", LzAsync::LogsMod::REALTIME);
    lJobDispatch();

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);
}

//=================================================================================
void TestTagSimple()
{
    // Create a new dispatch queue
    LzAsync::DispatchQ dispatchQ("TestTagSimple", 1);


    auto jobTag = [&]() 
    {
        for (int i = 0; i< 10 ; i++)
        {   
            // Simulate some work
            std::this_thread::sleep_for(std::chrono::milliseconds(10));


            TagStart(dispatchQ, "Tag");
            {
                // Simulate some work
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
            TagStop(dispatchQ, "Tag");

        }

    };


    // Start profiling
    dispatchQ.startProfiling("./json_files/TestTagSimple");

    // Dispatch a job with a tag
    dispatchQ.dispatch(jobTag, "Job with Tag");

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);
}


//=================================================================================
void TestTagInTag()
{
    // Create a new dispatch queue
    LzAsync::DispatchQ dispatchQ("TestTagInTag", 1);


    auto jobTag = [&]() 
    {
        TagStart(dispatchQ, "First tag");
        {   
            // Simulate some work
            std::this_thread::sleep_for(std::chrono::milliseconds(100));

            TagStart(dispatchQ, "Second tag");
            {
                // Simulate some work
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            TagStop(dispatchQ, "Second tag");

            // Simulate some work
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        TagStop(dispatchQ, "First tag");

    };


    // Start profiling
    dispatchQ.startProfiling("./json_files/TestTagInTag");

    // Dispatch a job with a tag
    dispatchQ.dispatch(jobTag, "Job with Tag");

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);
}

//=================================================================================

void SimpleMutexTest() 
{
    LzAsync::DispatchQ queue("SimpleMutex", 2);
    queue.startProfiling("./json_files/SimpleMutex");
    
    std::mutex myMutex;

    auto jobMutex = [&] 
    {
        // Log before lock
        MutexNeed(queue, "MyCounterMutex");  
        std::lock_guard<std::mutex> lock(myMutex);
        // Log after lock
        MutexHave(queue, "MyCounterMutex");  
        
        // Simulate some work
        std::this_thread::sleep_for( std::chrono::milliseconds(100) );
         
        // Log after the unlock
        MutexFree(queue, "MyCounterMutex"); 
    };


    // Jobs with mutex
    queue.dispatch(jobMutex, "IncrementTask 1");
    queue.dispatch(jobMutex, "IncrementTask 2");

    queue.waitForIdle(false);

}

//=================================================================================
void MutexTuto()
{
    std::mutex mutexA;
    LzAsync::DispatchQ dispatchQ("MutexTuto", 6);

    auto jobMutex = [&]() 
    {
        MutexNeed(dispatchQ, "mutex");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutex");

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(20));

        MutexFree(dispatchQ, "mutex");

    };

    // Start profiling
    dispatchQ.startProfiling("./json_files/TestMutexTuto");

    // Dispatch a job with a tag
    dispatchQ.dispatch(jobMutex, "Job with Mutex 1");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 2");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 3");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 4");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 5");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 6");

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);
}


//=================================================================================
void TestMutexInMutex()
{
    // Create a new dispatch queue
    LzAsync::DispatchQ dispatchQ("TestMutexInMutex", 1);

    std::mutex mutexA, mutexB;

    auto jobMutex = [&]() 
    {
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexB");
    };

    // Start profiling
    dispatchQ.startProfiling("./json_files/TestMutexInMutex");

    // Dispatch a job with a tag
    dispatchQ.dispatch(jobMutex, "Job with Mutex");

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);

}

//=================================================================================
void TestMutexInMutex2()
{
    // Create a new dispatch queue
    LzAsync::DispatchQ dispatchQ("TestMutexInMutex2", 2);

    std::mutex mutexA, mutexB;

    auto jobMutex = [&]() 
    {
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexB");
    };

    // Start profiling
    dispatchQ.startProfiling("./json_files/TestMutexInMutex2");

    // Dispatch a job with a tag
    dispatchQ.dispatch(jobMutex, "Job with Mutex 1");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 2");

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);

}


//=================================================================================
void TestMutexInMutex3()
{
    // Create a new dispatch queue
    LzAsync::DispatchQ dispatchQ("TestMutexInMutex3", 2);

    std::mutex mutexA, mutexB;

    auto jobMutex = [&]() 
    {
        {
            MutexNeed(dispatchQ, "mutexA");
            std::unique_lock<std::mutex> lockA(mutexA);
            MutexHave(dispatchQ, "mutexA");

            // Simulate some work
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            MutexFree(dispatchQ, "mutexA");

        }

        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        MutexFree(dispatchQ, "mutexB");
    };

    // Start profiling
    dispatchQ.startProfiling("./json_files/TestMutexInMutex3");

    // Dispatch a job with a tag
    dispatchQ.dispatch(jobMutex, "Job with Mutex 1");
    dispatchQ.dispatch(jobMutex, "Job with Mutex 2");

    // Wait for all jobs to finish
    dispatchQ.waitForIdle(false);
}

//=================================================================================

void TestDeadlocSomeTimes() 
{
    // Crée une DispatchQ avec 4 threads (pour exécuter les jobs en parallèle)
    LzAsync::DispatchQ dispatchQ("DeadlockTestQ", 3);

    // Mutex partagés qui causeront le deadlock
    std::mutex mutexA, mutexB, mutexC;

    // Job 1 : takes mutex A then B then C
    auto job1 = [&]() 
    {
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Simule un délai

        // Essaie de prendre mutexB (mais bloqué car job2 l'a déjà)
        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Simule un délai

        MutexNeed(dispatchQ, "mutexC");
        std::unique_lock<std::mutex> lockC(mutexC);
        MutexHave(dispatchQ, "mutexC");

        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Simule un délai

        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexB");
        MutexFree(dispatchQ, "mutexC");
    };

    // Job 2 : takes mutex A then C then B
    auto job2 = [&]() {
        MutexNeed(dispatchQ, "mutexA");
        std::unique_lock<std::mutex> lockA(mutexA);
        MutexHave(dispatchQ, "mutexA");

        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Simule un délai

        MutexNeed(dispatchQ, "mutexC");
        std::unique_lock<std::mutex> lockC(mutexC);
        MutexHave(dispatchQ, "mutexC");

        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Simule un délai


        // Essaie de prendre mutexB (mais bloqué car job2 l'a déjà)
        MutexNeed(dispatchQ, "mutexB");
        std::unique_lock<std::mutex> lockB(mutexB);
        MutexHave(dispatchQ, "mutexB");

        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Simule un délai

        MutexFree(dispatchQ, "mutexA");
        MutexFree(dispatchQ, "mutexC");
        MutexFree(dispatchQ, "mutexB");


    };


    dispatchQ.startProfiling("./json_files/DeadlocksomeTimes", LzAsync::LogsMod::REALTIME);
    // Envoie les jobs dans la DispatchQ (exécution en parallèle)
    dispatchQ.dispatch(job1, "Job1_Deadlock");
    dispatchQ.dispatch(job2, "Job2_Deadlock");
    dispatchQ.dispatch(job1, "Job1_Deadlock 2");
    // Attend que les jobs soient terminés (ne terminera jamais à cause du deadlock)
    dispatchQ.waitForIdle(false);
    std::cout << "Cette ligne ne sera jamais atteinte !" << std::endl;

}


} // namespace LzAsync
