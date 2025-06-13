#include "LzAsync/DispatchQ_UnitTest.h"
#include "LzAsync/DispatchQ.h"



int main()
{

    LzAsync::VerySimpleTest();
    LzAsync::SimpleTest();

    

    LzAsync::MutexTuto();
    LzAsync::SimpleMutexTest();
    LzAsync::TestMutexInMutex();
    LzAsync::TestMutexInMutex2();
    LzAsync::TestMutexInMutex3();
    

    LzAsync::TestTagSimple();
    LzAsync::TestTagInTag();


    // Comments because it exit the code
    // LzAsync::TestDeadlock1WithDispatchQ() ;
    // LzAsync::TestDeadloc2kWithDispatchQ() ;
    // LzAsync::TestDeadlocSomeTimes();

    LzAsync::DispatchQ_UnitTest_ALL();




}



