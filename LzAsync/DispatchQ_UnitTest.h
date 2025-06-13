#pragma once


namespace LzAsync
{
// Run all tests
void DispatchQ_UnitTest_ALL();

// Test 1: missing notification if not using mutex
void DispatchQ_UnitTest_1( bool pShouldLock );

// Test 2: test DispatchQ
void DispatchQ_UnitTest_2();

// Test recursive DQ
void TestRecursiveDQ__Counter();

void TestDeadlock1WithDispatchQ();
void TestDeadloc2kWithDispatchQ();

void VerySimpleTest();
void SimpleTest();

void TestTagSimple();
void TestTagInTag();

void SimpleMutexTest();
void MutexTuto();

void TestMutexInMutex();
void TestMutexInMutex2();
void TestMutexInMutex3();

void TestDeadlocSomeTimes();





}

