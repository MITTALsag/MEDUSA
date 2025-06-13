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

// Test recursive DQ -- MESHING
void TestRecursiveDQ__Meshing__5();
void TestRecursiveDQ__Meshing__267();
void TestRecursiveDQ__Meshing__793();
void TestRecursiveDQ__Meshing__2409();
void TestRecursiveDQ__Meshing__Circle_152();
void TestRecursiveDQ__Meshing__Circle_501();

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

