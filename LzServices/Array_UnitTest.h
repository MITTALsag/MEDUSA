#pragma once


namespace LzServices
{
// Run all tests
void Array_UnitTest_ALL();

// Test 1: init, set all, read, write, iterate, debug, IdxToCounter, CounterToIdx
void Array_UnitTest_1();

// Test 2: copy, move, compare
void Array_UnitTest_2();

// Test 3: lambda construct, IdxToCounter, CounterToIdx
void Array_UnitTest_3();

// Test 4: upper, lower bounds
void Array_UnitTest_4();

// Test 5: transformations
void Array_UnitTest_5();
}
