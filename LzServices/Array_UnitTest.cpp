#include "Array_UnitTest.h"
#include "Array.h"


namespace LzServices
{
//================================================================================
void Array_UnitTest_ALL()
{
    Array_UnitTest_1();
    Array_UnitTest_2();
    Array_UnitTest_3();
    Array_UnitTest_4();
    Array_UnitTest_5();
}

//================================================================================
void Array_UnitTest_1()
{
    // Log
    LzLogN("", "Test 1: init, set all, iterate, debug, IdxToCounter, CounterToIdx")

    // Empty array
    Array<string, short> A("A");

    // TEST
    if( !A.IsEmpty() )
        LzLogException("", "ASSERT FAILED")

    // Log
    A.LogInfo("-----------> Initial empty array A");

    // Test for Arrays with dims 1, 2, ..., lMaxDims
    const size_t lMaxDims = 7;
    for( size_t iNbDims=1 ; iNbDims<=lMaxDims ; iNbDims++ )
    {
        // Compute size
        vector<size_t> lSize( iNbDims );
        for( size_t d=0 ; d<iNbDims ; d++ )
            lSize[d] = d + 1;

        // Resize and set all
        A.Resize( lSize );

        // CHECK
        if( A.Dims() != iNbDims )
            LzLogException("", "ASSERT FAILED")

        // CHECK
        for( size_t d=0 ; d<iNbDims ; d++ )
        {
            if( A.Size(d) != lSize[d] )
                LzLogException("", "ASSERT FAILED")
        }

        // CHECK
        if( A.SizeVec() != lSize )
            LzLogException("", "ASSERT FAILED")

        // Log
        A.LogInfo("-----------> Resized "+A.Name());

        // CHECK
        if( A.IsEmpty() )
            LzLogException("", "ASSERT FAILED")

        // Set all
        A.SetAll( "--A--" );

        // Set all through ctor
        {
            Array<string, short> another_A( A.SizeVec(), "--A--", "another A" );

            // CHECK
            if( another_A != A )
                LzLogException("", "ASSERT FAILED")
        }

        // TEST
        for( size_t c=0 ; c<A.Count() ; c++ )
        {
            if( A.Buffer()[c] != "--A--" )
                LzLogException("", "ASSERT FAILED")
        }

        // Renum
        for( size_t c=0 ; c<A.Count() ; c++ )
            A.Buffer()[c] = A.Buffer()[c] + " " + std::to_string(c);

        // Log
        if( iNbDims <= 4 )
            A.Log();
        else
            LzLogM("", "... Array too big, not logging.")

        // Iterate
        for( vector<short> i=A.InitIdx() ; A.CheckIdx(i) ; A.NextIdx(i) )
        {
            // Compute counter from idx
            const size_t c = A.IdxToCounter( i );

            string lStr;
            {
                if( iNbDims == 1 )
                {
                    // Read
                    lStr = A( i[0] );

                    // Put it back
                    A( i[0] ) = lStr;
                }
                else
                if( iNbDims == 2 )
                {
                    // Read
                    lStr = A( i[0], i[1] );

                    // Put it back
                    A( i[0], i[1] ) = lStr;
                }
                else
                if( iNbDims == 3 )
                {
                    // Read
                    lStr = A( i[0], i[1], i[2] );

                    // Put it back
                    A( i[0], i[1], i[2] ) = lStr;
                }
                else
                {
                    // Read
                    lStr = A( i );

                    // Put it back
                    A( i ) = lStr;
                }
            }

            // CHECK
            if( lStr != "--A-- " + std::to_string(c) )
                LzLogException("", "ASSERT FAILED")

            // CHECK
            if( i != A.CounterToIdx(c) )
                LzLogException("", "ASSERT FAILED")
        }
    }

    // Log
    LzLogM("", "OK")
}

//================================================================================
void Array_UnitTest_2()
{
    // Log
    LzLogN("", "Test 2: copy, move, compare")

    // Move: constructors
    {
        // Create source array
        Array<string, short> A( [](Array<string, short> * pThis)
            {
                pThis->Resize( { 3, 4, 5 } );
                for( auto i=pThis->InitIdx() ; pThis->CheckIdx(i) ; pThis->NextIdx(i) )
                    (*pThis)( i ) = "Slice "+pThis->IdxToString( i );

            }, "A" );

        // Log
        A.Log();

        // Copy construct
        const Array<string, short> copyA( A );

        // Log
        copyA.Log();

        // CHECK
        if( copyA != A )
            LzLogException("", "ASSERT FAILED")

        // CHECK
        if( A.Buffer() == copyA.Buffer() )
            LzLogException("", "ASSERT FAILED")

        // Move construct
        const Array<string, short> moveA( std::move(A) );

        // Log
        moveA.Log();

        // Log
        A.Log();

        // CHECK
        if( !A.IsEmpty() )
            LzLogException("", "ASSERT FAILED")

        // CHECK
        if( moveA != copyA )
            LzLogException("", "ASSERT FAILED")
    }

    // Move: assign
    {
        // Create source array
        Array<string, short> A( [](Array<string, short> * pThis)
            {
                pThis->Resize( { 3, 4, 5 } );
                for( auto i=pThis->InitIdx() ; pThis->CheckIdx(i) ; pThis->NextIdx(i) )
                    (*pThis)( i ) = "Slice "+pThis->IdxToString( i );

            }, "A" );

        // Copy constructor
        Array<string, short> copyA = A;

        // Log
        copyA.Log();

        // Copy assign
        Array<string, short> assign_copyA;
        assign_copyA = A;

        // Log
        assign_copyA.Log();

        // CHECK
        if( assign_copyA != A )
            LzLogException("", "ASSERT FAILED")

        // CHECK
        if( A.IsEmpty() )
            LzLogException("", "ASSERT FAILED")

        // Move assign
        Array<string, short> assign_moveA;
        assign_moveA = std::move( A );

        // Log
        assign_moveA.Log();

        // Log
        A.Log();

        // CHECK
        if( assign_copyA != assign_moveA )
            LzLogException("", "ASSERT FAILED")

        // CHECK
        if( !A.IsEmpty() )
            LzLogException("", "ASSERT FAILED")
    }

    // Compare with empty
    {
        // Both are empty
        Array<long> A("naturally empty");
        Array<long> B( { 1, 2, 0, 4 }, "strangely empty");

        // CHECK
        if( A != B )
            LzLogException("", "ASSERT FAILED")

        // Only one is empty
        B.Resize( { 5 } );

        // CHECK
        if( A == B )
            LzLogException("", "ASSERT FAILED")
    }

    // Compare both non-empty
    {
        // A
        Array<int> A( { 3, 4 }, "array A" );
        for( vector<size_t> i=A.InitIdx() ; A.CheckIdx(i) ; A.NextIdx(i) )
        {
            size_t c = A.IdxToCounter( i );
            LzLogM("", "Setting value in '"<<A.Name()<<"', at "<<A.IdxToString(i)<<" = "<<c)
            A( i ) = c;
        }

        // Log
        A.Log();

        // B
        const Array<int> B( { 3, 4 },
                            { 0,  1,  2,
                              3,  4,  5,
                              6,  7,  8,
                              9, 10, 11 }, "array B" );
        // Log
        B.Log();

        // CHECK
        if( A != B )
            LzLogException("", "ASSERT FAILED")

        // Modify A
        A( 1, 2 ) = 666;

        // CHECK
        if( A == B )
            LzLogException("", "ASSERT FAILED")
    }

    // Log
    LzLogM("", "OK")
}

//================================================================================
void Array_UnitTest_3()
{
    // Log
    LzLogN("", "Test 3: lambda construct, IdxToCounter, CounterToIdx")

    // Stringify helper
    auto Stringify = []( char i, char j, char k, char l ) -> string
        {
            return "<< "+std::to_string(i)+", "
                        +std::to_string(j)+", "
                        +std::to_string(k)+", "
                        +std::to_string(l)+" >>";
        };

    // Lambda constructor
    const Array<string, char> LambdArray( [&Stringify](Array<string, char> * pThis)
        {
            // Resize
            pThis->Resize( { 2, 3, 4, 5 } );

            // Fill
            for( char i=0 ; i<pThis->Size(0) ; i++ )
            for( char j=0 ; j<pThis->Size(1) ; j++ )
            for( char k=0 ; k<pThis->Size(2) ; k++ )
            for( char l=0 ; l<pThis->Size(3) ; l++ )
            {
                // Compute string
                const string lStr = Stringify(i, j, k, l);

                // Set
                (*pThis)( {i, j, k, l} ) = lStr;

                // Idx to counter
                size_t c = pThis->IdxToCounter( {i, j, k, l} );

                // TEST
                if( pThis->Buffer()[c] != lStr )
                    LzLogException("", "ASSERT FAILED")
            }
        },
    "Lambda_ctor");

    // Log
    LambdArray.Log();

    // TEST
    for( size_t c=0 ; c<LambdArray.Count() ; c++ )
    {
        vector<char> lIdx;
        LambdArray.CounterToIdx( c, lIdx );

        if( LambdArray.Buffer()[ c ] != Stringify(lIdx[0], lIdx[1], lIdx[2], lIdx[3]) )
            LzLogException("", "ASSERT FAILED")
    }

    // Log
    LzLogM("", "OK")
}

//================================================================================
void Array_UnitTest_4()
{
    // Log
    LzLogN("", "Test 4: upper, lower bounds")

    // Lambda construct empty
    const Array<double> A( []( Array<double> * pThis ){} );

    // Log
    A.Log();

    // Min and Max of empty array should fail
    try
    {
        // Should raise exception
        A.GetMin();

        // CHECK
        LzLogException("", "ASSERT FAILED")
    }
    catch( ... ) { /* NOP */ }
    //
    try
    {
        // Should raise exception
        A.GetMax();

        // CHECK
        LzLogException("", "ASSERT FAILED")
    }
    catch( ... ) { /* NOP */ }

    // Lambda construct NON EMPTY
#if 1
    const Array<double> B( []( Array<double> * pThis )
        {
            // 2x4 matrix
            pThis->Resize( { 2, 4 } );

            // Init
            (*pThis)(0, 0) = 3.1415;
            (*pThis)(1, 0) = -3.1415;
            //
            (*pThis)(0, 1) = cos(3.1415);
            (*pThis)(1, 1) = sin(3.1415);
            //
            (*pThis)(0, 2) = exp(3.1415);
            (*pThis)(1, 2) = log(3.1415);
            //
            (*pThis)(0, 3) = tan(3.1415);
            (*pThis)(1, 3) = cosh(3.1415);
        } );
#else
    const Array<double> B = []() -> Array<double>
        {
            Array<double> lTmp("Coucou");

            // 2x4 matrix
            lTmp.Resize( { 2, 4 } );

            // Init
            lTmp(0, 0) = 3.1415;
            lTmp(1, 0) = -3.1415;
            //
            lTmp(0, 1) = cos(3.1415);
            lTmp(1, 1) = sin(3.1415);
            //
            lTmp(0, 2) = exp(3.1415);
            lTmp(1, 2) = log(3.1415);
            //
            lTmp(0, 3) = tan(3.1415);
            lTmp(1, 3) = cosh(3.1415);

LzLogM("", "*** lTmp.Buffer()= "<<lTmp.Buffer())

            return lTmp;
        }();
LzLogM("", "***    B.Buffer()= "<<B.Buffer())
#endif
    // Log
    B.Log();

    // Compute Min and Max
    const double lMin = B.GetMin();
    const double lMax = B.GetMax();

    // Log
    LzLogM("", "Min= "<<lMin);
    LzLogM("", "Max= "<<lMax);

    // CHECK
    if( lMin != -3.1415 )
        LzLogException("", "ASSERT FAILED")

    // CHECK
    if( lMax != exp(3.1415) )
        LzLogException("", "ASSERT FAILED")

    // Log
    LzLogM("", "OK")
}

//================================================================================
void Array_UnitTest_5()
{
    // Log
    LzLogN("", "Test 5: transformations")

    // Create source array
    const Array<string, short> A( [](Array<string, short> * pThis)
        {
            pThis->Resize( { 3, 4, 5 } );
            for( auto i=pThis->InitIdx() ; pThis->CheckIdx(i) ; pThis->NextIdx(i) )
                (*pThis)( i ) = "Slice "+pThis->IdxToString( i );

        }, "Source" );

    // Log
    A.Log();

    // Sub array
    Array<string, short> A_x2x("sub-array");
    A.SubArrayTo( {-1, 2, -1}, A_x2x );

    // Log
    A_x2x.Log();

    // CHECK
    if( A_x2x.SizeVec() != vector<size_t>{ 3, 5 } )
        LzLogException("", "ASSERT FAILED")

    // CHECK
    {
        for( size_t c=0 ; c<A_x2x.Count() ; c++ )
        {
            // Expected position in ALL strings= 9
            if( A_x2x.Buffer()[ c ].find(", 2,") != 9 )
                LzLogException("", "ASSERT FAILED")
        }
    }

    // Log
    LzLogM("", "OK")
}
}
