#include "List_UnitTest.h"
#include "List.h"


namespace LzServices
{
//================================================================================
void List_UnitTest_ALL()
{
    List_UnitTest_1();
    List_UnitTest_2();
    List_UnitTest_3();
}

//================================================================================
void List_UnitTest_1()
{
    //--------------------
    // Test: Sortable
    //--------------------

    // int
    {
        LzLogN("", "----- Testing with 'int'")

        using SrInt = Sortable<int, double>;

        const auto NullIfNoAdd    = LzServices::List<SrInt>::AddIntoReturnPolicy::NullIfNoAdd;
        const auto PosOfIdentical = LzServices::List<SrInt>::AddIntoReturnPolicy::PosOfIdentical;

        const SrInt I1 (1, 1.0);
        const SrInt I2a(2, 2.1);
        const SrInt I2b(2, 2.2);
        const SrInt I3 (3, 3.0);
        const SrInt I4a(4, 3.0);
        const SrInt I4b(4, 3.0);

        // Unique = true
        {
            List<SrInt> lList;
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I1,  true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I2a, true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I2b, true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I3,  true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I4a, true, NullIfNoAdd ))
            LzLogM("", "Should be     0: "<<lList.AddIntoIncrList( I4b, true, NullIfNoAdd ))
            LzLogM("", "")
        }

        // Unique = true
        {
            List<SrInt> lList;
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I1,  true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I2a, true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I2b, true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I3,  true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I4a, true, PosOfIdentical ))
            LzLogM("", "Should be same as above: "<<lList.AddIntoIncrList( I4b, true, PosOfIdentical ))
        }
    }

    // Point3D
    {
        LzLogN("", "----- Testing with 'Point3D'")

        using LzGeom::Point3D;
        using SrPnt = Sortable<Point3D, double>;

        const auto NullIfNoAdd    = LzServices::List<SrPnt>::AddIntoReturnPolicy::NullIfNoAdd;
        const auto PosOfIdentical = LzServices::List<SrPnt>::AddIntoReturnPolicy::PosOfIdentical;

        const Point3D A(1.0/3.0, 0.25, 0.36);
        const Point3D B(0.76, 1.0/7.0, 3.1415926);
        const Point3D C(11, 22, 33);

        const SrPnt I1 (A, 1.0);
        const SrPnt I2a(A, 2.1);
        const SrPnt I2b(A, 2.2);
        const SrPnt I3 (B, 3.0);
        const SrPnt I4a(C, 3.0);
        const SrPnt I4b(C, 3.0);

        // NullIfNoAdd
        {
            List<SrPnt> lList;
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I1,  true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I2a, true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I2b, true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I3,  true, NullIfNoAdd ))
            LzLogM("", "Should be NON 0: "<<lList.AddIntoIncrList( I4a, true, NullIfNoAdd ))
            LzLogM("", "Should be     0: "<<lList.AddIntoIncrList( I4b, true, NullIfNoAdd ))
            LzLogM("", "")
        }

        // PosOfIdentical
        {
            List<SrPnt> lList;
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I1,  true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I2a, true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I2b, true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I3,  true, PosOfIdentical ))
            LzLogM("", "Should be NON 0        : "<<lList.AddIntoIncrList( I4a, true, PosOfIdentical ))
            LzLogM("", "Should be same as above: "<<lList.AddIntoIncrList( I4b, true, PosOfIdentical ))
        }
    }
}

//================================================================================
// Google test bla bla bla
#define EXPECT_TRUE( COND ) \
    if( !(COND) ) \
        LzLogException("", "CHECK FAILED!")

//================================================================================
void List_UnitTest_2()
{
    // Test list
    List<double> L1;
    L1.AddIntoIncrList( 13, true );
    L1.AddIntoIncrList( 3, true );
    L1.AddIntoIncrList( 2.2, true );
    L1.AddIntoIncrList( 33, true );
    L1.AddIntoIncrList( 0.3, true );
    L1.AddIntoIncrList( -33, true );

    // Log
    LzLogM("", "L1= "<<L1.ListToString())

    // CHECK
    EXPECT_TRUE( L1.Count() == 6 );
    EXPECT_TRUE( L1.GetHead() == -33 );
    EXPECT_TRUE( L1.GetTail() == 33 );


    // Test list
    List<double> L2;
    L2.AddIntoIncrList( -22, true );
    L2.AddIntoIncrList( 2.2, true );
    L2.AddIntoIncrList( 2.2, /*NON UNIQUE*/false );
    L2.AddIntoIncrList( 222, true );
    L2.AddIntoIncrList( -0.2, true );
    L2.AddIntoIncrList( -2.22, true );

    // Log
    LzLogM("", "L2= "<<L2.ListToString())

    // CHECK
    EXPECT_TRUE( L2.Count() == 6 );
    EXPECT_TRUE( L2.GetHead() == -22 );
    EXPECT_TRUE( L2.GetTail() == 222 );


    // Test list
    List<double> L3 = L1;
    L3.MergeWithIncrList( L2, /*pUniqueFromBtoA=*/false );

    // Log
    LzLogM("", "L3= "<<L3.ListToString())

    // CHECK
    EXPECT_TRUE( L3.Count() == L1.Count() + L2.Count() );


    // Test list
    List<double> L4 = L1;
    L4.MergeWithIncrList( L2, /*pUniqueFromBtoA=*/true );

    // Log
    LzLogM("", "L4= "<<L4.ListToString())

    // CHECK
    EXPECT_TRUE( L4.Count() == 10 );
}

//================================================================================
void List_UnitTest_3()
{
    List<int> L1, L2;

    // Test empty lists
    {
        List<int> L1, L2({1, 5, 12});
        L1.MergeUniqueIncrLists( L2 );
        LzLogM("", "Empty+sorted= "<<L1);
    }
    {
        List<int> L1({1, 5, 12}), L2;
        L1.MergeUniqueIncrLists( L2 );
        LzLogM("", "Sorted+empty= "<<L1);
    }
    {
        List<int> L1({1}), L2;
        L1.MergeUniqueIncrLists( L2 );
        LzLogM("", "Single+empty= "<<L1);
    }
    {
        List<int> L1, L2({1});
        L1.MergeUniqueIncrLists( L2 );
        LzLogM("", "Empty+single= "<<L1);
    }
    {
        List<int> L1, L2;
        L1.MergeUniqueIncrLists( L2 );
        LzLogM("", "Empty+empty= "<<L1);
    }

    LzLogM("", "");

    // Test long lists
#define LEN 50
    std::uniform_int_distribution<int> lDist(-20,+20);

    for( int test=0 ; test<10 ; test++ )
    {
        List<int> L1;
        for( int i=0 ; i<LEN ; i++ )
            L1.AddIntoIncrList( lDist(LzServices::RandomEngine()), true );

        List<int> L2;
        for( int i=0 ; i<LEN ; i++ )
            L2.AddIntoIncrList( lDist(LzServices::RandomEngine()), true );

        LzLogM("", "L1= "<<L1);
        LzLogM("", "L2= "<<L1);

        // Good answer
        List<int> L3 = L1;
        BrowseList( i, L2 )
            L3.AddIntoIncrList( L2.GetAt(i), true );

        // Check
        L1.MergeUniqueIncrLists( L2 );
        if( L1 != L3 )
            LzLogM("", "*** ERROR !!!")
        else
            LzLogM("", "...ok")
    }
#undef LEN
}




}

























