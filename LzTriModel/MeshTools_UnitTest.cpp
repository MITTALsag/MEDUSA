#include "MeshTools_UnitTest.h"
#include "MeshTools.h"
#include <LzMath/Graph.h>


namespace LzTriModel
{
using LzGraph::Graph;


//==============================================================================
void MeshTools_UnitTest_ALL()
{
    MeshTools_UnitTest_1();
    MeshTools_UnitTest_2();
}

//==============================================================================
void MeshTools_UnitTest_1()
{
    // Log
    LzLogN("", "MeshTools_UnitTest_1")

    // Gen ell
    Mesh lEll;
    GenEllipsoid( lEll, Point3D(0,0,0), 10, 20, 30, 50, 100 );
    //
    lEll.Save(LzServices::StartUpPath()+"/_unit_test_MeshTools_Ellipsoid.obj");

    // Topo
    MeshTopology lTopo;
    lTopo.Set( lEll );

    for( int test=0 ; test<10 ; test++ )
    {
        // Log
        LzLogN("", "MeshTools_UnitTest_1: round "<<test)

        // Pick a random edge
        std::uniform_int_distribution<size_t> lRnd(0, lTopo.TopEdges().size()-1/*NB: '-1' because values in [a, b]*/);
        size_t lEdIdx = lRnd( LzServices::RandomEngine() );

        // The edge
        const TopEdge & lEd = lTopo.TopEdges( lEdIdx );
        const Point3D & lEdVer0 = lEll.mVertices[lEd.mV[0]];
        const Point3D & lEdVer1 = lEll.mVertices[lEd.mV[1]];
        const Point3D lEdPt = lEdVer0.MidPointTo( lEdVer1 );

        // Find opposite ver
        Point3D lOppVer;
        {
            double lMaxDist = 0.0;
            for( size_t v=0 ; v<lEll.mVertices.size() ; v++ )
            {
                const Point3D & iVer = lEll.mVertices[v];
                double lDist = iVer.DistanceTo( lEdPt );
                if( lMaxDist < lDist )
                {
                    lMaxDist = lDist;
                    lOppVer = iVer;
                }
            }
        }

        // Build cut plane
        const Plane3D lCut( lEdVer0, lEdVer1, lOppVer );

        // Find cut index
        size_t lFromIdx = 0; // To silence warnings about uninitialized values
        for( size_t v=0 ; v<lEll.mVertices.size() ; v++ )
        {
            double lAbsDst = std::abs( lCut.SignedDistanceTo( lEll.mVertices[v] ) );

            // Start vertex must lie outside the cut plane
            if( lAbsDst > 1.0 )
            {
                lFromIdx = v;
                break;
            }
        }

        // Log
        LzLogM("", "Cut plane= "<<lCut.ToString_ABCD())

        // Cut from vertex
        Mesh lVerMesh, lOppMesh;
        CutFromVertex( lEll, lTopo, lFromIdx, lCut, lVerMesh, lOppMesh );

        // Save
        lVerMesh.Save(LzServices::StartUpPath()+"/_unit_test_MeshTools_VerMesh.obj");
        lOppMesh.Save(LzServices::StartUpPath()+"/_unit_test_MeshTools_OppMesh.obj");

        // Check non empty
        if( lVerMesh.mTriangles.size()==0 )
            LzLogException("", "Empty ver mesh!")

        // Check non empty
        if( lOppMesh.mTriangles.size()==0 )
            LzLogException("", "Empty opp mesh!")

        // Merge again
        Mesh lMergeVerOpp = lVerMesh;
        lMergeVerOpp.Merge( lOppMesh );
        while( !lMergeVerOpp.RemoveDuplicateVerticesAndNormals() );

        // Compare area
        const double lRefArea = lEll.Area();
        const double lTestArea = lMergeVerOpp.Area();

        // Log
        LzLogM("", "Ref area= "<<lRefArea)
        LzLogM("", "Test area= "<<lTestArea)

        // Check
        if( std::abs(lRefArea - lTestArea) > 1e-5 )
            LzLogException("", "Area mismatch!")

        // Check topo
        MeshTopology lMergeTopo;
        lMergeTopo.Set( lMergeVerOpp );

        // Check
        if( lMergeTopo.Cracks().size() )
            LzLogException("", "Found cracks!")

        // Check
        if( lMergeTopo.MultiEdges().size() )
            LzLogException("", "Found multi edges!")
    }


}

//==============================================================================
void MeshTools_UnitTest_2()
{
    // Log
    LzLogN("", "MeshTools_UnitTest_2")

    // Gen ell
    Mesh lEll;
    GenEllipsoid( lEll, Point3D(0,0,0), 10, 20, 30, 50, 100 );
    //
    lEll.Save(LzServices::StartUpPath()+"/_unit_test_MeshTools_Test_2_Ellipsoid.obj");

    // Indices for outline computation
    const vector<size_t> lOutIdx{ 606, 646, 2347, 3655, 4501, 2204, 606 };
    const size_t lFrom_0 = 1676;
    const size_t lFrom_1 = 3531;

    // Topo
    MeshTopology lTopo;
    lTopo.Set( lEll );

    // Compute hermetic edges
    List<size_t> lHermEdIdx;
    {
        // Create graph
        Graph lGr;
        lGr.Init( lEll, &lTopo );

        for( size_t i0=0 ; i0+1<lOutIdx.size() ; i0++ )
        {
            const size_t lIdx0 = lOutIdx[i0 + 0];
            const size_t lIdx1 = lOutIdx[i0 + 1];

            List<size_t> lEdIdx;
            lGr.FindShortestPath( lIdx0, lIdx1, lEdIdx, LzGraph::PathMode::EdgeIndices );

            // Explicitly invoke move semantics (MSVC2019 won't optimize this call)
            lHermEdIdx.Append( std::move(lEdIdx) );
        }
    }

    // Extract patch 0
    Mesh lPatch_0;
    GetPatchFrom( lEll, lTopo, lFrom_0, lHermEdIdx, lPatch_0 );
    //
    lPatch_0.Save(LzServices::StartUpPath()+"/_unit_test_MeshTools_Test_2_Patch_0.obj");

    // Extract patch 1
    Mesh lPatch_1;
    GetPatchFrom( lEll, lTopo, lFrom_1, lHermEdIdx, lPatch_1 );
    //
    lPatch_1.Save(LzServices::StartUpPath()+"/_unit_test_MeshTools_Test_2_Patch_1.obj");

    // Check
    const double lRefArea = lEll.Area();
    const double lArea_0 = lPatch_0.Area();
    const double lArea_1 = lPatch_1.Area();
    //
    if( std::abs(lRefArea - lArea_0 - lArea_1) > 1e-6 )
        LzLogException("", "Assertion failed!")

    // Merge
    Mesh lMergeEll = lPatch_0;
    lMergeEll.Merge( lPatch_1 );
    //
    while( !lMergeEll.RemoveDuplicateVerticesAndNormals() );

    // Check topo
    MeshTopology lMergeTopo;
    lMergeTopo.Set( lMergeEll );

    // Check
    if( lMergeTopo.Cracks().size() || lMergeTopo.MultiEdges().size() )
        LzLogException("", "Assertion failed!")

    // Check
    const double lMergeArea = lMergeEll.Area();
    if( std::abs(lMergeArea - lRefArea) > 1e-6 )
        LzLogException("", "Assertion failed!")
}

}
