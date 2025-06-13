#include "Graph_UnitTest.h"
#include "Graph.h"
#include <LzTriModel/Mesh.h>
#include <LzTriModel/MeshTools.h>
#include <LzMath/ToolBox.h>
#include <LzGeom/RigidTr3D.h>


namespace LzGraph
{
using LzGeom::Point3D;
using LzGeom::Plane3D;


//==============================================================================
void Graph_UnitTest_ALL()
{
    Graph_UnitTest_1();
    Graph_UnitTest_2();
}

//==============================================================================
void Graph_UnitTest_1()
{
    // Log
    LzLogTimeN("", "Unit testing LzGraph::Graph")

    /*
     * Test cases
     * ==========
     *  - Initialize with multipatch
     *  - Compute distance using both vertex and edge indices
     */

    // Create sphere
    Mesh lSphere;
    LzTriModel::GenEllipsoid( lSphere, Point3D(), 5.0, 5.0, 5.0, 120, 60 );


    //--------------------------
    // Bad init
    //--------------------------
    {
        // Do not move lSphere
        Mesh lMultiMesh = lSphere;
        //
        lMultiMesh.RigidTransform( {0,0,0, {20.0, 0, 0}} );
        lMultiMesh.Merge( lSphere );
#if 0
        // Save
        lMultiMesh.Save(LzServices::StartUpPath()+"/_unit_test_Graph_sphere.obj");
#endif
        // Set graph
        Graph lGr;
        try
        {
            // Should fail
            lGr.Init( lMultiMesh );
        }
        catch( ... )
        {
            goto NoError;
        }

        // Error
        LzLogException("", "Init DID NOT FAIL!")

    NoError:
        // Ok
        LzLogM("", "Init failed, as it should.")
    }


    //--------------------------
    // Path finding
    //--------------------------
    {
        // Set topo
        MeshTopology lTopo;
        lTopo.Set( lSphere );

        // Set graph
        Graph lGr;
        lGr.Init( lSphere );

        // From North pole to South pole (on first sphere)
        const size_t lFrom = 0;
        const size_t lTo   = 7141;

        // Vertices
        List<size_t> lVerIdx;
        lGr.FindShortestPath( lFrom, lTo, lVerIdx, LzGraph::PathMode::VertexIndices );

        // Compute distance
        double lVerDist = 0;
        for( void * iV=lVerIdx.HeadPos() ; iV!=lVerIdx.TailPos() ; )
        {
            // A and B
            const Point3D & A = lSphere.mVertices[ lVerIdx.GetAtAndNext(iV) ];
            const Point3D & B = lSphere.mVertices[ lVerIdx.GetAt(iV) ];

            // Update distance
            lVerDist += A.DistanceTo( B );
        }

        // Log
        LzLogM("", "Distance from vertices= "<<lVerDist)

        // Edges
        List<size_t> lEdgIdx;
        lGr.FindShortestPath( lFrom, lTo, lEdgIdx, LzGraph::PathMode::EdgeIndices );

        // Compute distance
        double lEdgDist = 0;
        BrowseList( iE, lEdgIdx )
        {
            // Read edge
            const LzTriModel::TopEdge & lTEd = lTopo.TopEdges( lEdgIdx.GetAt(iE) );

            // A and B
            const Point3D & A = lSphere.mVertices[ lTEd.mV[0] ];
            const Point3D & B = lSphere.mVertices[ lTEd.mV[1] ];

            // Update distance
            lEdgDist += A.DistanceTo( B );
        }

        // Log
        LzLogM("", "Distance from edges= "<<lEdgDist)

        // Compute ideal distance
        const double lDiam = lSphere.mVertices[lFrom].DistanceTo( lSphere.mVertices[lTo] );
        const double lRefDist = 0.5 * LzMath::PI * lDiam;

        // Log
        LzLogM("", "Ideal distance= "<<lRefDist)

        // Check
        if( std::abs(lVerDist - lEdgDist) > 1e-6 )
            LzLogException("", "Found different distances, between ver and edge! "<<lVerDist<<" != "<<lEdgDist<<".")

        // Check
        if( std::abs(lVerDist - lRefDist) > 1e-2 )
            LzLogException("", "Found different distances, between ver and ref! "<<lVerDist<<" != "<<lRefDist<<".")
    }
}

//==============================================================================
void Graph_UnitTest_2()
{
    // Log
    LzLogTimeN("", "Unit testing LzGraph::Graph")

    // Create test mesh
    Mesh lMesh;
    LzTriModel::GenEllipsoid( lMesh, Point3D(), 50.0, 60.0, 80.0, 30, 40 );

    // Cut plane
    const Plane3D lCut(0.74, 0.21, 0.63, 5.85);

    // Cut mesh
    lMesh.Cut( lCut, 1e-6 );

    // Clean
    while( !lMesh.RemoveDuplicateVerticesAndNormals() );
        lMesh.RemoveUnusedVerticesAndNormals();

//lMesh.Save(LzServices::StartUpPath()+"/_tmp_ELL.obj");

    // Compute topology
    MeshTopology lTopo;
    lTopo.Set( lMesh );

    // Init graph
    Graph lG;
    lG.Init( lMesh, &lTopo );

    // Extract all possible EDGES and VERTICES and compare
//const size_t v = 246;
    for( size_t v=0 ; v<lMesh.mVertices.size() ; v++ )
    {
        // Log
        LzLogM("", "v= "<<v)

//const size_t w = 523;
        for( size_t w=0 ; w<lMesh.mVertices.size() ; w++ )
        {
//// Log
//LzLogN("", "v= "<<v<<", w= "<<w)

            // Lists of indices
            List<size_t> lVerIdx;
            List<size_t> lEdsIdx;

            // Get indices
            lG.FindShortestPath( v, w, lVerIdx, LzGraph::PathMode::VertexIndices );
            lG.FindShortestPath( v, w, lEdsIdx, LzGraph::PathMode::EdgeIndices );

            // Count
            const size_t lNbVer = lVerIdx.Count();
            const size_t lNbEds = lEdsIdx.Count();

//// Log
//LzLogM("", "Computed indices. "<<lNbVer<<" vertices, "<<lNbEds<<" edges.")

            // Same?
            if( v == w )
            {
                //-------------
                // Empty list
                //-------------

                if( lNbVer || lNbEds )
                    LzLogException("", "Count error: Vers= "<<lNbVer<<" OR "<<lNbEds<<", != 0!")
            }
            else
            {
                //-------------
                // Compare
                //-------------

                // Check length
                if( lNbVer != lNbEds+1 )
                    LzLogException("", "Count error: Vers= "<<lNbVer<<" vs "<<lNbEds<<" + 1!")

                // Extract vertices from edges
                List<size_t> lEdsVer;
                lEdsVer.AddTail( v );
                BrowseList( iE, lEdsIdx )
                {
                    const size_t e = lEdsIdx.GetAt( iE );
                    const LzTriModel::TopEdge & lTEd = lTopo.TopEdges( e );

                    const size_t last_v = lEdsVer.GetTail();
                    const size_t next_v = lTEd.OtherVer( last_v );

                    lEdsVer.AddTail( next_v );
                }

                // Must be same
                if( lVerIdx != lEdsVer )
                {
                    LzLogE("", "V= "<<lVerIdx.ListToString())
                    LzLogE("", "E= "<<lEdsVer.ListToString())
                    LzLogException("", "Vertex indices lists differ!")
                }
                else
                {
//LzLogM("", "V= "<<lVerIdx.ListToString())
//LzLogM("", "E= "<<lEdsVer.ListToString())
                }
            }
        }
    }

    // Log
    LzLogM("", "--- OK")
}
}
