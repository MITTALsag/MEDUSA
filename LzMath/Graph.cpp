#include "Graph.h"
#include <LzTriModel/Mesh.h>
#include <LzTriModel/MeshTopology.h>
////
////*** Unit test
//#include <LzMath/ToolBox.h>
//#include <LzTriModel/MeshTools.h>
//#include <LzGeom/RigidTr3D.h>


namespace LzGraph
{
using LzGeom::Point3D;
using LzTriModel::Mesh;
using LzTriModel::MeshTopology;


//==============================================================================
void Graph::Init( const Mesh & pMesh, const MeshTopology * ppTopo/*=nullptr*/ )
{
    // Set-up detonator
    LzServices::Detonator lDet( [&]{ mNodes.clear(); },
                                [&]{ mNodes.clear(); } );
    // Check
    if( pMesh.mVertices.size()==0 || pMesh.mTriangles.size()==0 )
        LzLogException("", "Cannot initialize graph from an empty mesh! Nb. vertices= "<<pMesh.mVertices.size()<<"; nb. triangles= "<<pMesh.mTriangles.size()<<".")

    // Check
    MeshTopology lTopo;
    const MeshTopology * lpTopo;
    if( ppTopo )
    {
        // Check
        if( !ppTopo->IsCompatible(pMesh) )
            LzLogException("", "Provided topology is not compatible with the provided mesh!")

        // Link
        lpTopo = ppTopo;
    }
    else
    {
        // Use local object
        lTopo.Set( pMesh );
        lpTopo = &lTopo;
    }


    //---------------------------------
    // Init data
    //---------------------------------

    // Get vertices
    const vector<Point3D> & lVers = pMesh.mVertices;

    // Clear previous, and resize graph
    mNodes.resize( lVers.size() );

    // Consider all edges
    for( size_t e=0 ; e<lpTopo->TopEdges().size() ; e++ )
    {
        // Get top edge
        const LzTriModel::TopEdge & lTE = lpTopo->TopEdges( e );

        // Set link in both nodes
        const double lLen = lVers[ lTE.mV[0] ].DistanceTo( lVers[ lTE.mV[1] ] );
        //
        mNodes[ lTE.mV[0] ].mNextVers.AddTail( lTE.mV[1] );
        mNodes[ lTE.mV[1] ].mNextVers.AddTail( lTE.mV[0] );
        //
        mNodes[ lTE.mV[0] ].mNextEds.AddTail( e );
        mNodes[ lTE.mV[1] ].mNextEds.AddTail( e );
        //
        mNodes[ lTE.mV[0] ].mLens.AddTail( lLen );
        mNodes[ lTE.mV[1] ].mLens.AddTail( lLen );
    }


    //---------------------------------
    // Check connectivity
    //---------------------------------

    //
    // Here: all mDist in all mNodes have max value (see Node's construction)
    //
    {
        // Value meaning 'node visited'
        constexpr double DONE_VALUE = 0.0;

        // Stack
        List<size_t> lStack;

        // Seed
        mNodes[0].mDist = DONE_VALUE;
        lStack.AddTail( 0 );

        // Process stack
        while( lStack.Count() )
        {
            // Pop node
            size_t lCurr_n = lStack.GetHead();
            lStack.DelHead();

            // Get node
            const Node & lNode = mNodes[ lCurr_n ];

            // Process neighbors
            BrowseList( iN, lNode.mNextVers )
            {
                // Read index of next node
                const size_t lNext_n = lNode.mNextVers.GetAt( iN );

                // Get next node's dist
                double & lNext_Dist = mNodes[ lNext_n ].mDist;

                // Already processed?
                if( lNext_Dist == DONE_VALUE )
                    continue;

                // Process
                lNext_Dist = DONE_VALUE;
                lStack.AddTail( lNext_n );
            }
        }

        // Check for unreachable nodes
        size_t lNbUnreach = 0;
        for( const Node & iN : mNodes )
        {
            if( iN.mDist != DONE_VALUE )
                lNbUnreach++;
        }

        // Check
        if( lNbUnreach )
            LzLogException("", "Found "<<lNbUnreach<<" node(s) which are unreachable from node 0! Total nb. of nodes= "<<mNodes.size()<<".")
    }
    //
    // Here: we can leave all mDist's set to DONE_VALUE, since they will be reset prior to path finding
    //


    // All is well
    lDet.defuse();
}

//==============================================================================
void Graph::FindShortestPath( size_t pFrom, size_t pTo, List<size_t> & pIdx, PathMode pMode )
{
//#define LOG_STUFF

    // Local type: index + distance
    using DistNode = LzServices::Sortable<size_t, double>;

    // Check
    if( pFrom>=mNodes.size() || pTo>=mNodes.size() )
        LzLogException("", "Incompatible indices! From= "<<pFrom<<", To= "<<pTo<<", Nodes= "<<mNodes.size()<<".")

    // Check
    if( pFrom == pTo )
    {
        // Log
        LzLogM("", "*** NOT AN ERROR: Check for possible unwanted behavior! Computing path from "<<pFrom<<" to itself, "<<pTo<<".")

        // Return an empty list
        pIdx.DelAll();

        return;
    }

    // Set-up detonator
    LzServices::Detonator lDet( [&]{ pIdx.DelAll(); }, [&]{ pIdx.DelAll(); } );

    //----------------------
    // Initialize queue
    //----------------------

    List< DistNode > lQ;
    {
#ifdef LOG_STUFF
        // Log
        LzLogTimeN("", "Init queue")
#endif
        // Seed start node
        mNodes[ pFrom ].mDist = 0.0;
        mNodes[ pFrom ].mPosInQ = lQ.AddTail( DistNode{pFrom, 0.0} );

        // Add all nodes
        for( size_t n=0 ; n<mNodes.size() ; n++ )
        {
            // Skip pFrom node
            if( n == pFrom )
                continue;

            // Read node
            Node & lN = mNodes[n];

            // Reset distance to max value, meaning 'uninitialized'
            lN.mDist = std::numeric_limits<double>::max();

            // Add to tail
            lN.mPosInQ = lQ.AddTail( DistNode{n, lN.mDist} );
        }
    }

    //
    // Here: all mPosInQ and all mDist have been set
    //

    //----------------------
    // Dijkstra loop
    //----------------------

    {
#ifdef LOG_STUFF
        // Log
        LzLogTimeN("", "Dijkstra loop")
#endif
        while( lQ.Count() )
        {
            // Pop head (smallest element)
            DistNode lDN = lQ.GetHead();
            lQ.DelHead();

            // Read node index and distance
            const size_t lCurr_n = lDN.mT;
            const double lCurr_d = lDN.mScore;

            // Get node
            Node & lNode = mNodes[ lCurr_n ];

            // Mark as removed from Q
            lNode.mPosInQ = nullptr;

            // Finished?
            if( lCurr_n == pTo )
            {
#ifdef LOG_STUFF
                // Log
                LzLogMsg("", "Finished! Curr d= "<<lNode.mDist)
#endif
                break;
            }

            // Update all nodes and push to Q
            /*
                for( std::list<int>::const_iterator x=L1.begin(), y=L2.begin(), z=L3.begin() ;
                     x!=L1.end() && y!=L2.end() && z!=L3.end() ;
                     x++, y++, z++ )
                {
                    ...
                }
            */
            BrowseThreeLists( iN, lNode.mNextVers,
                              iE, lNode.mNextEds,
                              iL, lNode.mLens )
            {
                // Read params for next node
                const size_t lNext_n = lNode.mNextVers.GetAt( iN );
                const size_t lNext_e = lNode.mNextEds.GetAt( iE );
                const double lNext_d = lNode.mLens.GetAt( iL );

                // Get next node
                Node & lNext_Node = mNodes[ lNext_n ];

                // Skip node if not in Q anymore
                if( lNext_Node.mPosInQ == nullptr )
                    continue;

                // New candidate distance
                const double lNew_d = lCurr_d + lNext_d;

                // Is better?
                if( lNext_Node.mDist > lNew_d )
                {
                    // Update
                    lNext_Node.mDist = lNew_d;
                    lNext_Node.mPrevVer = lCurr_n;

                    // Find edge to prev, if requested
                    if( pMode == PathMode::EdgeIndices )
                        lNext_Node.mPrevEd = lNext_e;

                    // Replace in Q
                    {
                        // Remove from previous position
                        lQ.DelAt( lNext_Node.mPosInQ );

                        // Reinsert in sorted Q with the updated distance
                        // NB: AddIntoIncrList won't return nullptr since the added element
                        //     is not already present in the list
                        lNext_Node.mPosInQ = lQ.AddIntoIncrList( DistNode{lNext_n, lNew_d}, true );
                    }
                }
            }
        }
    }

    //----------------------
    // Read path
    //----------------------

    for( size_t iCursor=pTo ; iCursor!=pFrom ; iCursor=mNodes[iCursor].mPrevVer )
    {
#ifdef LOG_STUFF
        // Log
        LzLogM("", iCursor<<" -> "<<mNodes[iCursor].mDist)
#endif
        // Stash vertices or edges
        if( pMode == PathMode::EdgeIndices )
            pIdx.AddHead( mNodes[iCursor].mPrevEd );
        else
            pIdx.AddHead( iCursor );
    }

    // Complete path (only in vertices mode)
    if( pMode == PathMode::VertexIndices )
        pIdx.AddHead(pFrom);

#ifdef LOG_STUFF
    // Log
    LzLogN("", "Path ("<<pIdx.Count()<<" nodes)")
    LzLogM("", "From= "<<pFrom)
    LzLogM("", "To=   "<<pTo)
    LzLogM("", "Path= "<<pIdx.ListToString())
#endif

    // All is well
    lDet.defuse();
}

//==============================================================================
void GetPathFromWaypoints( const vector<size_t> & pWaypoints,
                           List<size_t> & pPath,
                           const Mesh & pMesh,
                           const MeshTopology & pTopo )
{
    // Set-up detonator
    LzServices::Detonator lDet( [&] { pPath.DelAll(); },
                                [&] { pPath.DelAll(); } );
    // Init graph
    Graph lGraph;
    lGraph.Init( pMesh, &pTopo );

    // Go through the whole outline
    for( size_t p=0 ; p+1<pWaypoints.size() ; p++ )
    {
        // Current and next stop
        const size_t lFrom = pWaypoints[p + 0];
        const size_t lTo   = pWaypoints[p + 1];

        // Not moving?
        if( lTo == lFrom )
            continue;

        // Next of kin? (faster than path extraction)
        {
            // Get topo edges
            const List<size_t> & lTopEds = pTopo.TopVers( lFrom ).mE;

            // Check opposite vertices
            BrowseList( iE, lTopEds )
            {
                // Get opposite vertex index
                const size_t lOppVerIdx = pTopo.TopEdges( lTopEds.GetAt(iE) ).OtherVer( lFrom );

                // Is it lTo?
                if( lOppVerIdx == lTo )
                {
                    // Short circuit
                    pPath.AddTail( lTopEds.GetAt(iE) );
                    goto SkipGraphSearch;
                }
            }
        }

        // Nope, not next of kin. We need to look for a path.
        {
            // Find path
            List<size_t> lPath;
            lGraph.FindShortestPath( lFrom, lTo, lPath, PathMode::EdgeIndices );

            // Stash
            pPath.Append( std::move(lPath) );
        }

    SkipGraphSearch:
            ;
    }

    // All is well
    lDet.defuse();
}
}
