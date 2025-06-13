#include "MeshTools.h"
#include "Outliner.h"
#include "PatchComputer.h"
#include <LzMath/ToolBox.h>
#include <LzGeom/Line3D.h>
#include <map>


namespace LzTriModel
{
using std::map;
using std::tuple;



////================================================================================
//void ContractMesh( Mesh & pMesh, const List<size_t> & pMobileVers )
//{
//    MeshTopology lTopo;
//    lTopo.Set(pMesh);

//    vector<Vector3D> lForces;
//    lForces.assign( pMobileVers.Count(), Vector3D{0, 0, 0} );

//    vector<double> lLen0( pMobileVers.Count() );

//    const double contract_alpha = 0.8;
//    const double step_alpha = 0.2;
//    const double K = 1.0;
//    for( int iter=0 ; iter<1000 ; iter++ )
//    {
//        // Compute forces
//        size_t lIdxF = 0;
//        BrowseList( iV, pMobileVers )
//        {
//            const size_t lV = pMobileVers.GetAt(iV);
//            const TopVer & lTV = lTopo.TopVers(lV);


//            // Add contributions from all edges
//            BrowseList( iE, lTV.mE )
//            {
//                const TopEdge & lTE = lTopo.TopEdges( lTV.mE.GetAt(iE) );

//                const size_t lOtherV = lTE.OtherVer( lV );


//                Vector3D U = pMesh.mVertices[lOtherV] - pMesh.mVertices[lV];
//                U.Normalize();

//                const double lCurrLen = pMesh.mVertices[lOtherV].DistanceTo( pMesh.mVertices[lV] );
//                const double lNewLen = contract_alpha * lCurrLen;

//                if( iter == 0 )
//                    lLen0[lIdxF] = lCurrLen;

//                const Vector3D lForce = K*(lCurrLen/*lLen0[lIdxF]*/ - lNewLen)*U;

//                lForces[lIdxF] += lForce;
//            }

//            // Next force index
//            lIdxF++;
//        }

//        // Step down forces direction
//        lIdxF = 0;
//        BrowseList( iV, pMobileVers )
//        {
//            const size_t lV = pMobileVers.GetAt(iV);
//            pMesh.mVertices[lV] += step_alpha * lForces[lIdxF++];
//        }



//if( iter % 20 == 0 )
//    pMesh.Save(LzServices::StartUpPath()+"/_debug__CONTRACT_MESH__"+std::to_string(iter)+".obj");




//        // Reset forces
//        std::fill(lForces.begin(), lForces.end(), Vector3D{0, 0, 0});
//    }
//}


#pragma region "Patch extraction"
//================================================================================
size_t GetTriangleAwayFromCracks( const Mesh & pMesh, const MeshTopology & pTopo )
{
    // Check if there are triangles in the mesh
    if( pMesh.mTriangles.size() == 0 )
        LzLogException("", "Cannot find triangle most away from cracks in a mesh without triangles! You idiot.")

    // Check for cracks
    if( pTopo.Cracks().size() == 0 )
        LzLogException("", "Cannot find triangle most away from cracks in a mesh without cracks! Dumbass.")

    // Triangles table
    vector<int> lTriDist;
    lTriDist.assign( pTopo.TopTris().size(), -1 ); // -1 means 'uninitialized'

    // Init first round from cracks
    List<size_t> lTriFront;
    for( size_t krk=0 ; krk<pTopo.Cracks().size() ; krk++ )
    {
        // Get crack edge
        const TopEdge & lCrack = pTopo.TopEdges( pTopo.Cracks()[krk] );

        // Crack has a single triangle. Set its distance to cracks to 0
        lTriDist[ lCrack.mT.GetHead() ] = 0;

        // Update advancing tri front
        lTriFront.AddTail( lCrack.mT.GetHead() );
    }

    // Advance front
    int lCurrDist = 0;
    while( lTriFront.Count() )
    {
        // Create new front
        List<size_t> lNewFront;

        // Browse current front
        BrowseList( iIdx, lTriFront )
        {
            // Get triangle index in current front
            const size_t lTriIdx = lTriFront.GetAt( iIdx );

            // Get top tri
            const TopTri & lTTri = pTopo.TopTris( lTriIdx );

            // Check all edges
            for( int e=0 ; e<3 ; e++ )
            {
                // Get top edge
                const TopEdge & lTEd = pTopo.TopEdges( lTTri.mE[e] );

                // Check all triangles in edge
                BrowseList( iT, lTEd.mT )
                {
                    // Get triangle index
                    const int lEdTriIdx = lTEd.mT.GetAt( iT );

                    // Is this triangle unreached?
                    if( lTriDist[ lEdTriIdx ] == -1 )
                    {
                        // Yes, reach it
                        lTriDist[ lEdTriIdx ] = lCurrDist + 1;

                        // Store in new front
                        lNewFront.AddTail( lEdTriIdx );
                    }
                }
            }
        }

        // Move new front to front
        lTriFront = std::move( lNewFront );

        // Increase dist
        lCurrDist++;
    }

    // Find the triangle with the greatest distance
    const size_t lMaxTriIdx = distance(lTriDist.begin(), std::max_element(lTriDist.begin(), lTriDist.end()));
    return lMaxTriIdx;
}

//================================================================================
void GetPatchFrom( const Mesh & pMesh,
                   const MeshTopology & pTopo,
                   size_t pFromVerIdx,
                   const List<size_t> & pHermeticEdgIdx,
                   Mesh & pVerMesh )
{
    // Check
    if( &pMesh == &pVerMesh )
        LzLogException("", "Cannot extract a patch from a mesh into the same mesh!")

    // Check
    if( pFromVerIdx >= pMesh.mVertices.size() )
        LzLogException("", "Invalid vertex index "<<pFromVerIdx<<"! Max= "<<pMesh.mVertices.size())

    // Check
    if( !pTopo.IsCompatible(pMesh) )
        LzLogException("", "Provided topology is not compatible with the provided mesh!")

//    // Check
//    if( pHermeticEdgIdx.Count() == 0 )
//        LzLogException("", "No hermetic edges! The mesh would remain unchanged.")

    //
    // Here: we have at least a non-empty mesh since vertices.size > pFromVerIdx and pFromVerIdx >= 0
    //

    // Set-up the table of hermetic edges
    vector<bool> lHermeticEdge;
    lHermeticEdge.assign( pTopo.TopEdges().size(), false );

//#define CHECK_EDGE_TRI_COUNT

#ifdef CHECK_EDGE_TRI_COUNT
    // Set-up "Hermetic edge has tri" map
    // None of the hermetic edges has a tri
    map<size_t, tuple<bool,size_t>> lHermeticEdgeHasTri;
#endif

    BrowseList( iE, pHermeticEdgIdx )
    {
        // Edge index
        const size_t lEdIdx = pHermeticEdgIdx.GetAt(iE);

        // Check index
        if( lEdIdx >= lHermeticEdge.size() )
            LzLogException("", "Edge index "<<lEdIdx<<" exceeds max edge index in mesh! (size= "<<lHermeticEdge.size()<<")")

        // Hermetic edges table
        lHermeticEdge[ lEdIdx ] = true;

        /**
         * It is impossible to 'simply' check that each hermetic edge has exactly ONE triangle
         * in the extracted mesh since there are configurations where isoleted edge "spurs" can
         * appear, which will be "absorbed" by the extracted mesh. Such edges end up with both
         * triangles in the extracted mesh and it is not a problem at all.
         *
         * Example: "T shaped" edges contour: (A) --(e1)--|--(e3)-- (C)
         * ==================================             |
         *                                               (e2)
         *                                                |
         *                                               (B)
         * Contour vertices: (A) -> (B) -> (C)
         * Contour edges: (e1), (e2 (going down)), (e2 (going up)), (e3)
         *
         * After extraction from a vertex located below (B), the edge (e2) will have 2 triangles.
         *
         */

#ifdef CHECK_EDGE_TRI_COUNT
        // Hermetic edges flag and tri index
        lHermeticEdgeHasTri.emplace( lEdIdx, std::make_tuple(false, 0) );
#endif
    }

    // Setup detonator
    auto CLEAN_UP = [&] { pVerMesh.Free(); };
    LzServices::Detonator lDet( CLEAN_UP, CLEAN_UP );

    // Seed with init vertices
    pVerMesh.mVertices = pMesh.mVertices;

    // Seed with init normals
    pVerMesh.mNormals = pMesh.mNormals;
    //
    // This is ok since we will be pushing whole triangles
    // from the source mesh, whose triangles are assumed to
    // be well formed
    //

    // Flags for processed triangles
    vector<bool> lTriPushed;
    lTriPushed.assign( pMesh.mTriangles.size(), false );

    //---------------------------------
    // Stack
    //---------------------------------

    // Stack for processing triangles
    List<size_t> lTriStack;

    // Seed stack with starting vertex
    const TopVer & lTVer = pTopo.TopVers( pFromVerIdx );
    BrowseList( iT, lTVer.mT )
    {
        // Triangle index
        const size_t lTriIdx = lTVer.mT.GetAt(iT);

        // Push
        lTriStack.AddTail( lTriIdx );

        // Mark as processed
        lTriPushed[ lTriIdx ] = true;
    }


    //---------------------------------
    // Process tri stack
    //---------------------------------

    while( lTriStack.Count() )
    {
        // Pop a triangle from stack
        const size_t t = lTriStack.GetHead();
        lTriStack.DelHead();

        // Get triangle info
        const TopTri & lTTri = pTopo.TopTris( t );

//LzLogN("", "Triangle "<<t)

#ifdef CHECK_EDGE_TRI_COUNT
        // Check that any hermetic edge from this triangle
        // does not ALREADY have a triangle
        for( int e=0 ; e<3 ; e++ )
        {
            // Read edge index
            const size_t lEdIdx = lTTri.mE[ e ];

            // Try to find it within all hermetic edges
            map<size_t, tuple<bool,size_t>>::iterator lPos = lHermeticEdgeHasTri.find( lEdIdx );

            // Hermetic edge?
            if( lPos != lHermeticEdgeHasTri.end() )
            {
//LzLogM("", "---- found HE: "<<lEdIdx)


                // Already has a triangle?
                if( std::get<0>(lPos->second) )
                {
                    const size_t V0 = pTopo.TopEdges( lEdIdx ).mV[0];
                    const size_t V1 = pTopo.TopEdges( lEdIdx ).mV[1];

                    LzLogException("", "Too many triangles sharing hermetic edge "<<lEdIdx<<"(Vers= ["<<V0<<", "<<V1<<"])! "
                                       "Cannot add tri. "<<t<<"; hermetic edge already has tri. "<<std::get<1>(lPos->second)<<".")
                }

                // Ok. Store candidate triangle
                std::get<0>(lPos->second) = true;
                std::get<1>(lPos->second) = t;
            }

//***********************************************************************************
//**** ajouter a la doc de la methode !!!
//    for each hermetic edge, at most 1 triangle !!!
//
// multiedges NOT tolerated in the extracted patch, but ME tolerated in the input mesh !
//
//    exception if not proper mesh
//***********************************************************************************

        }
#endif
        // Stash triangle
        pVerMesh.mTriangles.push_back( pMesh.mTriangles[t] );

        // Check all edges
        for( int e=0 ; e<3 ; e++ )
        {
            // Edge index
            const size_t lEdgIdx = lTTri.mE[e];

            // Skip if hermetic edge
            if( lHermeticEdge[ lEdgIdx ] )
                continue;

            // Skip if CRACK or MULTIEDGE
            if( pTopo.TopEdges( lEdgIdx ).mT.Count() != 2 )
                continue;

            // Get opposite triangle
            const int lOppTri = pTopo.TopEdges( lEdgIdx ).OtherTri( t );

            // Skip if pushed
            if( lTriPushed[lOppTri] )
                continue;

            // Push
            lTriStack.AddTail( lOppTri );
            lTriPushed[ lOppTri ] = true;
        }
    }

#ifndef CHECK_EDGE_TRI_COUNT
    // Final check
    if( pVerMesh.mTriangles.size() == pMesh.mTriangles.size() )
        LzLogException("", "Extracted mesh and source mesh have the same number of triangles! The hermetic edges didn't do a good job splitting the mesh...")
#endif

    // All is well
    lDet.defuse();
}
#pragma endregion


#pragma region "Cutting"
//================================================================================
/**
 * Ajouter OPTIONS pour pouvoir passer FDC ou bien TOPO ... ou rien du tout et tout recalc en local
 *
 * a voir si algo marche avec mode "single" ... a priori oui,
 * il devrait (collection de patches)
 *
 * ajouter pSnapRange !!!
 *
 ******************** manage case : avec ou sans normales !!!
 ******************** manage case : avec ou sans normales !!!
 *
 ******************** plutot que renormalise, calculer les normales aux cut points (prendre les normales des triangles car bord du cut
 ******************** et on ne peut pas formuler l'hypothese de 1 ver = 1 nor pour interpoler les normales....
 * Snapper les vertices ou pas ???? car c'est genant de devoir modifier le mesh ...
 * ou alors copier le mesh localement !
 *
 * opp mesh et ver mesh optionnels car on ne veut pas necessairement les deux
 *
 */
void CutFromVertex( const Mesh & pMesh,
                    const MeshTopology & pTopo,
                    size_t pFromVerIdx,
                    const Plane3D & pCutPlane,
                    Mesh & pVerMesh,
                    Mesh & pOppMesh )
{
    /**
     *
     * This method extract an "edge connected" tri patch containing the starting vertex
     *
     * Each vertex has 3 possible states w.r.t. the cut plane: In, Out or Cut
     * In: the vertex is on the same side of the cut plane as the seed vertex
     * Out: on the opposite side
     * Cut: on cut
     *
     * A tolerance is used to define the threshold between the states. A vertex is Cut
     * iff its signed distance to the plane is within [-tolerance, +tolerance].
     *
     * The opposite mesh can be optionally generated.
     *
     * Edge-connectivity
     * -----------------
     *
     * Two triangles will belong to the same patch (Opp or Ver) iif it is possible
     * to cross from one triangle to the other through a point of the edge that
     * is strictly contained within the edge (i.e., excluding edge extremities),
     * while remaining on the same side of the cut plane.
     *
     * Consequence: for each triangle in a patch (Opp or Ver), each edge is assessed
     * and the edge-opposite triangles are pushed only in the following edge cases:
     * - I-I: all points on the edge are I, the triangle opposite this edge is pushed
     * - C-I: the center of the edge is I
     * - O-I: there is a strict cut on this edge and thus a "crossing" I point
     *
     * The following edge configurations do not result in a pushed opposite triangle:
     * - O-O: all points on the edge are O
     * - C-O: all points on the edge are O, and one extremity is C
     * - C-C: all points on the edge are C
     *
     **/

    //---------------------------------
    // Cut patterns and mappings
    //---------------------------------

    // Tolerance and vertex cut status
    enum VerStatus
    {
        In=0, // Vertex is on the +lSgnModifier side of the cut plane
        Out,  // Vertex is on the -lSgnModifier side of the cut plane
        Cut   // Vertex is within +/- lCutTol to the cut plane
    };

    // Patterns (unique patterns, using permutations below)
    static const VerStatus sCutPatterns[10][3] =
    {
        // 0 Out: 1+ In, 0+ Cut
/*0*/   { In, In, In   }, { Cut, In, In   }, { Cut, Cut, In  }, // Full triangle is IN
        // 0 In: 1+ Out, 0+ Cut,
/*3*/   { Out, Out, Out}, { Cut, Out, Out }, { Cut, Cut, Out }, // Full triangle is OUT
        // 0 In and 0 Out: 3 Cut
/*6*/   { Cut, Cut, Cut },                                      // Full triangle is OUT
        //
        // Remaining patterns must have 1+ In AND 1+ Out, screening by number of Cuts
        //
        // 0 Cuts: 1+ In and 1+ Out
/*7*/   { In, Out, Out }, { In, In, Out },                      // Triangle is cut (3 sub-tris)
        // 1 Cut: 1 In and 1 Out
/*9*/   { In, Out, Cut }                                        // Triangle is cut (2 sub-tris)
    };
    static const size_t sNbCutPatterns = sizeof( sCutPatterns ) / sizeof( *sCutPatterns );

    // Mappings
    static const int sMap[6][3] =
    {
        { 0, 1, 2 }, // +1: Identity
        { 0, 2, 1 }, // -1
        { 1, 0, 2 }, // -1
        { 1, 2, 0 }, // +1
        { 2, 0, 1 }, // +1
        { 2, 1, 0 }  // -1
    };
    static const int sMapDet[6] = { +1, -1, -1, +1, +1, -1 };
    static const int sInvMap[6][3] =
    {
        { 0, 1, 2 }, // Identity
        { 0, 2, 1 }, // Swaps 2 vertices only: inverse = self
        { 1, 0, 2 }, // Swaps 2 vertices only: inverse = self
        { 2, 0, 1 }, // Inverse of { 1, 2, 0 } (these 2 are swapped)
        { 1, 2, 0 }, // Inverse of { 2, 0, 1 } (these 2 are swapped)
        { 2, 1, 0 }  // Swaps 2 vertices only: inverse = self
    };

    /**
     * Mapping (1, 2, 0) transforms vertices (A, B, C) into (C, A, B)
     * The inverse mapping must transform (C, A, B) into (A, B, C)
     * C (at index 0) goes to index 2
     * A (at index 1) goes to index 0
     * B (at index 2) goes to index 1
     *
     * The inverse of (1, 2, 0) is thus (2, 0, 1)
     *
     * Or simply, if:
     * f(0) = 1
     * f(1) = 2
     * f(2) = 0
     *
     * then if g = inv(f):
     * g(0) = 2
     * g(1) = 0
     * g(2) = 1
     **/

#if 0
    //********************************
    //
    // BIG CHECK MAPPINGS
    //
    //********************************

    const Point3D ABC[] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0} };

    // Check mapping inversion
    for( int m=0 ; m<6 ; m++ )
    {
        for( int v=0 ; v<3 ; v++ )
        {
            if( sMap[m][ sInvMap[m][v] ] != v ) LzLogException("", "ERROR 1 m= "<<m)
            if( sInvMap[m][ sMap[m][v] ] != v ) LzLogException("", "ERROR 2 m= "<<m)
        }

        const Point3D mA = ABC[ sMap[m][0] ];
        const Point3D mB = ABC[ sMap[m][1] ];
        const Point3D mC = ABC[ sMap[m][2] ];
        const double lZ = ((mB - mA)^(mC - mA)).Z();

        LzLogM("", "Det= "<<lZ)
        LzLogM("", "Det= "<<sMapDet[m])

        if( lZ != sMapDet[m] )
            LzLogException("", "ERROR 3 m= "<<m)

        LzLogM("", "OK m= "<<m)
        LzLogM("", "-------------------")
    }
    return;
#endif

    //******************** ==> PARAMETERS (Options !)
    static const double sCutTol = 1e-6;

    // Check
    if( pFromVerIdx >= pMesh.mVertices.size() )
        LzLogException("", "Invalid vertex index "<<pFromVerIdx<<"! Max= "<<pMesh.mVertices.size())

    // Check
    if( !pTopo.IsCompatible(pMesh) )
        LzLogException("", "Provided topology is not compatible with the provided mesh!")

    //
    // Here: we have at least a non-empty mesh since vertices.size > 0
    //

    // Compute signe modifier
    const double lVerSgnDst = pCutPlane.SignedDistanceTo( pMesh.mVertices[pFromVerIdx] );

    // Check
    if( std::abs(lVerSgnDst) <= sCutTol )
        LzLogException("", "Cannot cut mesh from vertex! Starting vertex is ambiguous (signed distance= "
                           <<lVerSgnDst<<", tolerance= "<<sCutTol<<")")

    //
    // Here: the starting vertex will have a VerStatus either In or Out (but not Cut)
    //

    //---------------------------------
    // Init output meshes
    //---------------------------------

    // Setup detonator
    auto CLEAN_UP = [&]
        {
            pVerMesh.Free();
            pOppMesh.Free();
        };
    LzServices::Detonator lDet_SetUpCleanUpOutputMeshes( CLEAN_UP, CLEAN_UP );

    // Seed with init vertices
    pVerMesh.mVertices = pMesh.mVertices;
    pOppMesh.mVertices = pMesh.mVertices;

    // Store the number of initial vertices
    // This will be used to encore the indices of new vertices generated by the cut
    const size_t lInitVerSize = pMesh.mVertices.size();

    // We must be sure to have at least one normal since
    // we will be creating new triangles and assinging
    // normals to newly created vertices
    if( pMesh.mNormals.size() )
    {
        // Yes, use existing normals
        pVerMesh.mNormals = pMesh.mNormals;
        pOppMesh.mNormals = pMesh.mNormals;
    }
    else
    {
        // Push default normal to both meshes
        pVerMesh.mNormals = { Vector3D(1,0,0) };
        pOppMesh.mNormals = { Vector3D(1,0,0) };
    }

{ //-------- Artificial scope: freeing lEdgeInters via lDet_EdgeInters

    // List of new triangles resulting from cuts
    /**
     * These triangles will reference new vertices (resulting from cuts)
     * Each new vertex is associated to a unique existing edge
     * And since the final number of new vertices is unknown, a specific encoding
     * will be used new_ver_idx = lInitVerSize + cut_edge_index
     * Once all new vertices commited to both meshes, reindexation over all
     * new tris will be applied to recover the correct vertices indices.
     **/
    List<Triangle> lNewVerTris;
    List<Triangle> lNewOppTris;


    //---------------------------------
    // Compute cut orientation
    //---------------------------------

    // Sign modifier
    const double lSgnModifier = lVerSgnDst<0.0 ? -1.0 : +1.0 ;
//LzLogM("", "Sgn= "<<lSgnModifier)

    //
    // Factor applied to the signed distances of the vertices
    // to ensure that the VERTEX MESH is on the "positive" side of the cut plane
    //

    //---------------------------------
    // Triangles
    //---------------------------------

    // Triangle status (binary: could be replaced by a bool, however more explicit this way)
    enum TriStatus
    {
        Unreached=0,
        Pushed,
        // The Pushed state distinguishes Unreached triangles
        // from the ones on the stack but not yet processed
    };

    // Flags for processed triangles
    vector<TriStatus> lTriStatus;
    lTriStatus.assign( pMesh.mTriangles.size(), Unreached );

    //---------------------------------
    // Stack
    //---------------------------------

    // Stack for processing triangles
    List<size_t> lTriStack;

    // Seed stack with starting vertex
    const TopVer & lTV = pTopo.TopVers( pFromVerIdx );
    BrowseList( iT, lTV.mT )
    {
        // Triangle index
        const size_t lTriIdx = lTV.mT.GetAt(iT);

        // Push
        lTriStack.AddTail( lTriIdx );

        // Mark as processed
        lTriStatus[ lTriIdx ] = Pushed;
    }

    //---------------------------------
    // Edge intersections
    //---------------------------------

    struct InterPt
    {
        Point3D mPt;
        size_t mNewIdx;
    };
    vector<InterPt *> lEdgeInters;
    lEdgeInters.assign( pTopo.TopEdges().size(), nullptr );

    // Set up a non defusable detonator (will always detonate)
    LzServices::Detonator lDet_EdgeInters( [&]
        {
            for( size_t e=0 ; e<lEdgeInters.size() ; e++ )
                delete lEdgeInters[e];
        },
        LzServices::Detonator::Type::NonDefusable );

    //---------------------------------
    // Process tri stack
    //---------------------------------

    while( lTriStack.Count() )
    {
        // Pop a triangle from stack
        const size_t t = lTriStack.GetHead();
        lTriStack.DelHead();

//LzLogN("", "Triangle "<<t)

        // Get triangle info
        const Triangle & lTri = pMesh.mTriangles[t];
        //
        const Point3D * lpVers[3];
        for( int v=0 ; v<3 ; v++ )
            lpVers[v] = &pMesh.mVertices[ lTri.mIdxV[v] ];
        //
        const double lDists[3] =
        {
            lSgnModifier * pCutPlane.SignedDistanceTo( *lpVers[0] ),
            lSgnModifier * pCutPlane.SignedDistanceTo( *lpVers[1] ),
            lSgnModifier * pCutPlane.SignedDistanceTo( *lpVers[2] )
        };

        // Compute vertices status
        VerStatus lVerStats[3];
        for( int v=0 ; v<3 ; v++ )
        {
            if( lDists[v] > +sCutTol ) lVerStats[v] = In;
            else
            if( lDists[v] < -sCutTol ) lVerStats[v] = Out;
            else
                lVerStats[v] = Cut;
        }

        // Find the pattern and the mapping that fit the current triangle's configuration
        int lPatIdx, lMapIdx;
        for( size_t p=0 ; p<sNbCutPatterns ; p++ )
        {
            // Current pattern
            const VerStatus * const lpPat = sCutPatterns[p];

            // Try to find the mapping
            for( int m=0 ; m<6 ; m++ )
            {
                // Current mapping
                const int * const lpMap = sMap[m];

                // Found pattern?
                if( lVerStats[0] == lpPat[ lpMap[0] ]
                 && lVerStats[1] == lpPat[ lpMap[1] ]
                 && lVerStats[2] == lpPat[ lpMap[2] ] )
                {
                    // Yes: the pattern fits
                    lPatIdx = p;
                    lMapIdx = m;
                    goto FoundPattern;
                }

                // Nope...Try next mapping
            }
        }

        // Error
        LzLogException("", "Could not find suitable pattern and mapping for triangle "<<t<<"!")

FoundPattern:

//LzLogM("", "Triangle "<<t<<", Pat= "<<lPatIdx<<", Map= "<<lMapIdx)

        // Process pattern
        switch( lPatIdx )
        {
        //------------------------
        // Triangle IN
        //------------------------
        case 0:
        case 1:
        case 2:
            {
                // Store original triangle in Vertex Mesh
                pVerMesh.mTriangles.push_back( lTri );                

                // Get local topo
                const TopTri & lTTri = pTopo.TopTris( t );

                // Push triangles depending on specific case
                // 0: In ,  In, In --> Push all triangles
                // 1: Cut,  In, In --> Push all triangles
                // 2: Cut, Cut, In --> Push ONLY triangles opposite to I-C edges

                // Vertex indices C1 and C2: for pattern CCI only
                size_t lVerIdx_C1;
                size_t lVerIdx_C2;
                if( lPatIdx == 2 )
                {
                    // Find which index is I (at index 2 in the pattern)
                    const int lIdx_I = sInvMap[lMapIdx][ 2 ];

                    // Find the following C indices
                    const int lIdx_C1 = (lIdx_I + 1) % 3;
                    lVerIdx_C1 = lTri.mIdxV[ lIdx_C1 ];
                    //
                    const int lIdx_C2 = (lIdx_I + 2) % 3;
                    lVerIdx_C2 = lTri.mIdxV[ lIdx_C2 ];
                }

                // Check all 3 tri edges
                for( int e=0 ; e<3 ; e++ )
                {
                    // Get edge
                    const TopEdge & lTEd = pTopo.TopEdges( lTTri.mE[e] );

                    // Skip CC edge if pattern CCI
                    if( lPatIdx==2 && lTEd.HasVertex(lVerIdx_C1) && lTEd.HasVertex(lVerIdx_C2) )
                        continue;

                    // Multi-edge compatible edge browsing
                    // (since TopEdge::OtherTri throws an exception if an edge has more than 1 triangle)
                    const List<int> & lEdTris = lTEd.mT;

                    // Check all edge triangles
                    BrowseList( iEdT, lEdTris )
                    {

                        // No cracks here (at least the current triangle)
                        const int lOppTriIdx = lEdTris.GetAt(iEdT);

                        // Skip current triangle
                        if( lOppTriIdx == t )
                            continue;

                        // Consider only unreached triangles
                        // Ignore those who are: already Stacked or affected
                        // to one or another mesh
                        if( lTriStatus[lOppTriIdx] == Unreached  )
                        {
                            // Push to stack
                            lTriStatus[lOppTriIdx] = Pushed;
                            lTriStack.AddTail( lOppTriIdx );
                        }
                    }
                }
            }
            break;

        //------------------------
        // Triangle OUT
        //------------------------
        case 3:
        case 4:
        case 5:
            {
                // Store original triangle in Opposite Mesh
                pOppMesh.mTriangles.push_back( lTri );
            }
            break;

        //------------------------
        // Triangle OUT: C-C-C
        //------------------------
        case 6: // Cut-Cut-Cut
            {
                // Store original triangle in Opposite Mesh
                pOppMesh.mTriangles.push_back( lTri );
            }
            break;

        //------------------------
        // Triangle CUT: I-O1-O2
        //------------------------
        case 7:
            {
                // Find which index is I (at index 0 in the pattern)
                const int lIdx_I = sInvMap[lMapIdx][ 0 ];
                const size_t lVerIdx_I = lTri.mIdxV[ lIdx_I ];

                // Find the following O indices
                const int lIdx_O1 = (lIdx_I + 1) % 3;
                const size_t lVerIdx_O1 = lTri.mIdxV[ lIdx_O1 ];
                //
                const int lIdx_O2 = (lIdx_I + 2) % 3;
                const size_t lVerIdx_O2 = lTri.mIdxV[ lIdx_O2 ];

                // Find I-O1 and I-O2 edges in this triangle
                int lEd_IO1 = -1;
                int lEd_IO2 = -1;

                // Get local topo
                const TopTri & lTTri = pTopo.TopTris( t );

                // Check all 3 tri edges
                for( int e=0 ; e<3 ; e++ )
                {
                    // TopEdge
                    const TopEdge & lTEd = pTopo.TopEdges( lTTri.mE[e] );

                    // I-O1?
                    if( lTEd.HasVertex(lVerIdx_I) && lTEd.HasVertex(lVerIdx_O1) )
                        lEd_IO1 = lTTri.mE[e];
                    else
                    // I-O2?
                    if( lTEd.HasVertex(lVerIdx_I) && lTEd.HasVertex(lVerIdx_O2) )
                        lEd_IO2 = lTTri.mE[e];
                }

                // Didn't find both?
                if( lEd_IO1==-1 || lEd_IO1==-1 )
                    LzLogException("", "Could not find one of the cut edges in triangle "<<t<<"!"
                                       "(Pattern= "<<lPatIdx<<", Mapping= "<<lMapIdx<<")")

                // Compute intersections if they don't already exist

                // I-O1
                if( !lEdgeInters[lEd_IO1] )
                {
                    // Create new intersection point
                    InterPt * lpNewInter = new InterPt;

                    // Compute intersection position using signed distances
                    const Point3D & I  = *lpVers[ lIdx_I  ];
                    const Point3D & O1 = *lpVers[ lIdx_O1 ];
                    const Vector3D IO1 = O1 - I;
                    //
                    const double lDistRatio = (lDists[lIdx_I] / (lDists[lIdx_I] - lDists[lIdx_O1]));
                    lpNewInter->mPt = I + lDistRatio*IO1;

                    // Stash object; it will be deleted by detonator on exit
                    lEdgeInters[lEd_IO1] = lpNewInter;
                }

                // I-O2
                if( !lEdgeInters[lEd_IO2] )
                {
                    // Create new intersection point
                    InterPt * lpNewInter = new InterPt;

                    // Compute intersection position using signed distances
                    const Point3D & I  = *lpVers[ lIdx_I  ];
                    const Point3D & O2 = *lpVers[ lIdx_O2 ];
                    const Vector3D IO2 = O2 - I;
                    //
                    const double lDistRatio = (lDists[lIdx_I] / (lDists[lIdx_I] - lDists[lIdx_O2]));
                    lpNewInter->mPt = I + lDistRatio*IO2;

                    // Stash object; it will be deleted by detonator on exit
                    lEdgeInters[lEd_IO2] = lpNewInter;
                }

                //......................................
                // Create 1 new triangle on Ver side
                //......................................

                {
                    // NO NEED to take into account mapping orientation HERE
                    // since vertices order is recovered starting from the 'singular' vertex
                    // unique I or unique O
                    Triangle lVerTri(lVerIdx_I, lInitVerSize+lEd_IO1, lInitVerSize+lEd_IO2, (size_t)0, (size_t)0, (size_t)0 );

                    // Stash
                    lNewVerTris.AddTail( lVerTri );
                }

                // Push triangles on both sides of the cut edges

                // IO1
                {
                    // Check all triangles sharing this edge
                    const List<int> & lEdTris = pTopo.TopEdges(lEd_IO1).mT;
                    BrowseList( iET, lEdTris )
                    {
                        // Read triangle index
                        const int lTriIdx = lEdTris.GetAt(iET);

                        // Consider only unreached triangles
                        // Ignore those who are: already Stacked or affected
                        // to one or another mesh
                        if( lTriStatus[lTriIdx] == Unreached  )
                        {
                            // Push to stack
                            lTriStatus[lTriIdx] = Pushed;
                            lTriStack.AddTail( lTriIdx );
                        }
                    }
                }

                // IO2
                {
                    // Check all triangles sharing this edge
                    const List<int> & lEdTris = pTopo.TopEdges(lEd_IO2).mT;
                    BrowseList( iET, lEdTris )
                    {
                        // Read triangle index
                        const int lTriIdx = lEdTris.GetAt(iET);

                        // Consider only unreached triangles
                        // Ignore those who are: already Stacked or affected
                        // to one or another mesh
                        if( lTriStatus[lTriIdx] == Unreached  )
                        {
                            // Push to stack
                            lTriStatus[lTriIdx] = Pushed;
                            lTriStack.AddTail( lTriIdx );
                        }
                    }
                }

                //......................................
                // Create 2 new triangles on Opp side
                //......................................

                {
                    // NO NEED to take into account mapping orientation HERE
                    // since vertices order is recovered starting from the 'singular' vertex
                    // unique I or unique O
                    Triangle lOppTri_1(lVerIdx_O1, lInitVerSize+lEd_IO2, lInitVerSize+lEd_IO1, (size_t)0, (size_t)0, (size_t)0 );
                    Triangle lOppTri_2(lVerIdx_O1, lVerIdx_O2          , lInitVerSize+lEd_IO2, (size_t)0, (size_t)0, (size_t)0 );

                    // Stash
                    lNewOppTris.AddTail( lOppTri_1 );
                    lNewOppTris.AddTail( lOppTri_2 );
                }
            }
            break;

        //------------------------
        // Triangle CUT: I1-I2-O
        //------------------------
        case 8:
            {
                // Find which index is O (at index 2 in the pattern)
                const int lIdx_O = sInvMap[lMapIdx][ 2 ];
                const size_t lVerIdx_O = lTri.mIdxV[ lIdx_O ];

                // Find the following O indices
                const int lIdx_I1 = (lIdx_O + 1) % 3;
                const size_t lVerIdx_I1 = lTri.mIdxV[ lIdx_I1 ];
                //
                const int lIdx_I2 = (lIdx_O + 2) % 3;
                const size_t lVerIdx_I2 = lTri.mIdxV[ lIdx_I2 ];

                // Find O-I1 and O-I2 edges in this triangle
                int lEd_OI1 = -1;
                int lEd_OI2 = -1;

                // Get local topo
                const TopTri & lTTri = pTopo.TopTris( t );

                // Check all 3 tri edges
                for( int e=0 ; e<3 ; e++ )
                {
                    // TopEdge
                    const TopEdge & lTEd = pTopo.TopEdges( lTTri.mE[e] );

                    // I-O1?
                    if( lTEd.HasVertex(lVerIdx_O) && lTEd.HasVertex(lVerIdx_I1) )
                        lEd_OI1 = lTTri.mE[e];
                    else
                    // I-O2?
                    if( lTEd.HasVertex(lVerIdx_O) && lTEd.HasVertex(lVerIdx_I2) )
                        lEd_OI2 = lTTri.mE[e];
                }

                // Found both?
                if( lEd_OI1==-1 || lEd_OI2==-1 )
                    LzLogException("", "Could not find one of the cut edges in triangle "<<t<<"!"
                                       "(Pattern= "<<lPatIdx<<", Mapping= "<<lMapIdx<<")")

                // Compute intersections if they don't already exist

                // O-I1
                if( !lEdgeInters[lEd_OI1] )
                {
                    // Create new intersection point
                    InterPt * lpNewInter = new InterPt;

                    // Compute intersection position using signed distances
                    const Point3D & O  = *lpVers[ lIdx_O  ];
                    const Point3D & I1 = *lpVers[ lIdx_I1 ];
                    const Vector3D OI1 = I1 - O;
                    //
                    const double lDistRatio = (lDists[lIdx_O] / (lDists[lIdx_O] - lDists[lIdx_I1]));
                    lpNewInter->mPt = O + lDistRatio*OI1;

                    // Stash object; it will be deleted by detonator on exit
                    lEdgeInters[lEd_OI1] = lpNewInter;
                }

                // I-O2
                if( !lEdgeInters[lEd_OI2] )
                {
                    // Create new intersection point
                    InterPt * lpNewInter = new InterPt;

                    // Compute intersection position using signed distances
                    const Point3D & O  = *lpVers[ lIdx_O  ];
                    const Point3D & I2 = *lpVers[ lIdx_I2 ];
                    const Vector3D OI2 = I2 - O;
                    //
                    const double lDistRatio = (lDists[lIdx_O] / (lDists[lIdx_O] - lDists[lIdx_I2]));
                    lpNewInter->mPt = O + lDistRatio*OI2;

                    // Stash object; it will be deleted by detonator on exit
                    lEdgeInters[lEd_OI2] = lpNewInter;
                }

                //......................................
                // Create 2 new triangles on Ver side
                //......................................

                {
                    // NO NEED to take into account mapping orientation HERE
                    // since vertices order is recovered starting from the 'singular' vertex
                    // unique I or unique O
                    Triangle lVerTri_1(lVerIdx_I1, lInitVerSize+lEd_OI2, lInitVerSize+lEd_OI1, (size_t)0, (size_t)0, (size_t)0 );
                    Triangle lVerTri_2(lVerIdx_I2, lInitVerSize+lEd_OI2, lVerIdx_I1          , (size_t)0, (size_t)0, (size_t)0 );

                    // Stash
                    lNewVerTris.AddTail( lVerTri_1 );
                    lNewVerTris.AddTail( lVerTri_2 );
                }

                // Push triangles on both sides of all edges (2 'cut' edges and 1 'in' edge)

                // Check all 3 tri edges
                for( int e=0 ; e<3 ; e++ )
                {
                    // Multi-edge compatible edge browsing
                    // (since TopEdge::OtherTri throws an exception if an edge has more than 1 triangle)
                    const List<int> & lEdTris = pTopo.TopEdges( lTTri.mE[e] ).mT;

                    // Check all edge triangles
                    BrowseList( iEdT, lEdTris )
                    {
                        // No cracks here (at least the current triangle)
                        const int lOppTriIdx = lEdTris.GetAt(iEdT);

                        // Skip current triangle
                        if( lOppTriIdx == t )
                            continue;

                        // Consider only unreached triangles
                        // Ignore those who are: already Stacked or affected
                        // to one or another mesh
                        if( lTriStatus[lOppTriIdx] == Unreached  )
                        {
                            // Push to stack
                            lTriStatus[lOppTriIdx] = Pushed;
                            lTriStack.AddTail( lOppTriIdx );
                        }
                    }
                }

                //......................................
                // Create 1 new triangle on Opp side
                //......................................

                {
                    // NO NEED to take into account mapping orientation HERE
                    // since vertices order is recovered starting from the 'singular' vertex
                    // unique I or unique O
                    Triangle lOppTri(lVerIdx_O, lInitVerSize+lEd_OI1, lInitVerSize+lEd_OI2, (size_t)0, (size_t)0, (size_t)0 );

                    // Stash
                    lNewOppTris.AddTail( lOppTri );
                }
            }
            break;

        //------------------------
        // Triangle CUT: I-O-C
        //------------------------
        case 9:
            {
                // Find which index is I (at index 0 in the pattern)
                const int lIdx_I = sInvMap[lMapIdx][ 0 ];
                const size_t lVerIdx_I = lTri.mIdxV[ lIdx_I ];

                // Find which index is O (at index 1 in the pattern)
                const int lIdx_O = sInvMap[lMapIdx][ 1 ];
                const size_t lVerIdx_O = lTri.mIdxV[ lIdx_O ];

                // Find which index is C (at index 2 in the pattern)
                const int lIdx_C = sInvMap[lMapIdx][ 2 ];
                const size_t lVerIdx_C = lTri.mIdxV[ lIdx_C ];

                //
                // Here: the determinant of the mapping WILL influence
                // ----- the order of the vertices as I, O and C can be
                //       in either CW or CCW order
                //

                // Find I-O and I-C edges in this triangle
                int lEd_IO = -1;
                int lEd_IC = -1;

                // Get local topo
                const TopTri & lTTri = pTopo.TopTris( t );

                // Check all 3 tri edges
                for( int e=0 ; e<3 ; e++ )
                {
                    // TopEdge
                    const TopEdge & lTEd = pTopo.TopEdges( lTTri.mE[e] );

                    // I-O?
                    if( lTEd.HasVertex(lVerIdx_I) && lTEd.HasVertex(lVerIdx_O) )
                        lEd_IO = lTTri.mE[e];
                    else
                    // I-C?
                    if( lTEd.HasVertex(lVerIdx_I) && lTEd.HasVertex(lVerIdx_C) )
                        lEd_IC = lTTri.mE[e];
                }

                // Found both?
                if( lEd_IO==-1 || lEd_IC==-1 )
                    LzLogException("", "Could not find one of the cut edges in triangle "<<t<<"!"
                                       "(Pattern= "<<lPatIdx<<", Mapping= "<<lMapIdx<<")")

                // Compute I-O intersection if it doesn't already exist

                // I-O
                if( !lEdgeInters[lEd_IO] )
                {
                    // Create new intersection point
                    InterPt * lpNewInter = new InterPt;

                    // Compute intersection position using signed distances
                    const Point3D & I = *lpVers[ lIdx_I ];
                    const Point3D & O = *lpVers[ lIdx_O ];
                    const Vector3D IO = O - I;
                    //
                    const double lDistRatio = (lDists[lIdx_I] / (lDists[lIdx_I] - lDists[lIdx_O]));
                    lpNewInter->mPt = I + lDistRatio*IO;

                    // Stash object; it will be deleted by detonator on exit
                    lEdgeInters[lEd_IO] = lpNewInter;
                }

                //......................................
                // Create 1 new triangle on Ver side
                //......................................

                {
                    Triangle lVerTri;

                    // Take into account mapping orientation
                    if( sMapDet[lMapIdx] > 0 )
                        lVerTri = Triangle(lVerIdx_C, lVerIdx_I, lInitVerSize+lEd_IO, (size_t)0, (size_t)0, (size_t)0 );
                    else
                        lVerTri = Triangle(lVerIdx_I, lVerIdx_C, lInitVerSize+lEd_IO, (size_t)0, (size_t)0, (size_t)0 );

                    // Stash
                    lNewVerTris.AddTail( lVerTri );
                }

                // Push triangles on both sides of IO and IC (since only OC remains out

                // IO
                {
                    // Check all triangles sharing this edge
                    const List<int> & lEdTris = pTopo.TopEdges(lEd_IO).mT;
                    BrowseList( iET, lEdTris )
                    {
                        // Read triangle index
                        const int lTriIdx = lEdTris.GetAt(iET);

                        // Consider only unreached triangles
                        // Ignore those who are: already Stacked or affected
                        // to one or another mesh
                        if( lTriStatus[lTriIdx] == Unreached  )
                        {
                            // Push to stack
                            lTriStatus[lTriIdx] = Pushed;
                            lTriStack.AddTail( lTriIdx );
                        }
                    }
                }

                // IC
                {
                    // Check all triangles sharing this edge
                    const List<int> & lEdTris = pTopo.TopEdges(lEd_IC).mT;
                    BrowseList( iET, lEdTris )
                    {
                        // Read triangle index
                        const int lTriIdx = lEdTris.GetAt(iET);

                        // Consider only unreached triangles
                        // Ignore those who are: already Stacked or affected
                        // to one or another mesh
                        if( lTriStatus[lTriIdx] == Unreached  )
                        {
                            // Push to stack
                            lTriStatus[lTriIdx] = Pushed;
                            lTriStack.AddTail( lTriIdx );
                        }
                    }
                }

                //......................................
                // Create 1 new triangle on Opp side
                //......................................

                {
                    Triangle lOppTri;

                    // Take into account mapping orientation
                    if( sMapDet[lMapIdx] > 0 )
                        lOppTri = Triangle(lVerIdx_C, lInitVerSize+lEd_IO, lVerIdx_O, (size_t)0, (size_t)0, (size_t)0 );
                    else
                        lOppTri = Triangle(lVerIdx_O, lInitVerSize+lEd_IO, lVerIdx_C, (size_t)0, (size_t)0, (size_t)0 );

                    // Stash
                    lNewOppTris.AddTail( lOppTri );
                }
            }
            break;

        //------------------------
        default:
            LzLogException("", "Unexpected pattern index "<<lPatIdx<<"!")
        }
    }
    //--- end while( lTriStack.Count() )

    // Here: lTriStack is empty

    //---------------------------------
    // Finalize
    //---------------------------------

    // Add all new vertices to both meshes
    for( size_t e=0 ; e<lEdgeInters.size() ; e++ )
    {
        // Read intersection info
        InterPt * lpIP = lEdgeInters[e];

        // Have an intersection?
        if( lpIP )
        {
            // Store index for later reindexing
            lpIP->mNewIdx = pVerMesh.mVertices.size();

            // Push to both meshes
            pVerMesh.mVertices.push_back( lpIP->mPt );
            pOppMesh.mVertices.push_back( lpIP->mPt );
        }
    }
    LzLogM("", "Added "<<(pVerMesh.mVertices.size() - lInitVerSize)<<" new vertices to both meshes.")

    // Assemble Ver mesh
    LzLogM("", "Created "<<lNewVerTris.Count()<<" new triangles in VerMesh.")
    {
        // Reindex and push new triangles
        while( lNewVerTris.Count() )
        {
            // Pop triangle
            Triangle lTri = lNewVerTris.GetHead();
            lNewVerTris.DelHead();

            // Reindex
            for( int v=0 ; v<3 ; v++ )
            {
                // Get vertex index
                const size_t lVerIdx = lTri.mIdxV[v];

                // Need to reindex? (old or new vertex)
                if( lVerIdx > lInitVerSize )
                {
                    // Yes, new index -> need to reindex
                    const size_t lEdIdx = lVerIdx - lInitVerSize;
                    lTri.mIdxV[v] = lEdgeInters[lEdIdx]->mNewIdx;
                }
            }

            // Push to mesh
            pVerMesh.mTriangles.push_back( lTri );
        }

        // Here: lNewVerTris is empty
    }

    // Assemble Opp mesh
    LzLogM("", "Created "<<lNewOppTris.Count()<<" new triangles in OppMesh.")
    {
        // Add old triangles that have not been reached during the process
        for( size_t t=0 ; t<lTriStatus.size() ; t++ )
        {
            if( lTriStatus[t] == Unreached )
            {
                // Store triangle
                pOppMesh.mTriangles.push_back( pMesh.mTriangles[t] );
            }
        }

        // Here: we don't need lTriStatus anymore
        lTriStatus.clear();

        // Log
        LzLogM("", "Created "<<lNewOppTris.Count()<<" new triangles in OppMesh.")

        // Reindex and push new triangles
        while( lNewOppTris.Count() )
        {
            // Pop triangle
            Triangle lTri = lNewOppTris.GetHead();
            lNewOppTris.DelHead();

            // Reindex
            for( int v=0 ; v<3 ; v++ )
            {
                // Get vertex index
                const size_t lVerIdx = lTri.mIdxV[v];

                // Need to reindex? (old or new vertex)
                if( lVerIdx > lInitVerSize )
                {
                    // Yes, new index -> need to reindex
                    const size_t lEdIdx = lVerIdx - lInitVerSize;
                    lTri.mIdxV[v] = lEdgeInters[lEdIdx]->mNewIdx;
                }
            }

            // Push to mesh
            pOppMesh.mTriangles.push_back( lTri );
        }

        // Here: lNewOppTris is empty
    }

} //-------- Artificial scope: freeing lEdgeInters via lDet_EdgeInters

    //---------------------------------
    // Recompute normals
    //---------------------------------

//************ in options : pass topology if computed
    pVerMesh.SmoothNormalize();
    pOppMesh.SmoothNormalize();
//************ in options : pass topology if computed


    // All is well
    lDet_SetUpCleanUpOutputMeshes.defuse();
    // lDet_EdgeInters should not be defused!
}
#pragma endregion


#pragma region "Shape analysis"
/* BINARY FLAGS: https://m-peko.github.io/craft-cpp/posts/different-ways-to-define-binary-flags/
template <typename EnumT>
class Flags {
    static_assert(std::is_enum_v<EnumT>, "Flags can only be specialized for enum types");

    using UnderlyingT = typename std::make_unsigned_t<typename std::underlying_type_t<EnumT>>;

public:
    Flags& set(EnumT e, bool value = true) noexcept {
        bits_.set(underlying(e), value);
        return *this;
    }

    Flags& reset(EnumT e) noexcept {
        set(e, false);
        return *this;
    }

    Flags& reset() noexcept {
        bits_.reset();
        return *this;
    }

    [[nodiscard]] bool all() const noexcept {
        return bits_.all();
    }

    [[nodiscard]] bool any() const noexcept {
        return bits_.any();
    }

    [[nodiscard]] bool none() const noexcept {
        return bits_.none();
    }

    [[nodiscard]] constexpr std::size_t size() const noexcept {
        return bits_.size();
    }

    [[nodiscard]] std::size_t count() const noexcept {
        return bits_.count();
    }

    constexpr bool operator[](EnumT e) const {
        return bits_[underlying(e)];
    }

private:
    static constexpr UnderlyingT underlying(EnumT e) {
        return static_cast<UnderlyingT>(e);
    }

private:
    std::bitset<underlying(EnumT::size)> bits_;
};
*/
//================================================================================
bool CheckMeshProperties( const Mesh & pMesh, uint16_t pProp, const MeshTopology * ppTopo/*=nullptr*/ )
{
    // Pointer to topology: provided by client or my own personal one.
    const MeshTopology * lpTopo;
    //
    MeshTopology lMyOwnTopo;
    if( ppTopo )
    {
        // Check
        if( !ppTopo->IsCompatible(pMesh) )
            LzLogException("", "Cannot check mesh properties! The provided topology is not compatible with the mesh.")

        // Point to provided topology
        lpTopo = ppTopo;
    }
    else
    {
        // Set local topology
        lMyOwnTopo.Set( pMesh );

        // Point to local topology
        lpTopo = &lMyOwnTopo;
    }

    // Cracks
    if( pProp & MeshProperties::NO_CRACKS )
    {
        // Check
        if( lpTopo->Cracks().size() )
        {
            LzLogErr("", "Mesh should not have cracks! Found "<<lpTopo->Cracks().size()<<" crack(s).")
            return false;
        }
    }

    // Multi-edges
    if( pProp & MeshProperties::NO_MULTIEDGES )
    {
        // Check
        if( lpTopo->MultiEdges().size() )
        {
            LzLogErr("", "Mesh should not have multiedges! Found "<<lpTopo->MultiEdges().size()<<" multiedges(s).")
            return false;
        }
    }

    // Unused vertices
    if( pProp & MeshProperties::NO_UNUSED_VERS )
    {
        // Used vers table
        vector<bool> lUserVer;
        lUserVer.assign( pMesh.mVertices.size(), false );

        // Check all tris
        for( const Triangle & iTri : pMesh.mTriangles )
        {
            // Flag vertices as used
            lUserVer[ iTri.mIdxV[0] ] = true;
            lUserVer[ iTri.mIdxV[1] ] = true;
            lUserVer[ iTri.mIdxV[2] ] = true;
        }

        // Count unused vertices
        size_t lUnusedCount = 0;
        for( const bool iUsed : lUserVer )
            if( !iUsed ) lUnusedCount++;

        // Check
        if( lUnusedCount != 0 )
        {
            LzLogErr("", "Mesh should not have unused vertices! Found "<<lUnusedCount<<" unused vert(ex|ices).")
            return false;
        }
    }

    // Patches
    if( (pProp & MeshProperties::IS_SINGLE_PATCH)
     || (pProp & MeshProperties::IS_MULTI_PATCH) )
    {
        // Count patches
        const size_t lNbPatches = PatchComputer::CountPatches( *lpTopo );

        // Check
        if( (pProp & MeshProperties::IS_SINGLE_PATCH) && (lNbPatches != 1) )
        {
            LzLogErr("", "Mesh should have a single patch! Found "<<lNbPatches<<" patch(es).")
            return false;
        }

        // Check
        if( (pProp & MeshProperties::IS_MULTI_PATCH) && (lNbPatches == 1) )
        {
            LzLogErr("", "Mesh should be multi-patch! Found 1 patch.")
            return false;
        }
    }

    // Meters
    if( pProp & MeshProperties::IS_METERS )
    {
        // Check
        if( !pMesh.IsMeters() )
        {
            LzLogErr("", "Mesh should be in meters! It so appears that it is not.")
            return false;
        }
    }

    // Millimeters
    if( pProp & MeshProperties::IS_MILLIMETERS )
    {
        // Check
        if( pMesh.IsMeters() )
        {
            LzLogErr("", "Mesh should be in millimeters! It so appears that it is not.")
            return false;
        }
    }

    // Smooth normalized
    if( pProp & MeshProperties::IS_SMOOTH_NORM )
    {
        // Check
        if( pMesh.mVertices.size() != pMesh.mNormals.size() )
        {
            LzLogErr("", "Mesh should be smooth normlized! Nb. vertices= "<<pMesh.mVertices.size()<<" != nb. normals= "<<pMesh.mNormals.size()<<".")
            return false;
        }
    }

    return true;
}

//================================================================================
bool CheckTrianglesConsistency( const Mesh & pMesh,
                                double pMaxAngle_deg,
                                double pMinTriArea,
                                bool pListAllErrors,
                                const MeshTopology * ppTopo )
{
    // Result: true = everything is fine
    bool lCheckRes = true;

    // Pointer to topology: provided by client or my own personal one.
    const MeshTopology * lpTopo;
    //
    MeshTopology lMyOwnTopo;
    if( ppTopo )
    {
        // Check
        if( !ppTopo->IsCompatible(pMesh) )
            LzLogException("", "Cannot check mesh properties! The provided topology is not compatible with the mesh.")

        // Point to provided topology
        lpTopo = ppTopo;
    }
    else
    {
        // Set local topology
        lMyOwnTopo.Set( pMesh );

        // Point to local topology
        lpTopo = &lMyOwnTopo;
    }

    // Helper func
    auto GetAreaAndNormal = [&]( const Triangle & pTri, double & pArea, Vector3D & pNor )
        {
            // Get vertices
            const Point3D & A = pMesh.mVertices[ pTri.mIdxV[0] ];
            const Point3D & B = pMesh.mVertices[ pTri.mIdxV[1] ];
            const Point3D & C = pMesh.mVertices[ pTri.mIdxV[2] ];

            // Compute normal
            pNor = (B - A)^(C - A);

            // Compute area
            pArea = pNor.Norm() / 2.0;

            // Normalize normal
            pNor.Normalize();
        };


    // Check all triangles
    for( size_t t=0 ; t<pMesh.mTriangles.size() ; t++ )
    {
        // Get triangle
        const Triangle & lTri = pMesh.mTriangles[t];

        // Get area and normal
        double lArea;
        Vector3D lNormal;
        try
        {
            // Will throw exc if Normalize fails
            GetAreaAndNormal( lTri, lArea, lNormal );
        }
        catch(...)
        {
            // Log
            LzLogErr("", "Could not normalize normal in triangle "<<t<<"! The mesh seems to have poor quality.")
            if( !pListAllErrors )
                return false;
            lCheckRes = false;
        }

        //----------------------------------------------
        // Check area
        //----------------------------------------------

        if( lArea < pMinTriArea )
        {
            // Log
            LzLogErr("", "Triangle "<<t<<" is too small! Area= "<<lArea<<" < min area= "<<pMinTriArea<<".")
            if( !pListAllErrors )
                return false;
            lCheckRes = false;
        }

        //----------------------------------------------
        // Check angle with neighbor triangles
        //----------------------------------------------

        // Get topo tri
        const TopTri & lTTri = lpTopo->TopTris(t);

        // Check all neighbors across all edges (including multi-edges!)
        for( int e=0 ; e<3 ; e++ )
        {
            // Get top edge
            const TopEdge & lTEd = lpTopo->TopEdges( lTTri.mE[e] );

            // Check all opposite triangles across this edge
            BrowseList( iT, lTEd.mT )
            {
                // Get triangle index
                const size_t lOppIdx = lTEd.mT.GetAt( iT );

                // Skip current tri
                if( lOppIdx == t )
                    continue;

                // Get normal (and area... which we don't care about here)
                double lUnused;
                Vector3D lOppNormal;
                try
                {
                    // Will throw exc if Normalize fails
                    GetAreaAndNormal( pMesh.mTriangles[lOppIdx], lUnused, lOppNormal );
                }
                catch(...)
                {
                    // Log
                    LzLogErr("", "Could not normalize normal in (opposite) triangle "<<lOppIdx<<"! The mesh seems to have poor quality.")
                    if( !pListAllErrors )
                        return false;
                    lCheckRes = false;
                }

                // Here we have a normalized normal
                const double lAngle_deg = LzMath::RAD_2_DEG * lNormal.ShortPositiveAngleTo( lOppNormal );

                // Check angle
                if( lAngle_deg > pMaxAngle_deg )
                {
                    // Log
                    LzLogErr("", "Found an excessive angle between triangles "<<t<<" and "<<lOppIdx<<"! Angle= "<<lAngle_deg<<" > max angle= "<<pMaxAngle_deg<<" deg.")
                    if( !pListAllErrors )
                        return false;
                    lCheckRes = false;
                }
            }
        }
    }

    return lCheckRes;
}

//================================================================================
void PseudoSymPlaneOptions::EnableSubSamp( size_t pApproxSubSamps )
{
    // Check
    if( pApproxSubSamps < 3 )
        LzLogException("", "It seems silly to compute a symmetry plane using only "<<pApproxSubSamps<<" sub-samples!")

    // Set
    mSubSample = true;
    mApproxSubSamps = pApproxSubSamps;
}

//================================================================================
static Plane3D FindPseudoSymPlane( const Plane3D & pInitPlane, const vector<Point3D> & pPts,
                                   std::function<Point3D(const Point3D &)> pProjPt, const PseudoSymPlaneOptions & pOptions )
{
    // Check
    if( pPts.size() < 10 )
        LzLogException("", "Looking for a pseudo-symmetry plane in a points cloud or mesh with "
                           <<pPts.size()<<" point(s) is probably a bad idea!")

    // Adjust number of subsampled vertices
    vector<Point3D>         lLocalVers;
    const vector<Point3D> * lpUsedVers;
    if( pOptions.mSubSample )
    {
        //---> Subsample

        // Compute small step, which will yield a few samples more
        const size_t lStep = LzServices::GetSubSampSmallStep(pPts.size(), pOptions.mApproxSubSamps);

        // Resize array to contain actual number of samples
        lLocalVers.resize( LzServices::GetActualSamplesCount(pPts.size(), lStep) );

        // Read vertices
        for( size_t s=0 ; s<lLocalVers.size() ; s++ )
            lLocalVers[s] = pPts[ s * lStep ];

        // Point to used vertices
        lpUsedVers = &lLocalVers;
    }
    else
    {
        // Point to used vertices
        lpUsedVers = &pPts;
    }

    // Init guess
    Plane3D lSymPlane = pInitPlane;

    // Iterate
    vector<Point3D> lMidPts( lpUsedVers->size() );
    for( size_t i=0 ; i<pOptions.mMaxSteps ; i++ )
    {
        // Compute midpoints
        for( size_t v=0 ; v<lpUsedVers->size() ; v++ )
        {
            // Read vertex
            const Point3D & lV = (*lpUsedVers)[v];

            // Mirror point
            const Point3D lSymV = lSymPlane.Symmetrical( lV );

            // Project mirrored point onto object
            const Point3D lProjSymV = pProjPt( lSymV );

            // Stash midpoint
            lMidPts[v] = lV.MidPointTo( lProjSymV );
        }

        // Compute new plane
        Plane3D lNewPlane;
        if( pOptions.mConstrainSymNormal )
        {
            // Constrained symmetry plane

            // Project all midpoints onto constraint plane
            for( size_t m=0 ; m<lMidPts.size() ; m++ )
                lMidPts[m] = pOptions.mConstrainedNorPlane.Projection( lMidPts[m] );

            // Compute PCA in 2D
            Point3D lMean;
            Vector3D lV[3];
            LzMath::ToolBox::PCA( lMidPts, lMean, lV );

            // Assemble plane
            const Vector3D & lN = lV[1];
            const double A = lN.X();
            const double B = lN.Y();
            const double C = lN.Z();
            const double D = -( A * lMean.X() + B * lMean.Y() + C * lMean.Z() );
            //
            lNewPlane = Plane3D( A, B, C, D );
        }
        else
        {
            // Unconstrained symmetry plane
            lNewPlane = Plane3D( lMidPts );
        }

        // Check convergence
        if( SmallDiff(lNewPlane, lSymPlane, pOptions.mTol_mm, pOptions.mTol_deg) )
        {
            // Log
            LzLogM("", "Pseudo-symmetry plane found at step "<<i)
            LzLogM("", "Final plane: " << lSymPlane.ToString_ABCD())

            // return plane
            return lSymPlane;
        }

        // Move symmetry plane
        lSymPlane = lNewPlane;
    }

    // Error
    LzLogException("", "Symmetry plane not found after "<<pOptions.mMaxSteps<<" step(s)!")
}

//================================================================================
Plane3D FindPseudoSymPlane( const Plane3D & pInitPlane, const Mesh & pMesh, const FastDistComp & pFDC,
                            const PseudoSymPlaneOptions & pOptions/*=PseudoSymPlaneOptions()*/ )
{
    // Tri mesh-specific projector
    std::function<Point3D(const Point3D &)> ProjPt = [&]( const Point3D & pPt ) -> Point3D
        {
            // Project point onto object
            Point3D lProjPt;
            Vector3D lUnused;
            pFDC.SignedDistance(LzGeom::Tree3D::NearestMode::Accurate, pPt, lProjPt, lUnused);

            // Return projection
            return lProjPt;
        };

    // Call common implementation
    return FindPseudoSymPlane( pInitPlane, pMesh.mVertices, ProjPt, pOptions );
}

//================================================================================
Plane3D FindPseudoSymPlane( const Plane3D & pInitPlane, const vector<Point3D> & pPts, const Tree3D & pTree3D,
                            const PseudoSymPlaneOptions & pOptions/*=PseudoSymPlaneOptions()*/ )
{
    // Points-specific projector
    std::function<Point3D(const Point3D &)> ProjPt = [&]( const Point3D & pPt ) -> Point3D
        {
            // Project point onto object
            const unsigned int lNearestIdx = pTree3D.GetNearestVertex(pPt, Tree3D::NearestMode::Accurate);

            // Return projection
            return pPts[lNearestIdx];
        };

    // Call common implementation
    return FindPseudoSymPlane( pInitPlane, pPts, ProjPt, pOptions );
}


//================================================================================
#if 0 // Assessing differences between planes, based on rot and trans tolerance
for( size_t v=0 ; v<NB_VERTS ; v++ )
{
    LzLogTimeN("", "Processing vertebra "<<v)

    Mesh lVert_in_PCA = mSrcMeshes[v];

    const RigidTr3D & lAliasPCA_to_W = mSrc_AliasPCA_to_W[v];

    lVert_in_PCA.RigidTransform( lAliasPCA_to_W.Inverse() );

    // Compute FDC
    FastDistComp lFDC_in_PCA;
    lFDC_in_PCA.Set( lVert_in_PCA, mSrcTopos[v], LzTriModel::Patching::Split );

    // Options
    LzTriModel::PseudoSymPlaneOptions lOpt;

    // Find LOW RES sym plane
    Plane3D lP0 = LzTriModel::FindPseudoSymPlane(Plane3D(1,0,0,0), lVert_in_PCA, lFDC_in_PCA, lOpt);

    // Find HIHG RES sym plane
    lOpt.mTol_deg = 1e-6;
    lOpt.mTol_mm = 1e-6;
    lOpt.mMaxSteps = 500;
    Plane3D lP1 = LzTriModel::FindPseudoSymPlane(Plane3D(1,0,0,0), lVert_in_PCA, lFDC_in_PCA, lOpt);

    // Compute delta
    Point3D lProj0 = lP0.Projection( Point3D(0,0,0) );
    Point3D lProj1 = lP1.Projection( Point3D(0,0,0) );

    // Log
    {
        LzLogN("", "Delta")
        LzLogM("", "Distance delta= "<<lProj0.DistanceTo(lProj1)<<" mm")
        LzLogM("", "Angle delta= "<<LzMath::RAD_2_DEG * lP0.Normal().ShortPositiveAngleTo(lP1.Normal())<<" deg")
        LzLogM("", "")
        LzLogM("", lP0.ToString_ABCD())
        LzLogM("", lP1.ToString_ABCD())
    }
}
#endif

//================================================================================
static vector<double> ComputeSymScore( const Plane3D & pSymPlane, const vector<Point3D> & pPts,
                                       std::function<double(const Point3D &)> pPtDist, bool pLogValues )
{
    // Check
    if( pPts.size() < 10 )
        LzLogException("", "Scoring a pseudo-symmetry plane with a mesh with "<<pPts.size()<<" vertice(s) is probably a bad idea!")

    // Compute symmetry score
    vector<double> lMeanMaxStdDev( 3 );
    {
        // Collect errors
        vector<double> lErrs( pPts.size() );
        for( size_t v=0 ; v<pPts.size() ; v++ )
        {
            // Read vertex
            const Point3D & lV = pPts[v];

            // Mirror vertex
            const Point3D lSymV = pSymPlane.Symmetrical( lV );

            // Compute error for mirrored point
            lErrs[v] = pPtDist( lSymV );
        }

        // Stats
        LzMath::ToolBox::MeanMaxStdDev( lErrs, lMeanMaxStdDev[0], lMeanMaxStdDev[1], lMeanMaxStdDev[2] );

        // Log
        if( pLogValues )
        {
            LzLogN("", "Symmetry score")
            LzLogM("", "Mean=    "<<lMeanMaxStdDev[0])
            LzLogM("", "Max=     "<<lMeanMaxStdDev[1])
            LzLogM("", "Std dev= "<<lMeanMaxStdDev[2])
        }
    }

    return lMeanMaxStdDev;}

//================================================================================
vector<double> ComputeSymScore( const Plane3D & pSymPlane, const Mesh & pMesh,
                                const FastDistComp & pFDC, bool pLogValues )
{
    // Tri mesh-specific distance computer
    std::function<double(const Point3D &)> PtDist = [&]( const Point3D & pPt ) -> double
        {
            // Project point onto object
            return std::abs( pFDC.AccurateSignedDistance( pPt ) );
        };

    // Call common implementation
    return ComputeSymScore( pSymPlane, pMesh.mVertices, PtDist, pLogValues );
}

//================================================================================
vector<double> ComputeSymScore( const Plane3D & pSymPlane, const vector<Point3D> & pPts,
                                const Tree3D & pTree3D, bool pLogValues )
{
    // Points-specific distance computer
    std::function<double(const Point3D &)> PtDist = [&]( const Point3D & pPt ) -> double
        {
            // Project point onto object
            unsigned int lNearestIdx = pTree3D.GetNearestVertex(pPt, Tree3D::NearestMode::Accurate);
            return pPt.DistanceTo( pPts[lNearestIdx] );
        };

    // Call common implementation
    return ComputeSymScore( pSymPlane, pPts, PtDist, pLogValues );
}

//================================================================================
bool SmallDiff( const Plane3D & pP0, const Plane3D & pP1,
                double pTol_mm/*=1e-6*/, double pTol_deg/*=1e-6*/,
                const Point3D & pRefPt/*=Point3D(0,0,0)*/, bool pLog/*=false*/ )
{
    // Params 0
    const Point3D  lO_0 = pP0.Projection( pRefPt );
    const Vector3D lN_0 = pP0.Normal();

    // Param 1
    const Point3D lO_1 = pP1.Projection( pRefPt );
    Vector3D lN_1 = pP1.Normal();

    // Keep normal consistently oriented
    if( lN_0*lN_1 < 0 )
        lN_1 *= -1;

    // Compute deltas
    const double lDelta_mm  = lO_0.DistanceTo(lO_1);
    const double lDelta_deg = LzMath::RAD_2_DEG * lN_0.ShortPositiveAngleTo(lN_1);

    // Log
    if( pLog )
    {
        LzLogN("", "Difference between planes")
        LzLogM("", "Delta mm=  "<<lDelta_mm)
        LzLogM("", "Delta deg= "<<lDelta_deg)
    }

    // Return result
    return lDelta_mm<pTol_mm && lDelta_deg<pTol_deg;
}

//================================================================================
double ComputePartialVolume( const Mesh & pCutObj, const Point3D & pOrigin )
{
    // Const refs to stuff
    const vector<Triangle> & lTris = pCutObj.mTriangles;
    const vector<Point3D> & lVers = pCutObj.mVertices;

    // Init vol
    double lTotVol = 0.0;

    // Consider all triangles
    for( size_t t=0 ; t<lTris.size() ; t++ )
    {
        // Get triangle
        const Triangle & lTri = lTris[t];

        // Get points
        const Point3D & A = lVers[ lTri.mIdxV[0] ];
        const Point3D & B = lVers[ lTri.mIdxV[1] ];
        const Point3D & C = lVers[ lTri.mIdxV[2] ];

        // Compute tetrahedron vectors
        const Vector3D OA = A - pOrigin;
        const Vector3D OB = B - pOrigin;
        const Vector3D OC = C - pOrigin;

        // Add signed tet volume
        lTotVol += (OA^OB) * OC;
    }

    return lTotVol / 6.0;
}

//================================================================================
vector<double> ComputePartialVolAreaRatio( const Mesh & pObject,
                                           const Point3D & pOrigin, const Vector3D & pNormal )
{
    // Local usage
    using LzTriModel::Triangle;

//    //-------------------------
//    // Volume computer
//    //-------------------------

//    std::function<double(const Mesh &, const Point3D &)> CompPartVol =
//        []( const Mesh & pCutObj, const Point3D & pOrigin ) -> double
//        {
//            // Const refs to stuff
//            const vector<Triangle> & lTris = pCutObj.mTriangles;
//            const vector<Point3D> & lVers = pCutObj.mVertices;

//            // Init vol
//            double lTotVol = 0.0;

//            // Consider all triangles
//            for( size_t t=0 ; t<lTris.size() ; t++ )
//            {
//                // Get triangle
//                const Triangle & lTri = lTris[t];

//                // Get points
//                const Point3D & A = lVers[ lTri.mIdxV[0] ];
//                const Point3D & B = lVers[ lTri.mIdxV[1] ];
//                const Point3D & C = lVers[ lTri.mIdxV[2] ];

//                // Compute tetrahedron vectors
//                const Vector3D OA = A - pOrigin;
//                const Vector3D OB = B - pOrigin;
//                const Vector3D OC = C - pOrigin;

//                // Add signed tet volume
//                lTotVol += ((OA^OB) * OC) / 6.0;
//            }

//            return lTotVol;
//        };

//    //-------------------------
//    // Area computer
//    //-------------------------

//    std::function<double(const Mesh & pCutObj)> CompPartArea =
//        []( const Mesh & pCutObj ) -> double
//        {
//            // Const refs to stuff
//            const vector<Triangle> & lTris = pCutObj.mTriangles;
//            const vector<Point3D> & lVers = pCutObj.mVertices;

//            // Init area
//            double lTotArea = 0.0;

//            // Consider all triangles
//            for( size_t t=0 ; t<lTris.size() ; t++ )
//            {
//                // Get triangle
//                const Triangle & lTri = lTris[t];

//                // Get points
//                const Point3D & A = lVers[ lTri.mIdxV[0] ];
//                const Point3D & B = lVers[ lTri.mIdxV[1] ];
//                const Point3D & C = lVers[ lTri.mIdxV[2] ];

//                // Add area
//                lTotArea += ((B - A)^(C - A)).Norm();
//            }

//            return 0.5 * lTotArea;
//        };

    //-------------------------
    // Compute PVARs
    //-------------------------

    vector<double> lPVARs( 2 ); // [0] along +pNormal, [1] along -pNormal
    for( int i=0 ; i<2 ; i++ )
    {
        // Set cut plane
        const Plane3D lCutPlane( pOrigin, i==0 ? +pNormal : -pNormal );

        // Cut mesh
        Mesh lCutObj = pObject;
        lCutObj.Cut( lCutPlane, 1e-6 ); // No need to remove unused or duplicate

        // Compute Partial Volume and Partial Area
        const double lPV = ComputePartialVolume( lCutObj, pOrigin );
        const double lPA = lCutObj.Area(); //CompPartArea( lCutObj );

        // Check
        if( LzMath::ToolBox::IsZero(lPA) )
            LzLogException("", "Cannot compute partial volume to area ratio! Area is too small ("
                           <<lPA<<") on the "<<(i?"NEGATIVE":"POSITIVE")<<" side of the cut plane.")

        // Store ratio
        lPVARs[i] = lPV / lPA;
#if 1
        // Log
        LzLogM("", i<<": part vol=  "<<lPV)
        LzLogM("", i<<": part area= "<<lPA)
        LzLogM("", i<<": ratio=     "<<lPVARs[i])
#endif
    }

    return lPVARs;
}

//================================================================================
vector<double> ComputePartialVolAreaRatio( const vector<Point3D> & pVolPts, const vector<Point3D> & pSurfPts,
                                           const Point3D & pOrigin, const Vector3D & pNormal )
{
    //-------------------------
    // Compute PVARs
    //-------------------------

    // Compute partial volumes and areas
    vector<double> lPartVols{ 0, 0 };
    vector<double> lPartAreas{ 0, 0 };
    {
        // Cut plane
        const Plane3D lCutPlane( pOrigin, +pNormal );

        // Compute volumes
        for( const Point3D & iV : pVolPts )
        {
            if( lCutPlane.SignedDistanceTo(iV) >= 0 )
                lPartVols[0]++;
            else
                lPartVols[1]++;
        }

        // Compute areas
        for( const Point3D & iS : pSurfPts )
        {
            if( lCutPlane.SignedDistanceTo(iS) >= 0 )
                lPartAreas[0]++;
            else
                lPartAreas[1]++;
        }
    }

    // Check
    if( lPartAreas[0] == 0 )
        LzLogException("", "Cannot compute partial volume to area ratio! No voxels on the POSITIVE side of the cut.")

    // Check
    if( lPartAreas[1] == 0 )
        LzLogException("", "Cannot compute partial volume to area ratio! No voxels on the NEGATIVE side of the cut.")

    // Compute and return PVARs
    return { lPartVols[0]/lPartAreas[0],   // [0] along +pNormal
             lPartVols[1]/lPartAreas[1] }; // [1] along -pNormal
}

//================================================================================
void GetLineIntersections( const Mesh & pMesh,
                           const Point3D & pO, const Vector3D & pU,
                           ScoredPointsList & pSortedInters )
{
    // Remove previous
    pSortedInters.DelAll();

#if 1
    Eigen::Matrix2d M;
    Eigen::Vector2d B, XY;
#else
    Matrix M( 2, 2 ), M_inv( 2, 2 ), B( 2 ), XY( 2 ); // ==> 17 sec
#endif

    // Cut line
    const Line3D lLine( pO, pU );

    // Check all triangles
    for( size_t t=0 ; t<pMesh.mTriangles.size() ; t++ )
    {
        // Triangle and vertices
        const Triangle & lTri = pMesh.mTriangles[t];
        //
        const Point3D & lA = pMesh.mVertices[ lTri.mIdxV[0] ];
        const Point3D & lB = pMesh.mVertices[ lTri.mIdxV[1] ];
        const Point3D & lC = pMesh.mVertices[ lTri.mIdxV[2] ];

        // Compute intersection
        Point3D lI;
        if( lLine.Intersection_NoExc(Plane3D(lA, lB, lC), lI) )
        {
            // Compute local coordinates
            const Vector3D lX = lB - lA;
            const Vector3D lY = lC - lA;
            const Vector3D lU = lI - lA;
#if 1
            // Assemble matrix
            M( 0, 0 ) = lX * lX;
            M( 1, 1 ) = lY * lY;
            M( 0, 1 ) = M( 1, 0 ) = lX * lY;

            // Fill right hand side
            B( 0 ) = lX * lU;
            B( 1 ) = lY * lU;

            // Try to solve the system
            auto invM = M.fullPivLu(); // 2 to 3 times faster than auto invM = M.completeOrthogonalDecomposition();
            // in spite of what the doc says: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
            /* BENCHMARK
             *
                    const size_t SIZE = 10000;
                    vector<Eigen::Matrix2d> lMats;
                    vector<Eigen::Vector2d> lVecs;
                    lMats.resize( SIZE );
                    lVecs.reserve( SIZE );

                    for( size_t m=0 ; m<SIZE ; m++ )
                    {
                        Eigen::Vector2d A = Eigen::Vector2d::Random(2);
                        Eigen::Vector2d B = Eigen::Vector2d::Random(2);

                        lMats[m](0, 0) = A.dot(A);
                        lMats[m](1, 1) = B.dot(B);
                        lMats[m]( 0, 1 ) = lMats[m]( 1, 0 ) = A.dot(B);

                        lVecs[m] = Eigen::Vector2d::Random(2);
                    }

                    // fullPivLU
                    {
                        Eigen::Vector2d lRes;
                        lRes.setIdentity();

                        LzLogTimeN("", "fullPivLu")
                        for( size_t m=0 ; m<SIZE ; m++ )
                        {
                            auto invM = lMats[m].fullPivLu();
                            lRes = lRes + invM.solve(lVecs[m]);
                        }

                        LzLogM("", "Res=\n"<<lRes)
                    }

                    // completeOrthogonalDecomposition
                    {
                        Eigen::Vector2d lRes;
                        lRes.setIdentity();

                        LzLogTimeN("", "completeOrthogonalDecomposition")
                        for( size_t m=0 ; m<SIZE ; m++ )
                        {
                            auto invM = lMats[m].completeOrthogonalDecomposition();
                            lRes = lRes + invM.solve(lVecs[m]);
                        }

                        LzLogM("", "Res=\n"<<lRes)
                    }
             */
            if( invM.isInvertible() )
            {
                XY = invM.solve(B);

                // Check local coordinates
                const double xx = XY.x();
                const double yy = XY.y();
                if( xx>=-1e-6 && yy>=-1e-6 && (xx + yy)<=1+1e-6 )
                    pSortedInters.AddIntoIncrList( LzServices::Sortable<Point3D,double>(lI, (lI - pO)*pU), true );
            }
#else
            // Assemble matrix
            M.Elt( 0, 0 ) = lX * lX;
            M.Elt( 1, 1 ) = lY * lY;
            M.Elt( 0, 1 ) = M.Elt( 1, 0 ) = lX * lY;

            // Fill right hand side
            B.Elt( 0 ) = lX * lU;
            B.Elt( 1 ) = lY * lU;

            // Try to solve the system
            if( M.M2x2_InverseTo_NoExc(M_inv) )
            {
                M_inv.Mult( B, XY );

                // Check local coordinates
                if( XY.Elt( 0 )>=-1e-6 && XY.Elt( 1 )>=-1e-6 && (XY.Elt( 0 ) + XY.Elt( 1 ))<=1+1e-6 )
                    pSortedInters.AddIntoIncrList( LzServices::Sortable<Point3D,double>(lI, (lI - pO)*pU), true );
            }
#endif
        }
    }
}
#pragma endregion


#pragma region "Geometry generation"
//================================================================================
void GenEllipsoid( Mesh & pTo, const Point3D & pCenter,
                   double pXRadius, double pYRadius, double pZRadius,
                   size_t pMeridianSampling/*=30*/, size_t pParallelSampling/*=30*/ )
{
	// Check
	if( pMeridianSampling<2 || pParallelSampling<3 )
        LzLogException("", "This will never look like an ellipsoid!... Meridian sampling= "<<pMeridianSampling
                           <<", parallel sampling= "<<pParallelSampling<<".")

	// Free previous
	pTo.Free();


	//--------------------------------
#pragma region "Init"
	//--------------------------------

	// Sampling steps
    double lMeridianStep = LzMath::PI / pMeridianSampling;
    double lParallelStep = LzMath::DPI / pParallelSampling;

	// Count everybody
    const size_t lVerNorCount = 2/* poles */
                              + pParallelSampling * (pMeridianSampling-1)/* parallels */;
    const size_t lFaceCount = 2 * pParallelSampling /* poles */
                            + (2*pParallelSampling)*(pMeridianSampling-2) /* between parallels */;
	// Resize
	pTo.mVertices.resize( lVerNorCount );
	pTo.mNormals.resize( lVerNorCount );
	pTo.mTriangles.resize( lFaceCount );
#pragma endregion


	//--------------------------------
#pragma region "Vertices / normals"
	//--------------------------------

	// Pole X-
    pTo.mVertices[0] = Point3D( pCenter.X()-pXRadius, pCenter.Y(), pCenter.Z() );
	pTo.mNormals[0]  = Vector3D( -1.0, 0.0, 0.0 );

	// Vertices on parallels
    size_t iCurrNorVer = 1;
    for( double iM=+LzMath::PI-lMeridianStep ; iM>+lMeridianStep/2.0 ; iM-=lMeridianStep )
    for( double iP=0.0 ; iP<LzMath::DPI-lParallelStep/2.0 ; iP+=lParallelStep )
	{
		// Vertex
        pTo.mVertices[iCurrNorVer].X() = pCenter.X() + pXRadius*cos(iM);
        pTo.mVertices[iCurrNorVer].Y() = pCenter.Y() + pYRadius*sin(iM)*cos(iP);
        pTo.mVertices[iCurrNorVer].Z() = pCenter.Z() + pZRadius*sin(iM)*sin(iP);

		// Normal
		{
			// Partial derivatives
			Vector3D lDerM( -pXRadius*sin(iM), +pYRadius*cos(iM)*cos(iP), +pZRadius*cos(iM)*sin(iP) );
			Vector3D lDerP( +0.0, -pYRadius * sin(iM) * sin(iP), +pZRadius * sin(iM) * cos(iP) );

			// Cross prod
			pTo.mNormals[iCurrNorVer] = lDerM ^ lDerP;

			// Normalize
			double lNorm = pTo.mNormals[iCurrNorVer].Norm();
            if( !LzMath::ToolBox::IsZero(lNorm) )
				pTo.mNormals[iCurrNorVer] /= lNorm;
		}

		// Next patient please
		iCurrNorVer++;
	}

	// pole X+
    pTo.mVertices[iCurrNorVer] = Point3D( pCenter.X()+pXRadius, pCenter.Y(), pCenter.Z() );
	pTo.mNormals[iCurrNorVer]  = Vector3D( +1.0, 0.0, 0.0 );
#pragma endregion


	//--------------------------------
#pragma region "Faces"
	//--------------------------------

    size_t iCurrFace = 0;

	// Pole X-
    for( size_t iV=1 ; iV<=pParallelSampling ; iV++ )
	{
		Triangle & lTri = pTo.mTriangles[iCurrFace];

		lTri.mIdxV[0] = lTri.mIdxN[0] = iV;
		lTri.mIdxV[1] = lTri.mIdxN[1] = 0;

		int l_NextF = (iV==pParallelSampling) ? 1 : iV+1 ;
		lTri.mIdxV[2] = lTri.mIdxN[2] = l_NextF;

		iCurrFace++;
	}

	// Parallels
    for( size_t iP=0 ; iP<pMeridianSampling-2 ; iP++ )
    for( size_t i_V=1 ; i_V<=pParallelSampling ; i_V++ )
	{
        size_t l_Next_i_V			= (i_V==pParallelSampling) ? 1 : i_V+1 ;
        size_t l_CurrVertexOnMeri	= iP*pParallelSampling + i_V;
        size_t l_NextVertexOnMeri	= iP*pParallelSampling + l_Next_i_V;
        size_t l_OppVertex			= (iP+1)*pParallelSampling + i_V;
        size_t l_NextOppVertex		= (iP+1)*pParallelSampling + l_Next_i_V;

		// Tri 1
		{
			Triangle & lTri = pTo.mTriangles[iCurrFace];

			lTri.mIdxV[0] = lTri.mIdxN[0] = l_CurrVertexOnMeri;
			lTri.mIdxV[1] = lTri.mIdxN[1] = l_NextVertexOnMeri;
			lTri.mIdxV[2] = lTri.mIdxN[2] = l_OppVertex;
		}
		iCurrFace++;

		// Tri 2
		{
			Triangle & lTri = pTo.mTriangles[iCurrFace];

			lTri.mIdxV[0] = lTri.mIdxN[0] = l_OppVertex;
			lTri.mIdxV[1] = lTri.mIdxN[1] = l_NextVertexOnMeri;
			lTri.mIdxV[2] = lTri.mIdxN[2] = l_NextOppVertex;
		}
		iCurrFace++;
	}

	// Pole X+
    const size_t lLastStart = pParallelSampling*(pMeridianSampling-2)+1;
    const size_t lLastEnd   = pParallelSampling*(pMeridianSampling-1);
	//
    for( size_t iV=lLastStart ; iV<=lLastEnd ; iV++ )
	{
		Triangle & lTri = pTo.mTriangles[iCurrFace];

		lTri.mIdxV[0] = lTri.mIdxN[0] = iV;

		int lNextF = (iV==lLastEnd) ? lLastStart : iV+1 ;
		lTri.mIdxV[1] = lTri.mIdxN[1] = lNextF;
		lTri.mIdxV[2] = lTri.mIdxN[2] = lVerNorCount-1;

		iCurrFace++;
	}
#pragma endregion
}

//================================================================================
void GenCylinder( Mesh & pTo, const Point3D & pCenter1, double pRadius1, const Point3D & pCenter2, double pRadius2, size_t pSegments, size_t pLayers )
{
    // Check
    if( pSegments < 2 )
        LzLogException("", "Cannot build cylinder with less than 2 segments!")

    // Check
    if( pLayers < 1 )
        LzLogException("", "Cannot build cylinder with less than 1 layer!")

    // Set-up detonator
    LzServices::Detonator lDet( [&] { pTo.Free(); } );

    // Referential at Center1
    const Vector3D lVec12 = pCenter2 - pCenter1;
    Vector3D lX, lY, lZ;
    lVec12.BuildOrthonormalBase( lX, lY, lZ );

    // Free previous
    pTo.Free();

    // Steps
    const double lAStep = 2.0 * LzMath::PI / pSegments;
    const double lSStep = lVec12.Norm() / pLayers;


    //----------------------------
    // Vertices
    //----------------------------

    {
        // Alloc
        pTo.mVertices.resize( 2/*poles*/ + (pLayers + 1)*pSegments );

        // Set poles
        pTo.mVertices[0] = pCenter1;
        pTo.mVertices[1] = pCenter2;
        //
        size_t lCurrVer = 2;

        // Generate all remaining vertices and normals
        for( size_t iSlice=0 ; iSlice<=pLayers ; iSlice++ )
        {
            // Interpolate slice center and radius
            const Point3D lSliceO = pCenter1 + iSlice*lSStep*lZ;
            //
            const double lAlpha = iSlice / (double)pLayers; // 0 at iSlice=0 (Center 1) ; 1 at iSlice = pLayers (Center 2)
            const double lSliceRad = (1.0 - lAlpha)*pRadius1 + lAlpha*pRadius2;

            // Generate circle
            for( size_t iSeg=0 ; iSeg<pSegments ; iSeg++ )
            {
                const double lA = iSeg * lAStep;
                const Vector3D lU = std::cos(lA)*lX + std::sin(lA)*lY;

                pTo.mVertices[lCurrVer++] = lSliceO + lSliceRad*lU;
            }
        }

        // Check
        if( lCurrVer != pTo.mVertices.size() )
            LzLogException("", "Something is messed up!")
    }


    //----------------------------
    // Normals
    //----------------------------

    {
        // Radial normal vector
        // Looking from the side at a horizontal cylinder, Center1 ---> Center2; place a XY ref on a cyl's vertex
        // and compute the normal tilt in this "sagittal plane"
        Vector3D lCoordsN( pRadius1 - pRadius2, pCenter1.DistanceTo( pCenter2 ), 0 );
        lCoordsN.Normalize();
        const Vector3D lNZ_Z = lCoordsN.X() * lZ;
        const double lNU = lCoordsN.Y();

        // Alloc
        pTo.mNormals.resize( 2/*poles*/ + pSegments/*cyl edges and sides*/ );

        // Set poles
        pTo.mNormals[0] = -lZ;
        pTo.mNormals[1] = +lZ;
        //
        size_t lCurrNor = 2;

        // Generate circle
        for( size_t iSeg=0 ; iSeg<pSegments ; iSeg++ )
        {
            const double lA = iSeg * lAStep;
            const Vector3D lU = std::cos(lA)*lX + std::sin(lA)*lY;

            pTo.mNormals[lCurrNor++] = lNZ_Z + lNU*lU;
        }

        // Check
        if( lCurrNor != pTo.mNormals.size() )
            LzLogException("", "Something is messed up!")
    }


    //----------------------------
    // Triangles
    //----------------------------

    {
        // Alloc
        pTo.mTriangles.resize( 2*pSegments/*at poles*/ + 2*pLayers*pSegments );
        //
        size_t lCurrTri = 0;

        // Find first circle vertices
        const size_t lFirstVer_1 = 2;
        const size_t lFirstVer_2 = pTo.mVertices.size() - pSegments;

        // Disks
        for( size_t iSeg=0 ; iSeg<pSegments ; iSeg++ )
        {
            // Disk 1
            pTo.mTriangles[lCurrTri++] = Triangle( (size_t)0, lFirstVer_1+(iSeg+1)%pSegments, lFirstVer_1+iSeg, (size_t)0, (size_t)0, (size_t)0);

            // Disk 2
            pTo.mTriangles[lCurrTri++] = Triangle( (size_t)1, lFirstVer_2+iSeg, lFirstVer_2+(iSeg+1)%pSegments, (size_t)1, (size_t)1, (size_t)1);
        }

        // Sides
        for( size_t iSlice=0 ; iSlice<pLayers ; iSlice++ )
        {
            // Base vertex index on current slice
            const size_t lV0 = 2 + iSlice*pSegments;
            const size_t lV1 = lV0 + pSegments;

            for( size_t iSeg=0 ; iSeg<pSegments ; iSeg++ )
            {
                // Normal index
                const size_t lN0 = 2 + iSeg;
                const size_t lN1 = 2 + (iSeg+1)%pSegments;

                // Side triangles
                pTo.mTriangles[lCurrTri++] = Triangle( lV0+iSeg, lV1+(iSeg+1)%pSegments, lV1+iSeg, lN0, lN1, lN0 );
                pTo.mTriangles[lCurrTri++] = Triangle( lV0+iSeg, lV0+(iSeg+1)%pSegments, lV1+(iSeg+1)%pSegments, lN0, lN1, lN1 );
            }
        }

        // Check
        if( lCurrTri != pTo.mTriangles.size() )
            LzLogException("", "Something is messed up!")
    }

    // All is well
    lDet.defuse();
}
#pragma endregion


#pragma region "Error measurement"
//================================================================================
void AssessDelta_A_to_B( const Mesh & pA, const Mesh & pB, double & pMeanErr, double & pMaxErr, double & pStdDev,
                         const FastDistComp * ppFDC_B/*=nullptr*/, const MeshTopology * ppTopo_B/*=nullptr*/ )
{
    // Check
    if( ppFDC_B && ppTopo_B )
        LzLogException("", "Conflicting parameters! Specify a FDC or a topology for B, but not both.");

    // Local FDC and avatar
    FastDistComp lFDC;
    const FastDistComp * lpFDC;

    // Compute errors
    Vector<double> lErrs( pA.mVertices.size() );
    size_t i = 0;
    //
    // From A
    //
    {
        // Set temp FDC for B
        if( ppFDC_B )
            lpFDC = ppFDC_B;
        else
        {
            if( ppTopo_B )
                lFDC.Set( pB, *ppTopo_B, LzTriModel::Patching::Split );
            else
                lFDC.Set( pB, LzTriModel::Patching::Split );

            // Link to local FDC
            lpFDC = &lFDC;
        }

        // Compute error
        for( size_t v=0 ; v<pA.mVertices.size() ; v++ )
        {
            Point3D tmp1;
            Vector3D tmp2;
            double lDist = lpFDC->SignedDistance( Tree3D::NearestMode::Accurate, pA.mVertices[v], tmp1, tmp2);
            lErrs[i++] = fabs( lDist );
        }
    }

    //	appliquer a TxAnatomy::Registration
    LzMath::ToolBox::MeanMaxStdDev( lErrs, pMeanErr, pMaxErr, pStdDev );
}

//================================================================================
void AssessBilateralDelta( const Mesh & pA, const Mesh & pB, double & pMeanErr, double & pMaxErr, double & pStdDev,
						   const FastDistComp * ppFDC_A/*=nullptr*/, const FastDistComp * ppFDC_B/*=nullptr*/,
						   const MeshTopology * ppTopo_A/*=nullptr*/, const MeshTopology * ppTopo_B/*=nullptr*/ )
{
	// Check
	if( ppFDC_A && ppTopo_A )
        LzLogException("", "Conflicting parameters! Specify a FDC or a topology for A, but not both.");

	// Check
	if( ppFDC_B && ppTopo_B )
        LzLogException("", "Conflicting parameters! Specify a FDC or a topology for B, but not both.");

    // Check
    if( pA.mVertices.size()==0 || pB.mVertices.size()==0 )
        LzLogException("", "Cannot compute bilateral delta when one (or both) mesh(es) is(are) empty!")

    // Local FDC and avatar
	FastDistComp lFDC;
	const FastDistComp * lpFDC;

	// Compute errors
	Vector<double> lErrs( pA.mVertices.size() + pB.mVertices.size() );
    size_t i = 0;
	//
	// From A
	//
	{
		// Set temp FDC for B
		if( ppFDC_B )
			lpFDC = ppFDC_B;
		else
		{
			if( ppTopo_B )
                lFDC.Set( pB, *ppTopo_B, Patching::Split );
			else
                lFDC.Set( pB, Patching::Split );

			// Link to local FDC
			lpFDC = &lFDC;
		}

		// Compute error
        for( size_t v=0 ; v<pA.mVertices.size() ; v++ )
		{
            const double lDist = lpFDC->AccurateSignedDistance( pA.mVertices[v] );
            lErrs[i++] = fabs( lDist );
		}
	}
	//
	// From B
	//
	{
		// Set temp FDC for A
		if( ppFDC_A )
			lpFDC = ppFDC_A;
		else
		{
			if( ppTopo_A )
                lFDC.Set( pA, *ppTopo_A, Patching::Split );
			else
                lFDC.Set( pA, Patching::Split );

			// Link to local FDC
			lpFDC = &lFDC;
		}

		// Compute error
        for( size_t v=0 ; v<pB.mVertices.size() ; v++ )
		{
            const double lDist = lpFDC->AccurateSignedDistance( pB.mVertices[v] );
            lErrs[i++] = fabs( lDist );
		}
	}

	//	appliquer a TxAnatomy::Registration
    LzMath::ToolBox::MeanMaxStdDev( lErrs, pMeanErr, pMaxErr, pStdDev );
}
#pragma endregion


#pragma region "ICP fit"
//================================================================================
RigidTr3D Rigid_ICP( const Mesh & pSrcMesh, const Mesh & pDstMesh, double pMaxEulerRotDeg/*=0.1*/, double pMaxTransNorm/*=0.1*/, size_t pMaxSteps/*=300*/  )
{
	// No affine pose
	Matrix lIdPose( 4, 4 );
	lIdPose.LoadIdentity();

	// Compute reg
	Matrix lRigidMat = Rigid_ICP( pSrcMesh, pDstMesh, lIdPose, pMaxEulerRotDeg, pMaxTransNorm, pMaxSteps );

	// Convert matrix to rigid transform and return
	RigidTr3D lRes;
	lRes.FromMatrix4x4( lRigidMat );
	return lRes;
}

//================================================================================
Matrix Rigid_ICP( const Mesh & pSrcMesh, const Mesh & pDstMesh, const Matrix & pAffinePose, double pMaxEulerRotDeg/*=0.1*/, double pMaxTransNorm/*=0.1*/, size_t pMaxSteps/*=300*/ )
{
	// Copy source mesh
	Mesh lSrcMesh = pSrcMesh;

	// Apply pose
	lSrcMesh.AffineTransform( pAffinePose );

	// Compute ICP
	RigidTr3D lICP;

	// Create tree 3D
    LzGeom::Tree3D lDstTree;
	lDstTree.Create( pDstMesh.mVertices, false );
    lICP = LzMath::ToolBox::ArunLoop( lSrcMesh.mVertices, lDstTree, pMaxEulerRotDeg, pMaxTransNorm, pMaxSteps );

	// Append initial pose to ICP matrix and return
	Matrix M, M_Pose;
	lICP.ToMatrix4x4( M );
	M.Mult( pAffinePose, M_Pose );

	return M_Pose;
}
#pragma endregion


#pragma region "PCA fit"
//================================================================================
Matrix RegisterMinMax( const Mesh & pSrc, const Mesh & pDst, const RigidTr3D & pSrc2Wrld, const RigidTr3D & pDst2Wrld )
{
	// Compute source PCA bbox range
	double lSrcW[3];
	{
		Vector<double> lMin, lMax;

		const Vector3D lVs[] = { pSrc2Wrld*Vector3D(1,0,0), pSrc2Wrld*Vector3D(0,1,0), pSrc2Wrld*Vector3D(0,0,1) };
        LzMath::ToolBox::PCAMinMax( pSrc.mVertices, Point3D(0,0,0), lVs, lMin, lMax );

		for( int i=0 ; i<3 ; i++ )
			lSrcW[i] = lMax[i] - lMin[i];
	}

	// Compute destination PCA bbox range
	double lDstW[3];
	{
		Vector<double> lMin, lMax;

		const Vector3D lVs[] = { pDst2Wrld*Vector3D(1,0,0), pDst2Wrld*Vector3D(0,1,0), pDst2Wrld*Vector3D(0,0,1) };
        LzMath::ToolBox::PCAMinMax( pDst.mVertices, Point3D(0,0,0), lVs, lMin, lMax );

		for( int i=0 ; i<3 ; i++ )
			lDstW[i] = lMax[i] - lMin[i];
	}

	// Build transform matrix
	Matrix W2A, P2W, S(4,4);
	pSrc2Wrld.Inverse().ToMatrix4x4( W2A );
	pDst2Wrld.ToMatrix4x4( P2W );
	S.LoadIdentity();
	for( int i=0 ; i<3 ; i++ )
		S.Elt(i,i) = lDstW[i] / lSrcW[i];

	// Assembly
	Matrix S_W2A, P2W_S_W2A;
	S.Mult( W2A, S_W2A );
	P2W.Mult( S_W2A, P2W_S_W2A );

	return P2W_S_W2A;
}

//=============================================================================
Matrix AffinePCAFitPairedMeshes( Mesh * pSrc, const Mesh * pDst, size_t pMaxIter, double pTolToIdentity, List<size_t> * ppLockedScales/*=nullptr*/,
								 const FastDistComp * ppDstFDC/*=nullptr*/, const MeshTopology * ppDstTopo/*=nullptr*/ )
{
//Vector<Mesh*> lSrc( 1 );
//lSrc[0] = pSrc;
//
//Vector<const Mesh*> lDst( 1 );
//lDst[0] = pDst;

	// Check
	if( ppDstFDC && ppDstTopo )
        LzLogException("", "Conflicting parameters! Specify a FDC or a topology, but not both.");

    if( ppDstFDC )
	{
		Vector<const FastDistComp*> tmp({ ppDstFDC });
        return AffinePCAFitPairedMeshes( Vector<Mesh*>({ pSrc }), { pDst }, pMaxIter, pTolToIdentity, ppLockedScales, &tmp );
	}
	else
	if( ppDstTopo )
	{
		FastDistComp lDstFDC;
        lDstFDC.Set( *pDst, *ppDstTopo, LzTriModel::Patching::Split );
		Vector<const FastDistComp*> tmp({ &lDstFDC });
        return AffinePCAFitPairedMeshes( Vector<Mesh*>({ pSrc }), Vector<const Mesh*>({ pDst }), pMaxIter, pTolToIdentity, ppLockedScales, &tmp );
	}
	else
        return AffinePCAFitPairedMeshes( Vector<Mesh*>({ pSrc }), Vector<const Mesh*>({ pDst }), pMaxIter, pTolToIdentity, ppLockedScales );
}

//=============================================================================
Matrix AffinePCAFitPairedMeshes( const Vector<Mesh*> & pSrc, const Vector<const Mesh*> & pDst, size_t pMaxIter, double pTolToIdentity,
                                 List<size_t> * ppLockedScales/*=nullptr*/,
                                 const Vector<const FastDistComp*> * pDstFDCs/*=nullptr*/,
                                 const Vector<const MeshTopology*> * pDstTopos/*=nullptr*/ )
{
	// Check
	if( pSrc.Size() != pDst.Size() )
        LzLogException("", "Meshes count mismatch! Cannot register "<<pSrc.Size()<<" source and "<<pDst.Size()<<" destination mesh(es).");

	// Count vertices
    size_t lTotSrcVerCount = 0;
    for( size_t s=0 ; s<pSrc.Size() ; s++ )
        lTotSrcVerCount += (size_t)pSrc[s]->mVertices.size();

	// Check
	if( lTotSrcVerCount < 3 )
        LzLogException("", "Cannot rigidly register less than 3 points! Total vertices= "<<lTotSrcVerCount<<".");

	// Check
	if( pDstFDCs && pDstTopos )
        LzLogException("", "Conflicting parameters! Specify a FDC or a topology, but not both.");

	// Create Fast Distance Computers
	Vector<FastDistComp> lFDCs;
	if( pDstFDCs )
	{
		// Check only
		if( pDstFDCs->Size() != pDst.Size() )
            LzLogException("", "Mismatch between number of FDCs and number of destination meshes! "<<pDstFDCs->Size()<<" != "<<pDst.Size()<<".");
	}
	else
	{
		// Create FDCs from scratch
		lFDCs.Resize( pDst.Size() );

		// Use provided topologies?
		if( pDstTopos )
		{
			// Check
			if( pDstTopos->Size() != pDst.Size() )
                LzLogException("", "Mismatch between number of topos and number of destination meshes! "<<pDstTopos->Size()<<" != "<<pDst.Size()<<".");

			// Use topos
            for( size_t d=0 ; d<pDst.Size() ; d++ )
                lFDCs[d].Set( *pDst[d], *(*pDstTopos)[d], LzTriModel::Patching::Split );
		}
		else
		{
			// From scratch
            for( size_t d=0 ; d<pDst.Size() ; d++ )
                lFDCs[d].Set( *pDst[d], LzTriModel::Patching::Split );
		}
	}

	// Best deformation
//double lBestQuadErr = -1; // Uninitialized
double lBestQuadErr = std::numeric_limits<double>::max(); // Uninitialized
	Matrix lBestSrc2Dst;
    size_t lBestIteration;

	// Current deformation
	Matrix lCurrSrc2Dst( 4, 4 );
	lCurrSrc2Dst.LoadIdentity();

	// Source and paired points
	Vector<Point3D> lSrc( lTotSrcVerCount );
	Vector<Point3D> lDst( lTotSrcVerCount );

	// Init source points once and for all
	// ==> Will be updated at each iteration
	{
        size_t i = 0;
        for( size_t s=0 ; s<pSrc.Size() ; s++ )
		{
			// Get source mesh
            const Mesh & lSrcMesh = *pSrc[s];

			// Store source point
            for( size_t v=0 ; v<lSrcMesh.mVertices.size() ; v++ )
				lSrc[i++] = lSrcMesh.mVertices[v];
		}
	}

	// Compute PCA on source triangles: PCA is shape specific and not mesh specific (if PCA was done on vertices and not on weighted triangles)
	// ==> Will be updated at each iteration
	Point3D lSrcG;
	Vector3D lSrcV[3];
	{
		// Count all triangles
        size_t lTotTriCount = 0;
        for( size_t s=0 ; s<pSrc.Size() ; s++ )
			lTotTriCount += pSrc[s]->mTriangles.size();

		// Init tables
		Vector<Point3D> lAllCentroids( lTotTriCount );
		Vector<double> lAllAreas( lTotTriCount );

		// Compute all centroids
        size_t i = 0;
        for( size_t s=0 ; s<pSrc.Size() ; s++ )
		{
			// Get from mesh
			Vector<Point3D> lCentroids;
			Vector<double> lAreas;
			pSrc[s]->GetFacesCentroids( lCentroids, lAreas );

			// Commit
            for( size_t c=0 ; c<lCentroids.Size() ; c++ )
			{
				lAllCentroids[i + c] = lCentroids[c];
				lAllAreas    [i + c] = lAreas    [c];
			}

			// Update start
			i += lCentroids.Size();
		}

		// Compute weighted PCA
        LzMath::ToolBox::WeightedPCA( lAllCentroids, lAllAreas, lSrcG, lSrcV );
	}

	//----- The loop
    size_t iIter = 0;
	while( true )
	{
		// Pair points and compute error before registration
		double lErr = 0;
        size_t i = 0;
        for( size_t s=0 ; s<pSrc.Size() ; s++ )
		{
			// Get source mesh
            const Mesh & lSrcMesh = *pSrc[s];

			// Get relevant FDC
			const FastDistComp & lFDC = pDstFDCs ? *(*pDstFDCs)[s] : lFDCs[s] ;

			// Pair all source vertices with the corresponding destination mesh
            for( size_t v=0 ; v<lSrcMesh.mVertices.size() ; v++ )
			{
#if 1
				// Pair with surface point and not vertex
				Point3D lProjSrc;
                Vector3D tmp;
                double lD = fabs( lFDC.SignedDistance(LzGeom::Tree3D::NearestMode::Accurate, lSrcMesh.mVertices[v], lProjSrc, tmp) );

				// Update error
				lErr += lD*lD;

				// Store
				lDst[i] = lProjSrc;
#else
				*** DEBUG *** (in case of FastDistComp failure)
				*** DEBUG *** Associates each source vertex to a destination vertex (and not a surface point)

				// Get destination mesh
				const Mesh & lDstMesh = *pDst[s];

				// Pair with nearest vertex using exhaustive search
				double lMinDist = +std::numeric_limits<double>::max();
                for( size_t w=0 ; w<lDstMesh.mVertices.size() ; w++ )
				{
					double lDist = lSrcMesh.mVertices[v].DistanceTo( lDstMesh.mVertices[w] );
					if( lMinDist > lDist )
					{
						lMinDist = lDist;
						lDst[i] = lDstMesh.mVertices[w];
					}
				}

				// Store
				lSrc[i] = lSrcMesh.mVertices[v];

				//*** DEBUG *** (in case of FastDistComp failure)
#endif
				i++;
			}
		}

//// Log
//LzLogM("", "Affine fit error at "<<iIter<<" = "<<lErr)

		// Best so far?
        if( lBestQuadErr > lErr )
		{
			lBestQuadErr = lErr;
			lBestSrc2Dst = lCurrSrc2Dst;
			lBestIteration = iIter;
		}
//        else
//        {
//            // Log
//            LzLogM("", "Affine fit error stable or increasing at "<<iIter<<", from "<<lBestQuadErr<<" --> "<<lErr)
//
//            // Break from loop
//            break;
//        }

		// End of loop?
		if( iIter == pMaxIter )
			break;

		// Compute Affine increment
		Matrix lAffIncr;
		{
			RigidTr3D lArun;
			Matrix lAffMat;

			// Affine increment is identity?
			//
			//*** Here lSrc, lSrcG and lSrcV are updated by the iteration
			//
            if( LzMath::ToolBox::AffinePCAFitPairedPoints_Iteration(lSrc, lDst, lArun, lAffMat, pTolToIdentity, lSrcG, lSrcV, ppLockedScales) )
				break;

//// Log
//{
//    LzLogN("", "*** OLD INCREMENT")
//    lArun.Log("Arun");
//    lAffMat.Log("AffMat");
//}


			// Convert to matrix 4x4
			Matrix lArunMat;
			lArun.ToMatrix4x4( lArunMat );
			lAffMat.Mult( lArunMat, lAffIncr );
		}

		// Affine deform all source meshes
        for( size_t s=0 ; s<pSrc.Size() ; s++ )
			pSrc[s]->AffineTransform( lAffIncr );

		// Update affine matrix
		Matrix lNewSrc2Dst;
		lAffIncr.Mult( lCurrSrc2Dst, lNewSrc2Dst );
		lCurrSrc2Dst = lNewSrc2Dst;

		// Next iteration
		iIter++;
	}

	// Log
    LzLogMsg("", "Fit paired meshes after "<<iIter<<" iteration(s). Min. quadratic error= "<<lBestQuadErr<<" found at iteration #"<<lBestIteration<<".");

	return lBestSrc2Dst;
}
#pragma endregion


#pragma region "Elastic fit"
/* BACKUP

From Algorithmicus / ModelEditor

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
System::Void Form1::DEBUG15TSMI_Click(System::Object^  sender, System::EventArgs^  e)
{
	_TxLogTry
	_TxWaitCursor;

	// Split
	List<Mesh> lPatches;
	PatchComputer lPatchComp;
	lPatchComp.SplitMesh( mGLWin->mpMesh->mT, lPatches );

	// Check
	if( lPatches.Count() != 2 )
        LzLogException("", "Need 2 patches!");

	// Create aliases
	Mesh & lSrc = lPatches.GetHead();
	Mesh & lDst = lPatches.GetTail();

	// Normalize
	lSrc.SmoothNormalize();
	lDst.SmoothNormalize();

	// Save
    lSrc.Save( LzServices::StartUpPath() + "/DEBUG5_SRC.obj" );
    lDst.Save( LzServices::StartUpPath() + "/DEBUG5_DST.obj" );

	// Build source topology
	MeshTopology lSrcTopo;
	lSrcTopo.SetMesh( lSrc );

	// Vector of dual displacements
	Vector<Vector3D> lPrimalDisps( lSrcTopo.TopVers().size() );
	Vector<Vector3D> lDualDisps( lSrcTopo.TopTris().size() );

	// Build destination dist comp
	FastDistComp lDstFDC;
    lDstFDC.Set( lDst, LzTriModel::Patching::Split );

	// Iterate N times
const size_t N = 10;
const size_t SMOOTHER = 1;
//const size_t SMOOTHER = 10; une valeur trop grande induit des distorsions !!!!
    for( size_t n=0 ; n<N ; n++ )
	{
		// Compute primal displacements
        for( size_t v=0 ; v<lSrc.mVertices.size() ; v++ )
		{
			const Point3D & lV = lSrc.mVertices[v];
			Point3D lProjV;
            lDstFDC.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV, lProjV, Vector3D() );
			lPrimalDisps[v] = lProjV - lV;
		}

		// Smoothing loop: going back and forth between primal and dual
        for( size_t s=0 ; s<SMOOTHER ; s++ )
		{
			// Compute dual displacements
            for( size_t t=0 ; t<lSrc.mTriangles.size() ; t++ )
			{
				lDualDisps[t].Reset();
				for( int v=0 ; v<3 ; v++ )
					lDualDisps[t] += lPrimalDisps[ lSrc.mTriangles[t].mIdxV[v] ];

				lDualDisps[t] /= 3.0;
			}

			// Back to primal
            for( size_t v=0 ; v<lSrc.mVertices.size() ; v++ )
			{
				lPrimalDisps[v].Reset();
                const List<size_t> & lVerTris = lSrcTopo.TopVers()[v].mT;
				BrowseList( iT, lVerTris )
					lPrimalDisps[v] += lDualDisps[ lVerTris.GetAt(iT) ];

				lPrimalDisps[v] /= lVerTris.Count();
			}
		}

		// Apply resulting primal
        for( size_t v=0 ; v<lSrc.mVertices.size() ; v++ )
			lSrc.mVertices[v] += lPrimalDisps[v];

		// Normalize and save
		lSrc.SmoothNormalize();
        lSrc.Save( LzServices::StartUpPath() + "/DEBUG5_SRC_" + std::to_string(n) + ".obj" );
	}

	_TxLogFinalCatch("Could not DEBUG 15!",Manager::GetGlobal()->GetStdString("ERR Error"))
	finally
	{
		// Update log window
		Map::ActivePostTo( "LogFile", TxLinks::Messages::UPDATE );
	}
}



void ElasticFit( Mesh & pFrom, const Mesh & pTo, size_t pFitIter, size_t pSmoothIter, bool pBilateral )
{
	// Hommage a S. Kubrick
	Mesh & lFlyBone = pFrom;

	// Topo and distance computer
	MeshTopology lSrcTopo;
	lSrcTopo.SetMesh( lFlyBone );
	//
	FastDistComp lFDC_To;
    lFDC_To.Set( pTo, LzTriModel::Patching::Split );
	//
	// Tree3D for inverse registration
	Tree3D lTree_From;

	Vector<Vector3D> lFinalDisps;
	{
		// Reset final displacements
		lFinalDisps.Assign( lFlyBone.mVertices.size(), Vector3D(0,0,0) );

		// Apply incremental displacements
		Vector<Vector3D> lIncrDisps( lFlyBone.mVertices.size() );
        for( size_t n=0 ; n<pFitIter ; n++ )
		{
			//------------------------------
			// Compute displacements
			//------------------------------

			// Bilateral registration enabled?
			if( pBilateral )
			{
				// Bilateral displacements
				Vector< List<Vector3D> > lDispList( lFlyBone.mVertices.size() );

				// Src to Dst
                for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				{
					const Point3D & lV = lFlyBone.mVertices[v];

					Point3D lProjV;
                    lFDC_To.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV, lProjV, Vector3D() );
					lDispList[v].AddTail( lProjV - lV );
				}

				// Init tree
				lTree_From.Create( lFlyBone.mVertices );

				// Dst to Src
                for( size_t v=0 ; v<pTo.mVertices.size() ; v++ )
				{
					// Dst vertex
					const Point3D & lDstVer = pTo.mVertices[v];

                    const size_t lD2S = lTree_From.GetNearestVertex( lDstVer, Tree3D::NearestMode::Accurate );
					lDispList[lD2S].AddTail( lDstVer - lFlyBone.mVertices[lD2S] );
				}

				// Take average
                for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				{
					lIncrDisps[v].Reset();
					BrowseList( iD, lDispList[v] )
						lIncrDisps[v] += lDispList[v].GetAt(iD);

					lIncrDisps[v] /= lDispList[v].Count();
				}
			}
			else
			{
				// Src to Dst only
                for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				{
					const Point3D & lV = lFlyBone.mVertices[v];

					Point3D lProjV;
                    lFDC_To.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV, lProjV, Vector3D() );
					lIncrDisps[v] = lProjV - lV ;
				}
			}

			//------------------------------
			// Smooth displacements
			//------------------------------
			{
				Vector<Vector3D> lNewDisps( lIncrDisps.Size() );

                for( size_t s=0 ; s<pSmoothIter ; s++ )
				{
					// Smooth
                    for( size_t v=0 ; v<lIncrDisps.Size() ; v++ )
					{
						// Reset
						lNewDisps[v].Reset();

						// Average
                        const LzTriModel::TopVer & lTV = lSrcTopo.TopVers()[v];
						if( lTV.mE.Count() )
						{
							BrowseList( iE, lTV.mE )
							{
								int lOtherV = lSrcTopo.TopEdges()[ lTV.mE.GetAt(iE) ].OtherVer( v );
								lNewDisps[v] += lIncrDisps[ lOtherV ];
							}

							// Normalize
							lNewDisps[v] /= lTV.mE.Count();
						}
					}

					// Commit
					lIncrDisps = lNewDisps;
				}
			}

			// Apply resulting primal
            for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				lFlyBone.mVertices[v] += lIncrDisps[v];
		}
	}

	// Normalize
	lFlyBone.SmoothNormalize();
}

*/

//================================================================================
void ElasticFit( Mesh & pFrom, const Mesh & pTo, size_t pFitIter, size_t pSmoothIter, bool pBilateral,
                 List<size_t> * ppSrcIdx/*=nullptr*/, List<Point3D> * ppDstPos/*=nullptr*/ )
{
    // Check
    if( ppSrcIdx && ppDstPos )
    {
        if( ppSrcIdx->Count() != ppDstPos->Count() )
            LzLogException("", "Inconsistent displacement constraints! "<<ppSrcIdx->Count()<<" != "<<ppDstPos->Count()<<".")
    }
    else
    if( ppSrcIdx || ppDstPos )
        LzLogException("", "Inconsistent displacement constraints!")

    // Create constraints table
    Vector<bool> lIsConstrained;
    Vector<Point3D> lConstrainedPos;
    if( ppSrcIdx && ppDstPos )
    {
        // Size and set
        lIsConstrained.Assign( pFrom.mVertices.size(), false );
        lConstrainedPos.Resize( pFrom.mVertices.size() );

        BrowseTwoLists( iIdx, *ppSrcIdx, iPos, *ppDstPos )
        {
            size_t v = ppSrcIdx->GetAt(iIdx);

            lIsConstrained[v] = true;
            lConstrainedPos[v] = ppDstPos->GetAt(iPos);
        }
    }

    // Hommage a S. Kubrick
	Mesh & lFlyBone = pFrom;

	// Topo and distance computer
	MeshTopology lSrcTopo;
	lSrcTopo.Set( lFlyBone );
	//
	// FDC for inverse registration
	FastDistComp lFDC_Src;
	//
	FastDistComp lFDC_Dst;
    lFDC_Dst.Set( pTo, LzTriModel::Patching::Split );

	Vector<Vector3D> lFinalDisps;
	{
		// Reset final displacements
		lFinalDisps.Assign( lFlyBone.mVertices.size(), Vector3D(0,0,0) );

		// Apply incremental displacements
		Vector<Vector3D> lIncrDisps( lFlyBone.mVertices.size() );
        for( size_t n=0 ; n<pFitIter ; n++ )
		{
			//------------------------------
			// Compute displacements
			//------------------------------

			// Bilateral registration enabled?
			if( pBilateral )
			{
				// Bilateral displacements
				Vector< List<Vector3D> > lDispList( lFlyBone.mVertices.size() );

				// Src to Dst
                for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				{
					// Src vertex
					const Point3D & lV = lFlyBone.mVertices[v];

					// Project onto Dst
					Point3D lProjV;
					Vector3D tmp;
                    lFDC_Dst.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV, lProjV, tmp);
					lDispList[v].AddTail( lProjV - lV );
				}

				// Init FDC
                lFDC_Src.Set( lFlyBone, lSrcTopo, LzTriModel::Patching::Split );

				// Dst to Src
                for( size_t v=0 ; v<pTo.mVertices.size() ; v++ )
				{
					// Dst vertex
					const Point3D & lV = pTo.mVertices[v];

					// Project onto Src
					FastDistComp::SgnDstInfo lInfo;
					Point3D lProjV;
					Vector3D tmp;
                    lFDC_Src.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV, lProjV, tmp, &lInfo );

//*** DEBUG
//lInfo.Log();
//lInfo.LogProj( lFlyBone.mVertices, lProjV );
//*** DEBUG

					//-------------------------------
#pragma region "      Distribute displacement"
					//-------------------------------

					// Displacement
					const Vector3D lSrc2Dst = lV - lProjV;

					// Distribution
					switch( lInfo.mTarget )
					{
                    case FastDistComp::SgnDstInfo::Target::A:
                    case FastDistComp::SgnDstInfo::Target::B:
                    case FastDistComp::SgnDstInfo::Target::C:
					{
                        lDispList[ lInfo.mV[(size_t)lInfo.mTarget] ].AddTail( lSrc2Dst );
					}
					break;
					//
                    case FastDistComp::SgnDstInfo::Target::AB:
                    case FastDistComp::SgnDstInfo::Target::AC:
                    case FastDistComp::SgnDstInfo::Target::BC:
					{
						// Edge index
                        const size_t lEI = (size_t)lInfo.mTarget - 3;

						// Vertex indices for each edge
						// (uint)AB = 3
						// (uint)AC = 4
						// (uint)BC = 5
                        static const size_t sEdges[3][2] = { {0, 1}, {0, 2}, {1, 2} };

						// Read weights
						const double lW0 = lInfo.mW[ sEdges[lEI][0] ];
						const double lW1 = lInfo.mW[ sEdges[lEI][1] ];

						// (AAt)^-1 = 1.0/Alpha * Id
						const double lAlpha = lW0*lW0 + lW1*lW1;

						// Compute displacemants: At(AAt)^-1 U
						const Vector3D lDisp0 = (lW0 / lAlpha) * lSrc2Dst;
						const Vector3D lDisp1 = (lW1 / lAlpha) * lSrc2Dst;

						// Dispatch
						lDispList[ lInfo.mV[ sEdges[lEI][0] ] ].AddTail( lDisp0 );
						lDispList[ lInfo.mV[ sEdges[lEI][1] ] ].AddTail( lDisp1 );
					}
					break;
					//
                    case FastDistComp::SgnDstInfo::Target::ABC:
					{
						// Spread the butter
						Matrix A( 3, 9 );
						A.LoadValue( 0.0 );
						for( int w=0 ; w<3 ; w++ )
						for( int i=0 ; i<3 ; i++ )
							A.Elt( i, 3*w + i ) = lInfo.mW[w];

						// Compute pseudo inverse PI = At(AAt)^-1
						Matrix PI;
						{
							Matrix At;
							A.TransposeTo( At );
							Matrix AAt;
							A.Mult( At, AAt );
							Matrix inv_AAt;
							AAt.M3x3_InverseTo( inv_AAt );
							At.Mult( inv_AAt, PI );
						}

						// Solve
						Matrix X;
						{
							// Convert displacement to matrix
							Matrix Y( 3, 1, {lSrc2Dst.mV[0], lSrc2Dst.mV[1], lSrc2Dst.mV[2]} );

							// Solve
							PI.Mult( Y, X );
						}

						// Distribute displacement
						for( int w=0 ; w<3 ; w++ )
						{
							const Vector3D lDisp( X.Elt(w*3 + 0), X.Elt(w*3 + 1), X.Elt(w*3 + 2) );
							lDispList[ lInfo.mV[w] ].AddTail( lDisp );
						}
					}
					break;
					//
					default:
                        LzLogException("", "Unexpected signed distance target!");
					};
#pragma endregion
				}

				// Take average
                for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				{
					lIncrDisps[v].Reset();
					BrowseList( iD, lDispList[v] )
						lIncrDisps[v] += lDispList[v].GetAt(iD);

					lIncrDisps[v] /= lDispList[v].Count();
				}
			}
			else
			{
				// Src to Dst only
                for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				{
					const Point3D & lV = lFlyBone.mVertices[v];

                    // Constrained?
                    if( ppSrcIdx && ppDstPos && lIsConstrained[v] )
                    {
                        // Override with constraints
                        lIncrDisps[v] = lConstrainedPos[v] - lV;
                    }
                    else
                    {
                        // Free projection
                        Point3D lProjV;
                        Vector3D tmp;
                        lFDC_Dst.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lV, lProjV, tmp);
                        lIncrDisps[v] = lProjV - lV;
                    }
				}
			}

			//------------------------------
			// Smooth displacements
			//------------------------------
			{
				Vector<Vector3D> lNewDisps( lIncrDisps.Size() );

                for( size_t s=0 ; s<pSmoothIter ; s++ )
				{
					// Smooth
                    for( size_t v=0 ; v<lIncrDisps.Size() ; v++ )
					{
                        // Skip constrained
                        if( ppSrcIdx && ppDstPos && lIsConstrained[v] )
                        {
                            lNewDisps[v] = lIncrDisps[v];
                            continue;
                        }

						// Reset
						lNewDisps[v].Reset();

						// Average
                        const LzTriModel::TopVer & lTV = lSrcTopo.TopVers()[v];
						if( lTV.mE.Count() )
						{
							BrowseList( iE, lTV.mE )
							{
								int lOtherV = lSrcTopo.TopEdges()[ lTV.mE.GetAt(iE) ].OtherVer( v );
								lNewDisps[v] += lIncrDisps[ lOtherV ];
							}

							// Normalize
							lNewDisps[v] /= lTV.mE.Count();
						}
					}

                    // Commit all
                    lIncrDisps = lNewDisps;
				}
			}

			// Apply resulting primal
            for( size_t v=0 ; v<lFlyBone.mVertices.size() ; v++ )
				lFlyBone.mVertices[v] += lIncrDisps[v];
		}
	}

	// Normalize
	lFlyBone.SmoothNormalize();
}
#pragma endregion


#pragma region "Extending from cracks"
//================================================================================
class VerAndIdx
{
public:
    VerAndIdx( size_t pVer ) : mVer(pVer) {}
	bool operator>=( const VerAndIdx & pOther ) const { return mVer >= pOther.mVer; }
	bool operator!=( const VerAndIdx & pOther ) const { return mVer != pOther.mVer; }

    size_t mVer; // Vertex index
    size_t mIdx; // Index of this VerAndIdx within a list
};

//=============================================================================
void ExtendFromCracks( Mesh & pMesh, const MeshTopology * ppMeshTopo,
					   const Plane3D & pExtPlane,
                       size_t pLayersCount, double pLayerHeight,
                       size_t pSamplesCount, double pSampleHeight,
                       size_t pLaplaceRelaxSteps )
{
	// Check
	if( pLayersCount == 0 )
		return;

	// Check
	if( pSamplesCount < 2 )
        LzLogException("", "Cannot approximate local surface slopes using less than 2 samples!");

	//-------------------------------------------------
#pragma region "Create or use provided topology"
	//-------------------------------------------------

	MeshTopology lMeshTopo;
	if( ppMeshTopo )
		lMeshTopo = *ppMeshTopo;
	else
		lMeshTopo.Set( pMesh );
#pragma endregion


	//-------------------------------------------------
#pragma region "Find cracks at top using fat topology"
	//-------------------------------------------------

	List< VerAndIdx > lCrackVers;
	List< List<void *> > lCracks;
	{
        for( size_t c=0 ; c<lMeshTopo.Cracks().size() ; c++ )
		{
			// Crack edge
            const LzTriModel::TopEdge & lTopE = lMeshTopo.TopEdges()[ lMeshTopo.Cracks()[c] ];

			// Crack vertices
            const size_t lIdx0 = lTopE.mV[0];
            const size_t lIdx1 = lTopE.mV[1];

			// Keep or reject edge
			if( pExtPlane.SignedDistanceTo( pMesh.mVertices[lIdx0] ) > 0
			 && pExtPlane.SignedDistanceTo( pMesh.mVertices[lIdx1] ) > 0 )
			{
				// Map this crack to vertex indices
				List<void *> lCTVs;
				lCTVs.AddTail( lCrackVers.AddIntoIncrList(VerAndIdx(lIdx0), true, List<VerAndIdx>::AddIntoReturnPolicy::PosOfIdentical) );
				lCTVs.AddTail( lCrackVers.AddIntoIncrList(VerAndIdx(lIdx1), true, List<VerAndIdx>::AddIntoReturnPolicy::PosOfIdentical) );

				// Stash this crack
				lCracks.AddTail( lCTVs );
			}
		}

		// Update indices in VerAndIdx
        size_t lIdx = 0;
		BrowseList( iV, lCrackVers )
		{
			lCrackVers.GetAt(iV).mIdx = lIdx;
			lIdx++;
		}
	}
#pragma endregion


	//-------------------------------------------------
#pragma region "Project all crack vertices on least squares plane"
	//-------------------------------------------------

	Plane3D lExtPlane_LSQ;
	{
		// Compute weighted centroid on crack
		Point3D lCrackCenter;
		{
			Vector<Point3D> lMidCracks( lCracks.Count() );
			Vector<double> lCrackLengths( lCracks.Count() );
            size_t i = 0;
			BrowseList( iC, lCracks )
			{
				const List<void *> & lCrack = lCracks.GetAt( iC );

				const Point3D & lA = pMesh.mVertices[ lCrackVers.GetAt( lCrack.GetHead() ).mVer ];
				const Point3D & lB = pMesh.mVertices[ lCrackVers.GetAt( lCrack.GetTail() ).mVer ];

				lMidCracks[i] = lA.MidPointTo( lB );
				lCrackLengths[i] = lA.DistanceTo( lB );

				i++;
			}

			// Compute centroid
            lCrackCenter = LzMath::ToolBox::WeightedCentroid( lMidCracks, lCrackLengths );
		}

		// Compute least squares plane
		lExtPlane_LSQ = Plane3D( lCrackCenter, pExtPlane.Normal() );

		// Project vertices
		BrowseList( iV, lCrackVers )
		{
            const size_t lIdx = lCrackVers.GetAt(iV).mVer;
			pMesh.mVertices[ lIdx ] = lExtPlane_LSQ.Projection( pMesh.mVertices[ lIdx ] );
		}
	}
#pragma endregion


	//-------------------------------------------------
#pragma region "Create new vertices"
	//-------------------------------------------------

	Vector<Point3D> lNewVertices;
    for( size_t iL=0 ; iL<pLayersCount ; iL++ )
	{
		BrowseList( iV, lCrackVers )
		{
			const Point3D & lBaseVer = pMesh.mVertices[ lCrackVers.GetAt(iV).mVer ];
			lNewVertices.PushBack( lBaseVer + (iL+1)*pLayerHeight*pExtPlane.Normal() );
		}
	}
#pragma endregion


	// A few constants...
    const size_t lOldVerCount = pMesh.mVertices.size();
    const size_t lVersPerLayer = lCrackVers.Count();


	//-------------------------------------------------
#pragma region "Create new triangles"
	//-------------------------------------------------

    Vector<LzTriModel::Triangle> lNewTriangles;
    for( size_t iL=0 ; iL<pLayersCount ; iL++ )
	{
		BrowseList( iC, lCracks )
		{
			// Pointer to edge VerAndIdx
			const List<void *> & lCrackPos = lCracks.GetAt(iC);
			//
			const VerAndIdx & lHead_VI = lCrackVers.GetAt( lCrackPos.GetHead() );
			const VerAndIdx & lTail_VI = lCrackVers.GetAt( lCrackPos.GetTail() );

			// Find quad vertices
            size_t lCCWQuadIdx[4];
			if( iL == 0 )
			{
				lCCWQuadIdx[ 0 ] = lTail_VI.mVer;
				lCCWQuadIdx[ 1 ] = lHead_VI.mVer;

				lCCWQuadIdx[ 2 ] = lOldVerCount + lHead_VI.mIdx;
				lCCWQuadIdx[ 3 ] = lOldVerCount + lTail_VI.mIdx;
			}
			else
			{
				lCCWQuadIdx[ 0 ] = lOldVerCount + (iL-1)*lVersPerLayer + lTail_VI.mIdx;
				lCCWQuadIdx[ 1 ] = lOldVerCount + (iL-1)*lVersPerLayer + lHead_VI.mIdx;

				lCCWQuadIdx[ 2 ] = lOldVerCount + (iL+0)*lVersPerLayer + lHead_VI.mIdx;
				lCCWQuadIdx[ 3 ] = lOldVerCount + (iL+0)*lVersPerLayer + lTail_VI.mIdx;
			}

			// Create two triangles
            lNewTriangles.PushBack( LzTriModel::Triangle(lCCWQuadIdx[0],lCCWQuadIdx[1],lCCWQuadIdx[3],(size_t)0,(size_t)0,(size_t)0) );
            lNewTriangles.PushBack( LzTriModel::Triangle(lCCWQuadIdx[1],lCCWQuadIdx[2],lCCWQuadIdx[3],(size_t)0,(size_t)0,(size_t)0) );
		}
	}
#pragma endregion


	//-------------------------------------------------
#pragma region "Stash new vertices and triangles"
	//-------------------------------------------------

	// Vertices
    for( size_t nv=0 ; nv<lNewVertices.Size() ; nv++ )
		pMesh.mVertices.push_back( lNewVertices[nv] );

	// Triangles
    for( size_t nt=0 ; nt<lNewTriangles.Size() ; nt++ )
		pMesh.mTriangles.push_back( lNewTriangles[nt] );
#pragma endregion


	//-------------------------------------------------
	//
	// Here we have a new pMesh
	//
	//-------------------------------------------------

	// Update topology
	lMeshTopo.Set( pMesh );


	//-------------------------------------------------
#pragma region "Correct slope"
	//-------------------------------------------------

	{
		// Set outliner
		Outliner lNewOutliner;
		lNewOutliner.Set( pMesh, lMeshTopo );

		// Sample horizontal outlines: all in one
		List< List<Point3D> > lSampleHorOuts;
		Point3D lExtPlane_LSQ_Centroid;
        for( size_t s=0 ; s<pSamplesCount ; s++ )
		{
			// Sample plane
			Plane3D lSamplePlane = lExtPlane_LSQ;
			lSamplePlane -= s * pSampleHeight * lExtPlane_LSQ.Normal();

			// Compute outline
			List< List<Point3D> > lOuts;
			lNewOutliner.Cut( lSamplePlane, lOuts );

			// Check
			if( lOuts.Count() != 1 )
                LzLogException("", "Found 0 or more than one outline at sampling plane #"<<s<<"! "<<lOuts.Count()<<" outlines found.");

			// Check
			if( !Outliner::IsClosed( lOuts.GetHead() ) )
                LzLogException("", "Outline at extension plane is not closed!");

			// Stash
			lSampleHorOuts.Append( lOuts );

			// Compute centroid at least square extension plane
			if( s == 0 )
				lExtPlane_LSQ_Centroid = Outliner::OutlineWeightedCentroid( lOuts.GetHead() );
		}

		// Least squares lines for each crack vertex
		Vector< Line3D > lCrackVerLines( lCrackVers.Count() );
		{
            size_t cv = 0;
			BrowseList( iV, lCrackVers )
			{
				// Crack vertex and crack vertex plane
				const Point3D & lCrackVer = pMesh.mVertices[ lCrackVers.GetAt(iV).mVer ];
				const Plane3D lCrackCutPlane( lExtPlane_LSQ_Centroid, lExtPlane_LSQ.Normal(), lCrackVer-lExtPlane_LSQ_Centroid );

				// Cut sample outlines
				List< List<Point3D> > lInterList;
				Outliner::OutlineIntersections( lSampleHorOuts, lCrackCutPlane, lInterList );

				// Select inters
				List<Point3D> lSelectedInters;
				const Plane3D lSelector( lExtPlane_LSQ_Centroid, lCrackVer-lExtPlane_LSQ_Centroid );
				BrowseList( iIL, lInterList )
				BrowseList( iI, lInterList.GetAt(iIL) )
				{
					// Intersection point
					const Point3D & lInter = lInterList.GetAt(iIL).GetAt(iI);

					// Selected?
					if( lSelector.SignedDistanceTo(lInter) > 0 )
						lSelectedInters.AddTail( lInter );
				}

				// Compute least square line
				lCrackVerLines[cv] = Line3D( lSelectedInters );

				// Next
				cv++;
			}
		}

		// Displace all extended vertices to the corresponding least squares line
        for( size_t v=lOldVerCount ; v<pMesh.mVertices.size() ; v++ )
		{
			// Vertex index in crack ver index numeration
            size_t lCrackIdx = (v - lOldVerCount) % lVersPerLayer;

			// Vertex extension plane
			const Plane3D lVerPlane( pMesh.mVertices[v], lExtPlane_LSQ.Normal() );

			// New position
			const Point3D lNewVerPos = lCrackVerLines[lCrackIdx].Intersection( lVerPlane );

			// Move
			pMesh.mVertices[v] = lNewVerPos;
		}
	}
#pragma endregion


	//-------------------------------------------------
#pragma region "Laplacian relaxation"
	//-------------------------------------------------

	if( pLaplaceRelaxSteps )
	{
		// Fix previous skin nodes
        vector<bool> lFixedVer;
        lFixedVer.assign( pMesh.mVertices.size(), false );

		// Fix all old vertices below or on (original) extension plane
        for( size_t old_v=0 ; old_v<lOldVerCount ; old_v++ )
		{
			if( pExtPlane.SignedDistanceTo( pMesh.mVertices[old_v] ) <= 0 )
				lFixedVer[old_v] = true;
		}
#if 0
		Nope... ugly results

		// Fix all new vertices on top crack
        for( size_t crack_v=0 ; crack_v<lVersPerLayer ; crack_v++ )
			lFixedVer[pMesh.mVertices.size() - 1 - crack_v] = true;

		// Relax
		pMesh.MoveVerticesToMeanPos_HC( pLaplaceRelaxSteps, &lMeshTopo, &lFixedVer );
#else
		// Allow only sliding on new crack
        vector<const Plane3D *> lSlidePlanes;
        lSlidePlanes.assign( pMesh.mVertices.size(), (const Plane3D *)nullptr );
		//
        for( size_t crack_v=0 ; crack_v<lVersPerLayer ; crack_v++ )
			lSlidePlanes[pMesh.mVertices.size() - 1 - crack_v] = &lExtPlane_LSQ;

		// Relax
		pMesh.MoveVerticesToMeanPos_HC( pLaplaceRelaxSteps, &lMeshTopo, &lFixedVer, &lSlidePlanes );
#endif
	}
#pragma endregion
}

//================================================================================
void RepulseMesh( Mesh & pMesh, const Mesh & pRepulsor, double pRepDist, double pRelaxRadius )
{
	// Log
    LzLogN("", "RepulseMesh");

	// Process repulsor
	Mesh lRepulsor = pRepulsor;
	MeshTopology lRepTopo;
	lRepTopo.Set( lRepulsor );
	//
//**** NO NEED FOR SMOOTH NORMALIZATION : Normal extrude recomputes normals
//// Smooth normalize
//lRepulsor.SmoothNormalize( &lRepTopo );
//**** NO NEED FOR SMOOTH NORMALIZATION : Normal extrude recomputes normals
	//
	// Normal extrude
	lRepulsor.NormalExtrude( pRepDist, true );
//*** DEBUG
lRepulsor.Save( LzServices::StartUpPath() + "/_tmp_Rep_repulsor.obj");
//*** DEBUG
	//
    // Create FDC
	FastDistComp lRepFDC;
    lRepFDC.Set( lRepulsor, lRepTopo, LzTriModel::Patching::Split );

	// Repulse loop
    size_t iCount = 0;
const size_t lMaxCount = 20;
const size_t lRelaxSteps = 20;
	//
	// Moving vertices
    vector<bool> lIsVerFixed( pMesh.mVertices.size() );
	//
	// Mesh topology
	MeshTopology lMeshTopo;
	lMeshTopo.Set( pMesh );
	//
	// Projected indices
    List<size_t> lProjIndices;
	do
	{
#if 1
		// Find collision vers and repulse
		List<Point3D> lProjVertices;
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Check distance
			Point3D lProjV;
			Vector3D tmp;
            double lDist = lRepFDC.SignedDistance( Tree3D::NearestMode::Accurate,
                                                   pMesh.mVertices[v],
                                                   lProjV, tmp );
			if( lDist < -1e-6 )
			{
                // Store initial vertex as projected
				lProjVertices.AddTail( pMesh.mVertices[v] );
				lProjIndices.AddIntoIncrList( v, true );

				// Move vertex
				pMesh.mVertices[v] = lProjV;
			}
		}

		// Log
        LzLogM("", "Round "<<iCount<<": repulsed "<<lProjVertices.Count()<<" vertice(s).")

		// Check if any projected
		if( lProjVertices.Count() == 0 )
		{
            LzLogMsg("", "Finished repulsing mesh after "<<iCount<<" iteration(s).");
			return;
		}

// *** RELAX PROJECTED VERTICES ON REPULSOR??
// *** RELAX PROJECTED VERTICES ON REPULSOR??
// *** RELAX PROJECTED VERTICES ON REPULSOR??
// *** RELAX PROJECTED VERTICES ON REPULSOR??

		// Reset fixed flags
        std::fill( lIsVerFixed.begin(), lIsVerFixed.end(), false );
//		lIsVerFixed.SetAll( false );

		// Commit projected vertices so far, over ALL iterations
		BrowseList( iPI, lProjIndices )
			lIsVerFixed[ lProjIndices.GetAt(iPI) ] = true;

		// Fix vertices relax radius-away from projected vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Skip already fixed vertices
			if( lIsVerFixed[v] )
				continue;

			// Considered vertex
			const Point3D & lVer = pMesh.mVertices[v];

			// Find min distance to proj vers, in THIS iteration
			double lMinDist = std::numeric_limits<double>::max();
			BrowseList( iPV, lProjVertices )
			{
				double lDist = lVer.DistanceTo( lProjVertices.GetAt(iPV) );
				if( lMinDist > lDist )
					lMinDist = lDist;
			}

			// Check
			if( lMinDist > pRelaxRadius )
				lIsVerFixed[v] = true;
		}

		// Relax
        // PREVIOUSLY: pMesh.MoveVerticesToMeanPos( lRelaxSteps, &lMeshTopo, &lIsVerFixed );
        // But probably bad idea since this makes the mesh shrink!
        pMesh.MoveVerticesToMeanPos_HC( lRelaxSteps, &lMeshTopo, &lIsVerFixed );
#else
		// Find collision vers and repulse
		List<Point3D> lProjVertices;
		lIsVerFixed.SetAll( false );
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Check distance
			Point3D lProjV;
			double lDist = lRepFDC.SignedDistance( Tree3D::NearestMode::Accurate, pMesh.mVertices[v], lProjV, Vector3D() );
//if( lDist < 0 )
if( lDist < -1e-6 )
			{
				// Store initial vertex as projected
				lProjVertices.AddTail( pMesh.mVertices[v] );

				// Move vertex
				pMesh.mVertices[v] = lProjV;

				// Set as fixed
				lIsVerFixed[v] = true;
			}
		}

//*** DEBUG
{
	// Check
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		// Check distance
        double lDist = lRepFDC.AccurateSignedDistance( pMesh.mVertices[v] );
		if( lDist < -1e-6 )
            LzLogErr("", "Found dist < -1e-6 ! Dist= "<<lDist);
	}
}
//*** DEBUG

		// *** RELAX PROJECTED VERTICES ON REPULSOR??
		// *** RELAX PROJECTED VERTICES ON REPULSOR??
		// *** RELAX PROJECTED VERTICES ON REPULSOR??
		// *** RELAX PROJECTED VERTICES ON REPULSOR??

		// Check if any projected
		if( lProjVertices.Count() == 0 )
		{
            LzLogMsg("", "Finished repulsing mesh after "<<iCount<<" iteration(s).");
			return;
		}

//*** DEBUG
LzLogMsg("", "Projected "<<lProjVertices.Count()<<" vertices.");
//*** DEBUG

//*** DEBUG
pMesh.Save( LzServices::StartUpPath() + "/_tmp_Rep_after_proj_" + std::to_string(iCount) + ".obj");
//*** DEBUG

		// Fix vertices relax radius-away from projected vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			const Point3D & lVer = pMesh.mVertices[v];

			// Find min distance to proj vers
			double lMinDist = std::numeric_limits<double>::max();
			BrowseList( iPV, lProjVertices )
			{
				double lDist = lVer.DistanceTo( lProjVertices.GetAt(iPV) );
				if( lMinDist > lDist )
					lMinDist = lDist;
			}

			// Check
			if( lMinDist > pRelaxRadius )
				lIsVerFixed[v] = true;
		}

		// Relax
pMesh.MoveVerticesToMeanPos( lRelaxSteps, &lMeshTopo, &lIsVerFixed );
//pMesh.MoveVerticesToMeanPos_HC( lRelaxSteps, &lMeshTopo, &lIsVerFixed );

//*** DEBUG
{
	// Check
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		// Check distance
        double lDist = lRepFDC.AccurateSignedDistance( pMesh.mVertices[v] );
		if( lDist < -1e-6 )
            LzLogErr("", "AFTER RELAX Found dist < -1e-6 ! Dist= "<<lDist);
	}
}
//*** DEBUG

//*** DEBUG
pMesh.Save( LzServices::StartUpPath() + "/_tmp_Rep_after_relax_" + std::to_string(iCount) + ".obj");
//*** DEBUG
#endif
		// Next
		iCount++;
	}
	while( iCount < lMaxCount );

	// Error
    LzLogException("", "Could not repulse mesh after "<<iCount<<" iterations!");
}

//================================================================================
size_t SolveBilatConflict( Mesh & pA, Mesh & pB,
                           const SolveBilatConflictOptions & pOptions/*=SolveBilatConflictOptions()*/ )
{
    // Compute topologies
    MeshTopology lTopo_A;
    lTopo_A.Set( pA );
    //
    MeshTopology lTopo_B;
    lTopo_B.Set( pB );

    // FDCs
    FastDistComp lFDC_A;
    FastDistComp lFDC_B;

    //------------------------
    // Mesh relaxation
    //------------------------

    auto RelaxPostProjMesh = [&pOptions]( Mesh & pMesh,
                                          const MeshTopology & pMeshTopo,
                                          vector<bool> & pIsVerFixed,
                                          const List<size_t> & pAllProjIdx,
                                          const List<Point3D> & pProjVers )
        {
            // Reset fixed flags
            std::fill( pIsVerFixed.begin(), pIsVerFixed.end(), false );

            // Commit projected vertices so far, over ALL iterations
            BrowseList( iPI, pAllProjIdx )
                pIsVerFixed[ pAllProjIdx.GetAt(iPI) ] = true;

            // Fix vertices relax radius-away from projected vertices
            for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
            {
                // Skip already fixed vertices
                if( pIsVerFixed[v] )
                    continue;

                // Considered vertex
                const Point3D & lVer = pMesh.mVertices[v];

                // Find min distance to proj vers, in THIS iteration
                double lMinDist = std::numeric_limits<double>::max();
                BrowseList( iPV, pProjVers )
                {
                    double lDist = lVer.DistanceTo( pProjVers.GetAt(iPV) );
                    if( lMinDist > lDist )
                        lMinDist = lDist;
                }

                // Check
                if( lMinDist > pOptions.mRelaxRadius )
                    pIsVerFixed[v] = true;
            }

            // Relax
            if( pOptions.mHCSmoothSteps > 0 )
                pMesh.MoveVerticesToMeanPos_HC( pOptions.mHCSmoothSteps, &pMeshTopo, &pIsVerFixed );
        };

    //------------------------
    // Main loop
    //------------------------

    // Buffers for post-projection mesh relaxation
    vector<bool> lIsVerFixed_A;
    vector<bool> lIsVerFixed_B;
    if( pOptions.mHCSmoothSteps > 0 )
    {
        lIsVerFixed_A.resize( pA.mVertices.size() );
        lIsVerFixed_B.resize( pB.mVertices.size() );
    }

    // Projected indices for both meshes to avoid having to move them
    List<size_t> lAllProjIdx_A;
    List<size_t> lAllProjIdx_B;

    // Loop
    size_t iIter = 0;
    do
    {
        // Set or update FDCs
        lFDC_A.Set( pA, lTopo_A, pOptions.mPatching_A );
        lFDC_B.Set( pB, lTopo_B, pOptions.mPatching_B );

        // Check collision from A to B
        List<Point3D> lProjVers_A;
        for( size_t v=0 ; v<pA.mVertices.size() ; v++ )
        {
            // Considered vertex
            Point3D & lVer = pA.mVertices[v];

            // Check distance
            Point3D lProjVer;
            Vector3D tmp;
            const double lSgnDst = lFDC_B.SignedDistance( Tree3D::NearestMode::Accurate,
                                                          lVer, lProjVer, tmp );
            // Actual collision?
            if( lSgnDst < 0 )
            {
                // Store initial vertex as projected (store pre-projection position)
                lProjVers_A.AddTail( lVer );
                lAllProjIdx_A.AddIntoIncrList( v, true );

                // One-shot correction?
                if( lSgnDst > pOptions.mOneShotMinDst )
                {
                    // Move vertex all the way
                    lVer = lProjVer;
                }
                else
                {
                    // Move vertex at conflict solving rate
                    lVer += pOptions.mSolveRate*(lProjVer - lVer);
                }
            }
        }

        // Relax mesh A
        if( pOptions.mHCSmoothSteps > 0 )
            RelaxPostProjMesh( pA, lTopo_A, lIsVerFixed_A, lAllProjIdx_A, lProjVers_A );

        // Check collision from B to A
        List<Point3D> lProjVers_B;
        for( size_t v=0 ; v<pB.mVertices.size() ; v++ )
        {
            // Considered vertex
            Point3D & lVer = pB.mVertices[v];

            // Check distance
            Point3D lProjVer;
            Vector3D tmp;
            const double lSgnDst = lFDC_A.SignedDistance( Tree3D::NearestMode::Accurate,
                                                          lVer, lProjVer, tmp );
            // Actual collision?
            if( lSgnDst < 0 )
            {
                // Store initial vertex as projected (store pre-projection position)
                lProjVers_B.AddTail( lVer );
                lAllProjIdx_B.AddIntoIncrList( v, true );

                // One-shot correction?
                if( lSgnDst > pOptions.mOneShotMinDst )
                {
                    // Move vertex all the way
                    lVer = lProjVer;
                }
                else
                {
                    // Move vertex at conflict solving rate
                    lVer += pOptions.mSolveRate*(lProjVer - lVer);
                }
            }
        }

        // Relax mesh B
        if( pOptions.mHCSmoothSteps > 0 )
            RelaxPostProjMesh( pB, lTopo_B, lIsVerFixed_B, lAllProjIdx_B, lProjVers_B );

        // Debug save
        if( pOptions.mDebugSave )
        {
            const string lIter = std::to_string( iIter );
            pA.Save(LzServices::StartUpPath()+"/_DEBUG_SAVE__SolveBilatConflict__A__"+lIter+".obj");
            pB.Save(LzServices::StartUpPath()+"/_DEBUG_SAVE__SolveBilatConflict__B__"+lIter+".obj");
        }

        // End loop if no more projections
        if( lProjVers_A.Count() + lProjVers_B.Count() == 0 )
        {
            // Log
            LzLogM("", "Finished solving conflict at iteration "<<iIter<<".")
            return iIter;
        }

        // Increment iter
        iIter++;
    }
    while( iIter < pOptions.mMaxIter );

    // Error
    LzLogException("", "Could not solve bilateral conflict after "<<iIter<<" iterations!");
}

//================================================================================
void UNIT_TEST__SolveBilatConflict()
{
    // Log
    LzLogN("", "UNIT TEST: SolveBilatConflict")

    // Create ellipsoids
    Mesh lE0;
    GenEllipsoid( lE0, {0,0,0}, 10.0, 30.0, 30.0 );
//lE0.Load(R""(D:\LinkZ\Tech\_Code\Twinsight\Inferator\build-MSVC2019_64bit\release\_output_TibCart_Lat_CLOSED_CART__NO_CONFLICT.obj)"");
    //
    Mesh lE1;
    GenEllipsoid( lE1, {0,0,0}, 10.0, 30.0, 30.0 );
//lE1.Load(R""(D:\LinkZ\Tech\_Code\Twinsight\Inferator\build-MSVC2019_64bit\release\_output_FemCart_CLOSED_CART__NO_CONFLICT.obj)"");

    // Crash ellipsoids
    SolveBilatConflictOptions lOptions;
    lOptions.mDebugSave = true;  // Enable debug save
//lOptions.mHCSmoothSteps = 0;
#if 1
    // Small conflict
    lE1.RigidTransform( {0,0,0, Vector3D{17.0, 0.0, 0.0}} );
#elif 0
    // Big conflict
    lOptions.mSolveRate = 0.1;
    lOptions.mMaxIter = 100;
    //
    lE1.RigidTransform( {0,0,0, Vector3D{10.5, /*1*/0.0, 0.0}} );
#endif

    // Save
    lE0.Save(LzServices::StartUpPath()+"/_UNIT_TEST__SolveBilatConflict__E0.obj");
    lE1.Save(LzServices::StartUpPath()+"/_UNIT_TEST__SolveBilatConflict__E1.obj");

    // Solve bilateral conflict
    SolveBilatConflict( lE0, lE1, lOptions );

    // Save solution
    lE0.Save(LzServices::StartUpPath()+"/_UNIT_TEST__SolveBilatConflict__E0__solved.obj");
    lE1.Save(LzServices::StartUpPath()+"/_UNIT_TEST__SolveBilatConflict__E1__solved.obj");
}
#pragma endregion


#pragma region "Mesh refinement"
//================================================================================
//void RefineTriangle( RefineMode pMode, const Point3D & pV0, const Point3D & pV1, const Point3D & pV2, List<size_t> & pRefEdges )
void RefineTriangle( RefineMode pMode, double pThreshold, const Point3D pVers[3], List<size_t> & pRefEdges )
{
	if( pMode == RefineMode::AbsoluteLength )
	{
		//-----------------------------
		// AbsoluteLength
		//-----------------------------

		// Check all edges
		for( int v=0 ; v<3 ; v++ )
		{
			// Edge length
			double lLen = pVers[v].DistanceTo( pVers[(v+1) % 3] );

			// If edge violates criterion, stash
			if( lLen > pThreshold )
				pRefEdges.AddTail( (v+2) % 3 );
		}
	}
	else
	{
		//-----------------------------
		// RelativeLength
		//-----------------------------

		// Edge lengths
		double lEdgLens[3];

		// Find min edge length
        size_t lMinIdx;
		{
			// Init with edge 2
			lEdgLens[2] = pVers[0].DistanceTo( pVers[1] );
			lMinIdx = 2;

			// Edge 0 shorter?
			lEdgLens[0] = pVers[1].DistanceTo( pVers[2] );
			if( lEdgLens[lMinIdx] > lEdgLens[0] )
				lMinIdx = 0;

			// Edge 1 shorter?
			lEdgLens[1] = pVers[2].DistanceTo( pVers[0] );
			if( lEdgLens[lMinIdx] > lEdgLens[1] )
				lMinIdx = 1;
		}

		// Check next edge
		{
            size_t lIdx1 = (lMinIdx + 1) % 3;
			if( lEdgLens[lIdx1] > pThreshold * lEdgLens[lMinIdx] )
				pRefEdges.AddTail( lIdx1 );
		}

		// Check next next edge
		{
            size_t lIdx2 = (lMinIdx + 2) % 3;
			if( lEdgLens[lIdx2] > pThreshold * lEdgLens[lMinIdx] )
				pRefEdges.AddTail( lIdx2 );
		}
	}
}

//================================================================================
size_t RefineMesh( RefineMode pMode, double pThreshold, Mesh & pMesh )
{
	// Compute topology
	MeshTopology lTopo;
	lTopo.Set( pMesh );

	// Handy shortcut
	vector<Point3D> & lVertices = pMesh.mVertices;

	//----------------------------
	// Create new vertices
	//----------------------------

	// Init look-up table
	Vector<int> lTopEdgeToNewVers;
	lTopEdgeToNewVers.Assign( lTopo.TopEdges().size(), -1/* means: this edge is not split (yet?) */ );

	// Check all triangles
    for( size_t t=0 ; t<lTopo.TopTris().size() ; t++ )
	{
		// Top triangle
		const TopTri & lTTri = lTopo.TopTris()[t];

		// Check if edges need refinement
		const Point3D lVers[3] = { lVertices[ lTTri.mV[0] ], lVertices[ lTTri.mV[1] ], lVertices[ lTTri.mV[2] ] };
        List<size_t> lRefEdges;
		RefineTriangle( pMode, pThreshold, lVers, lRefEdges );

		// Refine all edges
		BrowseList( iRE, lRefEdges )
		{
			// Local and global edge indices
            size_t lLocalEdIdx = lRefEdges.GetAt( iRE );
            size_t lGlobEdIdx = lTTri.mE[ lLocalEdIdx ];

			// Not already subdivided?
			if( lTopEdgeToNewVers[ lGlobEdIdx ] == -1 )
			{
				// Create new vertex
				const Point3D lMidVer = lVers[ (lLocalEdIdx+1) % 3 ].MidPointTo( lVers[ (lLocalEdIdx+2) % 3 ] );

				// Update index in look-up
				lTopEdgeToNewVers[ lGlobEdIdx ] = lVertices.size();

				// Stash vertex
				lVertices.push_back( lMidVer );
			}
		}
	}

	//----------------------------
	// Create new triangles
	//----------------------------

	// Final result
    size_t lNbBadTris = 0;

	// New triangles
	vector<Triangle> lNewTriangles;

    for( size_t t=0 ; t<lTopo.TopTris().size() ; t++ )
	{
		// Top triangle
		const TopTri & lTTri = lTopo.TopTris()[t];

		// Count split edges
        List<size_t> lLocalSplitEdges;
        for( size_t e=0 ; e<3 ; e++ )
		{
			// Split edges are also stored in CCW order
			if( lTopEdgeToNewVers[ lTTri.mE[e] ] != -1 )
				lLocalSplitEdges.AddTail( e );
		}

		// Deal with cases
		switch( lLocalSplitEdges.Count() )
		{
		// No split
		case 0:
			{
				lNewTriangles.push_back( pMesh.mTriangles[t] );
			}
			break;

		// Single edge split
		case 1:
			{
				// Add 2 triangles
                size_t lLocalSplitEd = lLocalSplitEdges.GetHead();

				Triangle lTs[] =
				{
					Triangle( lTTri.mV[lLocalSplitEd], lTTri.mV[(lLocalSplitEd + 1)%3], lTopEdgeToNewVers[lTTri.mE[lLocalSplitEd]], 0, 0, 0 ),
					Triangle( lTopEdgeToNewVers[lTTri.mE[lLocalSplitEd]], lTTri.mV[(lLocalSplitEd + 2)%3], lTTri.mV[lLocalSplitEd], 0, 0, 0 )
				};

				for( int i=0 ; i<2 ; i++ )
				{
					// Stash new triangle
					lNewTriangles.push_back( lTs[i] );

					// Check quality
					const Point3D lVers[3] = { lVertices[ lTs[i].mIdxV[0] ], lVertices[ lTs[i].mIdxV[1] ], lVertices[ lTs[i].mIdxV[2] ] };
                    List<size_t> lRefEdges;
					RefineTriangle( pMode, pThreshold, lVers, lRefEdges );

					// Bad triangle?
					if( lRefEdges.Count() )
						lNbBadTris++;
				}
			}
			break;

		// Two edge splits
		case 2:
			{
				// Add 3 triangles
                size_t lNotEd = 3 - lLocalSplitEdges.GetHead() - lLocalSplitEdges.GetTail();

				Triangle lTs[] =
				{
					Triangle(lTTri.mV[(lNotEd + 1)%3], lTTri.mV[(lNotEd + 2)%3], lTopEdgeToNewVers[lTTri.mE[(lNotEd + 1)%3]], 0, 0, 0),
					Triangle(lTTri.mV[(lNotEd + 1)%3], lTopEdgeToNewVers[lTTri.mE[(lNotEd + 1)%3]], lTopEdgeToNewVers[lTTri.mE[(lNotEd + 2)%3]], 0, 0, 0),
					Triangle(lTopEdgeToNewVers[lTTri.mE[(lNotEd + 1)%3]], lTTri.mV[lNotEd], lTopEdgeToNewVers[lTTri.mE[(lNotEd + 2)%3]], 0, 0, 0)
				};

				for( int i=0 ; i<3 ; i++ )
				{
					// Stash new triangle
					lNewTriangles.push_back( lTs[i] );

					// Check quality
					const Point3D lVers[3] = { lVertices[ lTs[i].mIdxV[0] ], lVertices[ lTs[i].mIdxV[1] ], lVertices[ lTs[i].mIdxV[2] ] };
                    List<size_t> lRefEdges;
					RefineTriangle( pMode, pThreshold, lVers, lRefEdges );

					// Bad triangle?
					if( lRefEdges.Count() )
						lNbBadTris++;
				}
			}
			break;

		// Three edge splits
		case 3:
			{
				// Add 4 triangles
				Triangle lTs[] =
				{
					Triangle(lTTri.mV[0], lTopEdgeToNewVers[lTTri.mE[2]], lTopEdgeToNewVers[lTTri.mE[1]], 0, 0, 0),
					Triangle(lTTri.mV[1], lTopEdgeToNewVers[lTTri.mE[0]], lTopEdgeToNewVers[lTTri.mE[2]], 0, 0, 0),
					Triangle(lTTri.mV[2], lTopEdgeToNewVers[lTTri.mE[1]], lTopEdgeToNewVers[lTTri.mE[0]], 0, 0, 0),
					Triangle(lTopEdgeToNewVers[lTTri.mE[0]], lTopEdgeToNewVers[lTTri.mE[1]], lTopEdgeToNewVers[lTTri.mE[2]], 0, 0, 0)
				};

				for( int i=0 ; i<4 ; i++ )
				{
					// Stash new triangle
					lNewTriangles.push_back( lTs[i] );

					// Check quality
					const Point3D lVers[3] = { lVertices[ lTs[i].mIdxV[0] ], lVertices[ lTs[i].mIdxV[1] ], lVertices[ lTs[i].mIdxV[2] ] };
                    List<size_t> lRefEdges;
					RefineTriangle( pMode, pThreshold, lVers, lRefEdges );

					// Bad triangle?
					if( lRefEdges.Count() )
						lNbBadTris++;
				}
			}
			break;

		default:
            LzLogException("", "Holy smoking toledos! Found "<<lLocalSplitEdges.Count()<<" edges in a triangle...");
		}
	}

	// Commit new triangles to mesh
	pMesh.mTriangles = lNewTriangles;

	// Returns number of bad triangles
	return lNbBadTris;
}
#pragma endregion


#pragma region "Mesh simplification"
//================================================================================
#if 1 /// NEW ! better !!

// Allow or not flipped triangles (lorsque normale est calculable) : effet "cornet de frites" !!! ou 2 triangles sont en compet de flippage avec une ribambelle de petits edges en eventail qui vont etre collapses

//
//Preserve cracks, crack-neutral, remove cracks

// test anti triangle flipping !!!

// ajouter crack preservation au ver on edge collapse
//
//
//
//************ gerer cas interblocage triangles : ENCORE UN EDGE A VIRER---TRI 47.obj
//
void TopoClean( Mesh & pMesh, double pVerTol, double pEdgTol, bool pPreserveCracks, bool pCheckMEs )
{
	// Log
    LzLogN("", "LzTriModel::TopoClean: vertex tolerance= "<<pVerTol<<", edge tolerance= "<<pEdgTol<<".");

	//===============
	// Infinite loop
	//===============

    size_t iStep = 0;
	while( true )
	{
LzLogN("", "Iteration "<<iStep);

		//-------------------------------------
#pragma region "Reset iteration"
		//-------------------------------------

		// Reset topology
		MeshTopology lTopo;
        LzServices::DisableLog();
		{
			lTopo.Set( pMesh );
		}
        LzServices::EnableLog();
//**** DEBUG
if( pCheckMEs && lTopo.MultiEdges().size() )
{
    pMesh.Save(LzServices::StartUpPath() + "/_tmp_MEs_in_TopoClean.obj");
    LzLogException("", "Found MEs at iteration "<<iStep<<"!")
}
//**** DEBUG
		//
		// Reset impact zones
		vector<bool> lLostTopo;
		lLostTopo.assign( pMesh.mVertices.size(), false );

		// Reset rejected triangles
		vector<bool> lRejectTri;
		lRejectTri.assign( pMesh.mTriangles.size(), false );
        size_t lRejectTrisCount = 0;

		// Build cracked vers table
		vector<bool> lCrackedVer;
		if( pPreserveCracks )
		{
			// Create table of cracked vertices
			lCrackedVer.assign( pMesh.mVertices.size(), false );

			// Fill table
            for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
			{
				// Read edges stemming from v
                const List<size_t> & lEds = lTopo.TopVers()[ v ].mE;
				BrowseList( iE, lEds )
				{
					if( lTopo.TopEdges( lEds.GetAt(iE) ).mT.Count() == 1 )
					{
						lCrackedVer[ v ] = true;
						break;
					}
				}
			}
		}
#pragma endregion


		//-------------------------------------
#pragma region "Collapse vertices on vertices"
		//-------------------------------------

		// Reset vertex merge table and counter
        vector<size_t> lMergesWith;
		lMergesWith.assign( pMesh.mVertices.size(), pMesh.mVertices.size() ); // If idx = mVertices.size() then vertex is not merging
        size_t lMergers = 0;

		// Try to merge all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
//LzLogN("", "---- v= "<<v);

			// Skip vertices for which the topology has been lost in this iteration
			if( lLostTopo[v] )
				continue;

			// Determine if vertex needs to be merged, its merger and lost triangles

			// Read edges stemming from v
            const List<size_t> & lEds = lTopo.TopVers()[ v ].mE;

			// Skip floating vertices
			if( lEds.Count() == 0 )
				continue;

			// Skip if crack and need to skip cracks
			if( pPreserveCracks && lCrackedVer[v] )
				continue;

			// Find all potential collapsable edges
			class OppEdge
			{
			public:
                OppEdge( size_t pOpp, size_t pEdg ) : mOpp(pOpp), mEdg(pEdg) {}

                size_t mOpp;
                size_t mEdg;

                // To make it possible to use OppEdge within a LzServices::Sortable wrapper
                bool operator!=( const OppEdge & pOther ) const { return mOpp!=pOther.mOpp || mEdg!=pOther.mEdg; }
			};
            typedef LzServices::Sortable<OppEdge, double> PotColl;
			List<PotColl> lPotColls;
			//
			BrowseList( iE, lEds )
			{
				// Edge index
                const size_t lEdg = lEds.GetAt(iE);

				// Opposite ver idx
                const size_t lOpp = lTopo.TopEdges( lEdg ).OtherVer( v );

				// Skip if lost topology at lOpp
				if( lLostTopo[lOpp] )
					continue;

				// Skip if crack and need to skip cracks
				if( pPreserveCracks && lCrackedVer[lOpp] )
					continue;

				// Skip if lOpp < v
				if( lOpp < v )
					continue;

				// Compute distance and update min
				const double lLen = pMesh.mVertices[v].DistanceTo( pMesh.mVertices[lOpp] );

				// Skip if length above tolerance
				if( lLen >= pVerTol )
					continue;

				// Stash
				lPotColls.AddIntoIncrList( PotColl( OppEdge(lOpp, lEdg), lLen ), /*pUnique= */true );
			}

			// Skip vertex if no collapsable edges found
			if( lPotColls.Count() == 0 )
				continue;


			// And the Winner is....
			//
			// List of rejected triangles for chosen edge collapse
            List<size_t> lRej;
			//
			// Index of winning opposite vertex
            size_t w;
			//
			// List of triangles sharing V, and the winner W
            const List<size_t> & lTrisSharing_V = lTopo.TopVers( v ).mT;
            const List<size_t> * lpTrisSharing_W;


//**** DEBUG
// Check: rejected tri must share V or W
auto checkRejTri = [&]( size_t pTri, int pXX )
{
	const TopTri & lTri = lTopo.TopTris()[ pTri ];
	if( !lTri.HasVertex(v) && !lTri.HasVertex(w) )
        LzLogErr("", "******** CANNOT REJECT THIS TRI "<<pXX<<" !")
};
//**** DEBUG


			// Check all candidates
			BrowseList( iPC, lPotColls )
			{
				// Read candidate and reset rejected triangles
				lRej.DelAll();

				// Read index of opposite vertex
				w = lPotColls.GetAt(iPC).mT.mOpp;

				// Point tot list of triangles sharing w
				lpTrisSharing_W = &lTopo.TopVers( w ).mT;

				//---------------------------------------------------------
				// Find triangles sharing both V AND W
				//---------------------------------------------------------

				// ---> REJECT these triangles as they collapse
				//
				const List<int> & lTrisSharing_V_AND_W = lTopo.TopEdges( lPotColls.GetAt(iPC).mT.mEdg ).mT;
				BrowseList( iT, lTrisSharing_V_AND_W )
				{
					// Read triangle index
                    const size_t lTri = lTrisSharing_V_AND_W.GetAt( iT );

//*** DEBUG
if( lRejectTri[ lTri ] )
    LzLogException("", "Triangle already rejected 1 = "<<lTri<<" while processing vertex v= "<<v<<", w= "<<w<<" !")

//*** DEBUG
checkRejTri(lTri, 1);
					// Reject triangle
					lRej.AddTail( lTri );
				}

				//---------------------------------------------------------
				// Find triangles sharing V AND NOT W (and vice versa)
				//---------------------------------------------------------

				// Among tese triangles, we will look for annihilating triangles: X Y V and X Y W
				//
				// ---> REJECT these triangles as they collapse one onto another
				{
					// V but NOT W
					List<int> lTrisSharing_V_NOT_W;
					BrowseList( iT, lTrisSharing_V )
					{
                        const size_t lTri = lTrisSharing_V.GetAt(iT);
						if( lRej.FindElt(lTri) == nullptr )
							lTrisSharing_V_NOT_W.AddTail( lTri );
					}

					// W but NOT V
					List<int> lTrisSharing_W_NOT_V;
					BrowseList( iT, *lpTrisSharing_W )
					{
                        const size_t lTri = lpTrisSharing_W->GetAt(iT);
						if( lRej.FindElt(lTri) == nullptr )
							lTrisSharing_W_NOT_V.AddTail( lTri );
					}

					// Look for X-Y-V and X-Y-W (or Y-X-W) triangles
					BrowseList( iV, lTrisSharing_V_NOT_W )
					{
                        const size_t lIdx_V = lTrisSharing_V_NOT_W.GetAt( iV );
						const TopTri & lTri_V = lTopo.TopTris()[ lIdx_V ];

						BrowseList( iW, lTrisSharing_W_NOT_V )
						{
                            const size_t lIdx_W = lTrisSharing_W_NOT_V.GetAt( iW );
							const TopTri & lTri_W = lTopo.TopTris()[ lIdx_W ];

							// If the remaining 2 vertices are the same, reject both triangles
							if( lTri_V.EdgeWithoutVertex(v) == lTri_W.EdgeWithoutVertex(w) )
							{
//*** DEBUG
if( lRejectTri[ lIdx_V ] || lRejectTri[ lIdx_W ] )
    LzLogException("", "Triangle already rejected 2!")

//*** DEBUG
checkRejTri(lIdx_V, 2);
checkRejTri(lIdx_W, 3);

								// Reject triangles
								lRej.AddTail( lIdx_V );
								lRej.AddTail( lIdx_W );
							}
						}
					}
				}

				//---------------------------------------------------------
				// Check for possible resulting multiedges
				//---------------------------------------------------------

				// Find all vertices edge-connected with both V AND W
				//
				// Cannot use edge triangles here as X can be connected to V and W
				// through 2 different triangles i.e. triangle XVW might not exist
				//
                List<size_t> lVersCnx_V_AND_W;
				{
					// Vertices connected to V
                    List<size_t> lVersCnx_V;
					{
						// Edges from V
                        const List<size_t> & lEds = lTopo.TopVers()[ v ].mE;
						BrowseList( iE, lEds )
							lVersCnx_V.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( v ), /*pUnique= */true );
					}

					// Vertices connected to W
                    List<size_t> lVersCnx_W;
					{
						// Edges from w
                        const List<size_t> & lEds = lTopo.TopVers()[ w ].mE;
						BrowseList( iE, lEds )
							lVersCnx_W.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( w ), /*pUnique= */true );
					}

					// Vertices connected to V AND W
					lVersCnx_V_AND_W = lVersCnx_V.InterIncrLists( lVersCnx_W, /*pUniqueInInter= */true );
				}

//LzLogM("", "lVersCnx_V_AND_W= "<<lVersCnx_V_AND_W.ListToString())

				// Each vertex X in VersCnx_V_AND_W will produce a single edge X-V (== X-W) once V and W are merged
				// Check that the number of triangles sharing this new edge is not greater than 2
				// to avoid multi-edges
				//
				bool lCanCollapse = true;
				BrowseList( iX, lVersCnx_V_AND_W )
				{
					// X and edges from X
                    const size_t lX = lVersCnx_V_AND_W.GetAt(iX);
                    const List<size_t> & lEds = lTopo.TopVers()[ lX ].mE;

					// Find edge towards V
                    size_t lEd_X_to_V;
					BrowseList( iE, lEds )
					{
                        size_t lEd = lEds.GetAt( iE );
						if( lTopo.TopEdges()[ lEd ].OtherVer( lX ) == v )
						{
							lEd_X_to_V = lEd;
							goto Found_Ed_X_to_V;
						}
					}

					// Error
                    LzLogException("", "Could not find Ed_X_to_V!")

				Found_Ed_X_to_V:

    //LzLogM("", "Ed_X_to_V= "<<lTopo.TopEdges()[ lEd_X_to_V ].mV[0]<<" - "<<lTopo.TopEdges()[ lEd_X_to_V ].mV[1])
    //LzLogM("", "Ed_X_to_V: triangles= "<<lTopo.TopEdges()[ lEd_X_to_V ].mT.ListToString())

					// Find edge towards W
                    size_t lEd_X_to_W;
					BrowseList( iE, lEds )
					{
                        size_t lEd = lEds.GetAt( iE );
						if( lTopo.TopEdges()[ lEd ].OtherVer( lX ) == w )
						{
							lEd_X_to_W = lEd;
							goto Found_Ed_X_to_W;
						}
					}

					// Error
                    LzLogException("", "Could not find Ed_X_to_W!")

				Found_Ed_X_to_W:
					;

    //LzLogM("", "Ed_X_to_W= "<<lTopo.TopEdges()[ lEd_X_to_W ].mV[0]<<" - "<<lTopo.TopEdges()[ lEd_X_to_W ].mV[1])
    //LzLogM("", "Ed_X_to_W: triangles= "<<lTopo.TopEdges()[ lEd_X_to_W ].mT.ListToString())
    //LzLogM("", "lRej: triangles= "<<lRej.ListToString())

					// Now merge not discarded triangles from these two edges
                    List<size_t> lEdTris;
					{
						// Lambda
                        auto addTris = [&]( size_t e )
						{
							const List<int> & lTris = lTopo.TopEdges()[ e ].mT;
							BrowseList( iT, lTris )
							{
								const int lTri = lTris.GetAt(iT);
								if( lRej.FindElt(lTri) == nullptr )
									lEdTris.AddIntoIncrList( lTri, /*pUnique= */true );
							}
						};

						// Add not discarded triangles from edges X to W, and W to V
						addTris( lEd_X_to_V );
						addTris( lEd_X_to_W );
					}

					// And count them
					if( lEdTris.Count() > 2 )
					{
    //LzLogM("", "Collapsed edge tris= "<<lEdTris.Count()<<" if X= "<<lX)
    //LzLogM("", "Collapsed edge tris= "<<lEdTris.ListToString())
						lCanCollapse = false;
						break;
					}
				}

				// If this edge can be collapsed then goto collapse
				if( lCanCollapse )
					goto Collapse_current_V_W;
			}

			// Could not find any collapsable edge: consider next vertex
			continue;

		Collapse_current_V_W:

			//--------------------------------
			//===> Collapse edge
			//--------------------------------

			// Mark triangles as rejected
			BrowseList( iT, lRej )
			{
//if( lRej.GetAt(iT) == 3854 )
//	TxLogE("*********** rejecting 3854")

				lRejectTri[ lRej.GetAt(iT) ] = true;
			}

			// Update count
			lRejectTrisCount += lRej.Count();

//LzLogM("", "Rejecting triangles "<<lRej.ListToString())

			// Mark vertices as merged and update count
			lMergesWith[v] = w;
			lMergesWith[w] = v; //** even if merge is done while commiting v, need to flag w as merged
			lMergers++;

			// Find best merged position to preserve original surface
			{
				// Sum square volumes of local tets depending on merged vertex v'
				//
				// Minimize a L^2 + b L + c, for L in [0, 1]
				//
				// Then v' = v + L vw, i.e. vv' = L vw
				//
				double a = 0;
				double b = 0;

				// V, W, and VW
				const Point3D & V = pMesh.mVertices[v];
				const Point3D & W = pMesh.mVertices[w];
				const Vector3D VW = W - V;

				// Contribution from V
				BrowseList( iT, lTrisSharing_V )
				{
                    const size_t lTri = lTrisSharing_V.GetAt(iT);

					// Skip if rejected
					if( lRejectTri[ lTri ] )
						continue;

					// Edges
					const TopTri & lTT = lTopo.TopTris()[ lTri ];
					const int lE0 = lTT.EdgeWithVertex( v );
					const int lE1 = lTT.EdgeWithVertex( v, lE0 );

					// Points
					const Point3D & A = pMesh.mVertices[ lTopo.TopEdges()[lE0].OtherVer(v) ];
					const Point3D & B = pMesh.mVertices[ lTopo.TopEdges()[lE1].OtherVer(v) ];

					// Vol
					const double lVol = ((A - V)^(B - V)) * VW;

					// Accum
					a += lVol*lVol;
				}

				// Contribution from W
				BrowseList( iT, *lpTrisSharing_W )
				{
                    const size_t lTri = lpTrisSharing_W->GetAt(iT);

					// Skip if rejected
					if( lRejectTri[ lTri ] )
						continue;

					// Edges
					const TopTri & lTT = lTopo.TopTris()[ lTri ];
					const int lE0 = lTT.EdgeWithVertex( w );
					const int lE1 = lTT.EdgeWithVertex( w, lE0 );

					// Points
					const Point3D & A = pMesh.mVertices[ lTopo.TopEdges()[lE0].OtherVer(w) ];
					const Point3D & B = pMesh.mVertices[ lTopo.TopEdges()[lE1].OtherVer(w) ];

					// Vol
					const double lVol = ((A - W)^(B - W)) * VW;

					// Accum
					a += lVol*lVol;
					b += -2 * lVol * lVol;
				}

				// Minimize
				double L;
                if( LzMath::ToolBox::IsZero(2*a) )
					L = 1.0;
				else
					L = -b / (2*a);

				// Clamp
				if( L < 0.0 ) L = 0.0;
				if( L > 1.0 ) L = 1.0;

				// Set optimal position
				pMesh.mVertices[v] = V + L * VW;
pMesh.mVertices[w] = pMesh.mVertices[v]; //*** useless, since w > v: see remapping of indices below
			}

			// Mark vertices that lost their topo information
			{
				// List of all contaminated vertices
                List<size_t> lLost;

				// Lambda
                auto addVers = [&]( List<size_t> & pTo, size_t i )
				{
                    const List<size_t> & lEds = lTopo.TopVers()[ i ].mE;
					BrowseList( iE, lEds )
						pTo.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( i ), /*pUnique= */true );
				};

				// Add vertices that are edge connected to either v or w
				// v will add w, and w will add v
				addVers( lLost, v );
				addVers( lLost, w );

				// Flag all
				BrowseList( iL, lLost )
					lLostTopo[ lLost.GetAt(iL) ] = true;
			}
		}
#pragma endregion

		//
		// Finished iterating over all vertices for VER on VER collapse
		//
LzLogM("", "-----> Mergers= "<<lMergers)


		//-------------------------------------
#pragma region "Collapse vertices on edges"
		//-------------------------------------

		// Reset splitters and movers
		List<Triangle> lSplitTris;

		// Reset movers
		class Mover
		{
		public:
            Mover( size_t pIdx, const Point3D & pPos ) : mIdx(pIdx), mPos(pPos) {}

            size_t mIdx;
			Point3D mPos;
		};
		List<Mover> lMovers;

		// Try to project all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Skip if lost topology at v: we need the list of triangles sharing v
			if( lLostTopo[v] )
				continue;

			// Candidate vertex position
			const Point3D & lVer = pMesh.mVertices[ v ];

			// List candidate collapse edges
            List< LzServices::Sortable<unsigned, double> > lPotEdges;
			//
            const List<size_t> & lTrisAt_V = lTopo.TopVers( v ).mT;
			BrowseList( iT, lTrisAt_V )
			{
				// Candidate edge index
                const size_t lEdgIdx = lTopo.TopTris( lTrisAt_V.GetAt(iT) ).EdgeWithoutVertex( v );
				const TopEdge & lTEdg = lTopo.TopEdges( lEdgIdx );

				// Skip if any of the edge's vertices has lost topo
				if( lLostTopo[ lTEdg.mV[0] ] || lLostTopo[ lTEdg.mV[1] ] )
					continue;

				// Edge points
				const Point3D & A = pMesh.mVertices[ lTEdg.mV[0] ];
				const Point3D & B = pMesh.mVertices[ lTEdg.mV[1] ];

				// Check projection vs tolerance
				Line3D lEdgLine;
				try { lEdgLine = Line3D( A, B ); }
				catch( ... )
				{
					// Could not build line: next candidate
					continue;
				}

				// Acceptable distance is above or equal to edge tolerance
const Point3D lProjV = lEdgLine.Projection( lVer );
if( (lProjV - A)*(lProjV - B) >= 0 )
	continue;

				// Distance
				const double lDist = lProjV.DistanceTo( lVer );
				if( lDist >= pEdgTol )
					continue;

				// Store
                lPotEdges.AddIntoIncrList( LzServices::Sortable<unsigned, double>( lEdgIdx, lDist ), /*pUnique= */true );
			}

			// Skip vertex if no landing edges found
			if( lPotEdges.Count() == 0 )
				continue;

			// And the Winner is....
			//
			// List of rejected and split triangles for chosen edge collapse
            List<size_t> lRej;
			List<Triangle> lMySplitTris;
			//
			// Projected position of v
			Point3D lNewPos_V;

// List of all Apexes in split triangles to define Great Ball of Fire (see below)
//List<int> lSplitApexes;
const TopEdge * lpE; //****************************************************************** still needed to flag vertices in great ball of fire : see below

			// Read edges at v
            const List<size_t> & lEdgsAt_V = lTopo.TopVers( v ).mE;

			// Check all candidates
			BrowseList( iPE, lPotEdges )
			{
				// Reset winner edge information
				lRej.DelAll();
				lMySplitTris.DelAll();

				// Read candidate
const size_t e = lPotEdges.GetAt( iPE ).mT;
const double lDist = lPotEdges.GetAt( iPE ).mScore;
/*const TopEdge **/ lpE = &lTopo.TopEdges( e ); //***************************************************************

#if 1
////// Mid point of landing edge, because why not...
//////lNewPos_V = pMesh.mVertices[ lpE->mV[0] ].MidPointTo( pMesh.mVertices[ lpE->mV[1] ] );

					// New position = projection on edge (falls between edge vertices, as checked previously)
//************ keep this information within Sortable< > ...?
					lNewPos_V = Line3D( pMesh.mVertices[ lpE->mV[0] ], pMesh.mVertices[ lpE->mV[1] ] ).Projection( lVer );
#else
					// ***** MINIMIZATION TETRAS
					// ***** MINIMIZATION TETRAS
					// ***** MINIMIZATION TETRAS
#endif
				//
				// NB: if the only triangle around the edge is the one having v
				//     then this triangle will simply be rejected, and v will be
				//     projected somewhere on its opposite edge.
				//

				// Triangles around the edge
				const List<int> & lEdgeTris = lpE->mT;
				BrowseList( iT, lEdgeTris )
				{
//if( LOGGG )
//LzLogM("", "V-E: v= "<<v<<", e= "<<lpE->mV[0]<<" - "<<lpE->mV[1]<<", e.t= "<<lpE->mT.ListToString())

					// Read triangle
					const int lTriIdx = lEdgeTris.GetAt( iT );
					const TopTri & lTT = lTopo.TopTris( lTriIdx );

					// Reject triangle
					lRej.AddTail( lTriIdx );

					// Skip triangle if collapsed triangle
					if( lTT.HasVertex( v ) )
						continue;

					// Find apex
					//
					// ApexIdx:    tri index (0, 1 or 2) of vertex opposing edge E in this triangle
					// ApexVerIdx: index of apex in mVertices
					//
					// If apex is v, then this collapse will produce a multi-edge
					int lApexIdx;
                    const size_t lApexVerIdx = lTT.VertexOppositeToEdge(e, &lApexIdx);

					// Check if that apex already has an edge to v
					BrowseList( iE, lEdgsAt_V )
					{
						if( lTopo.TopEdges( lEdgsAt_V.GetAt(iE) ).OtherVer( v ) == lApexVerIdx )
						{
//if( LOGGG )
//LzLogM("", "V-E: REJECTED")
							goto Consider_next_candidate_edge;
						}
					}

					// Split into 2 triangles
					const Triangle lSpTris[] =
					{
                        Triangle( v, (size_t)lTT.mV[lApexIdx], (size_t)lTT.mV[(lApexIdx + 1)%3], (size_t)0, (size_t)0, (size_t)0 ),
                        Triangle( v, (size_t)lTT.mV[(lApexIdx + 2)%3], (size_t)lTT.mV[lApexIdx], (size_t)0, (size_t)0, (size_t)0 )
					};

					// Check heights in triangles split beyond e using new position for v
					//
					// All heights in new triangles must be above threshold
					for( int t=0 ; t<2 ; t++ )
					{
						// Triangle
						const Triangle & T = lSpTris[ t ];

//if( LOGGG )
//LzLogM("", "V-E: Checking split triangle :"<<v<<", "<<T.mIdxV[1]<<", "<<T.mIdxV[2])

						// Triangle points
						const Point3D * ABC[] = { &lNewPos_V, &pMesh.mVertices[ T.mIdxV[1] ], &pMesh.mVertices[ T.mIdxV[2] ] };

						// Check all three heights
						for( int h=0 ; h<3 ; h++ )
						{
							const Point3D & A = *ABC[  h        ];
							const Point3D & B = *ABC[ (h + 1)%3 ];
							const Point3D & C = *ABC[ (h + 2)%3 ];

							// Build edge line
							Line3D BC;
							try { BC = Line3D( B, C ); }
							catch( ... )
							{
								goto Consider_next_candidate_edge;
							}

							// Check if new height can only make the overall situation better
							// to avoid infinite loops
							if( BC.DistanceTo( A ) <= lDist )
								goto Consider_next_candidate_edge;
						}
					}

					// Check heights of triangles behind mover vertex v
					//
					// Consider all triangles sharing v that are not collapsed
					BrowseList( iT, lTrisAt_V )
					{
                        const size_t lIdx = lTrisAt_V.GetAt( iT );
						const TopTri & lTT = lTopo.TopTris( lIdx );

// Skip collapsed
if( lTT.HasVertex( v ) && lTT.HasVertex( lpE->mV[0] ) && lTT.HasVertex( lpE->mV[1] ) ) //********************************************** retenir index
	continue;

						// Triangle points
						const Point3D * ABC[3];
						for( int i=0 ; i<3 ; i++ )
						{
							if( lTT.mV[i] == v )
								ABC[ i ] = &lNewPos_V;
							else
								ABC[ i ] = &pMesh.mVertices[ lTT.mV[i] ];
						}

						// Check all three heights
						for( int h=0 ; h<3 ; h++ )
						{
							const Point3D & A = *ABC[  h        ];
							const Point3D & B = *ABC[ (h + 1)%3 ];
							const Point3D & C = *ABC[ (h + 2)%3 ];

							// Build edge line
							Line3D BC;
							try { BC = Line3D( B, C ); }
							catch( ... )
							{
								goto Consider_next_candidate_edge;
							}


							// Check if new height can only make the overall situation better
							// to avoid infinite loops
							if( BC.DistanceTo( A ) <= lDist )
								goto Consider_next_candidate_edge;
						}
					}

					// Stash
					for( int t=0 ; t<2 ; t++ )
						lMySplitTris.AddTail( lSpTris[ t ] );
				}

				//
				//--- Additional checks here --- (cracks preservation etc.)
				//

				// Vertex v can be projected on this edge
				goto Project_current_V_on_E;

			Consider_next_candidate_edge:
				;
			}

//if( LOGGG )
//LzLogM("", "Next vertex...............")

			// Could not find any edge to project onto: consider next vertex
			continue;

Project_current_V_on_E:

//if( LOGGG )
LzLogM("", "Projecting current "<<v) //<<" on edge "<<lpE->mV[0]<<" - "<<lpE->mV[1])
//pMesh.Save(LzServices::StartUpPath() + "/__tmp_MEs_in_TopoClean__v= "+std::to_string(v)+".obj");

			//--------------------------------
			//===> Project vertex on edge
			//--------------------------------

			// Build split tris
			lSplitTris.Append( lMySplitTris );




////**** DEBUG
//LzLogM("", "Rejected triangles= "<<lRej.ListToString())
//LzLogM("", "Split apexes   = "<<lSplitApexes.ListToString())
////**** DEBUG


			// Mark triangles as rejected
			BrowseList( iT, lRej )
				lRejectTri[ lRej.GetAt(iT) ] = true;

			// Update count
			lRejectTrisCount += lRej.Count();

			// Mark vertex as mover
#if 1
			//// Mid point of landing edge, because why not...
			//const Point3D lNewPos = pMesh.mVertices[ lpE->mV[0] ].MidPointTo( pMesh.mVertices[ lpE->mV[1] ] );

			lMovers.AddTail( Mover(v, lNewPos_V ) );
#else
			// Find best merged position to preserve original surface
			{
				// Consider tetraherons! Consider its tree!


				lMovers.AddTail( Mover(v, some optimal 3D point ) );
			}
#endif

			// Mark vertices that lost their topo information
			{
				// Conjure Great Ball of Fire around rejected and split triangles
#if 0
				// Flag vertex
				lLostTopo[ v ] = true;

				// Flag edge ends
				lLostTopo[ lpE->mV[0] ] = true;
				lLostTopo[ lpE->mV[1] ] = true;

				// Flag all apexes
				BrowseList( iA, lSplitApexes )
					lLostTopo[ lSplitApexes.GetAt(iA) ] = true;
#else
				// List of all contaminated vertices
                List<size_t> lLost;

				// Lambda
                auto addVers = [&]( List<size_t> & pTo, size_t i )
				{
                    const List<size_t> & lEds = lTopo.TopVers()[ i ].mE;
					BrowseList( iE, lEds )
						pTo.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( i ), /*pUnique= */true );
				};

				// Add vertices that are edge connected to either v or w
				// v will add w, and w will add v
				addVers( lLost, v );
				addVers( lLost, lpE->mV[0] );
				addVers( lLost, lpE->mV[0] );

				// Flag all
				BrowseList( iL, lLost )
					lLostTopo[ lLost.GetAt(iL) ] = true;
#endif
			}
		}
#pragma endregion

		//
		// Finished iterating over all vertices for VER on EDGE collapse
		//
{
    LzLogN("", "lSplitTris= "<<lSplitTris.Count());
	BrowseList( iT, lSplitTris )
        LzLogM("", lSplitTris.GetAt(iT).mIdxV[0]<<" - "<<lSplitTris.GetAt(iT).mIdxV[1]<<" - "<<lSplitTris.GetAt(iT).mIdxV[2])
}

string lDEBUG = "";
static size_t DEBUG_I = 0;
{
    LzLogN("", "lMovers= "<<lMovers.Count());
	BrowseList( iM, lMovers )
        LzLogM("", lMovers.GetAt(iM).mIdx);

	if( lMovers.Count() == 1 )
	{
		if( lMovers[0].mIdx == 2055 )
		{
			lDEBUG = "After move 2055 - "+std::to_string( DEBUG_I );
            LzLogM("", "Saving result to "<<lDEBUG)
			DEBUG_I++;
		}
		else
		if( lMovers[0].mIdx == 2004 )
		{
			lDEBUG = "After move 2004 - "+std::to_string( DEBUG_I );
            LzLogM("", "Saving result to "<<lDEBUG)
			DEBUG_I++;
		}
	}
}


		//-------------------------------------
		// Finish if not merged anything
		//-------------------------------------

		if( lMergers==0 && lSplitTris.Count()==0 && lMovers.Count()==0 )
			break;


		//-------------------------------------
#pragma region "Commit"
		//-------------------------------------

		// New vertices and old-to-new mapping
		vector<Point3D> lNewVertices( pMesh.mVertices.size() - lMergers );
LzLogM("", "-----> New vertices= "<<lNewVertices.size())

        vector<size_t> lOldToNewVerIdx( pMesh.mVertices.size() );
		{
            size_t iNewIdx = 0;
            for( size_t old_idx=0 ; old_idx<lMergesWith.size() ; old_idx++ )
			{
				// Read merger
                const size_t w = lMergesWith[ old_idx ];

				// Vertex is merging?
				if( w < pMesh.mVertices.size() )
				{
					// Yes, merging

					// Do once
					if( old_idx < w )
					{
						lNewVertices[ iNewIdx ] = pMesh.mVertices[ old_idx ];
						lOldToNewVerIdx[ old_idx ] = iNewIdx;
						lOldToNewVerIdx[       w ] = iNewIdx;
						iNewIdx++;
					}
				}
				else
				{
					// Nope, not merging

					// Copy previous position
					lNewVertices[ iNewIdx ] = pMesh.mVertices[ old_idx ];
					lOldToNewVerIdx[ old_idx ] = iNewIdx;
					iNewIdx++;
				}
			}
		}

		//--------------------------------
		// Normals are lost
		//--------------------------------

		pMesh.mNormals = { Vector3D(1, 0, 0) };

		//--------------------------------
		// Create new triangles
		//--------------------------------

		// Create new array
		vector<Triangle> lNewTriangles( pMesh.mTriangles.size() - lRejectTrisCount + lSplitTris.Count() );
		{
            size_t iIdx = 0;

			// Renum old triangles, skip rejected
            for( size_t old=0 ; old<pMesh.mTriangles.size() ; old++ )
			{
				// Keeping this triangle?
				if( !lRejectTri[ old ] )
				{
					// Create new triangle
					const Triangle & lOld = pMesh.mTriangles[ old ];
					lNewTriangles[ iIdx ] = Triangle( lOldToNewVerIdx[ lOld.mIdxV[0] ],
													  lOldToNewVerIdx[ lOld.mIdxV[1] ],
                                                      lOldToNewVerIdx[ lOld.mIdxV[2] ], (size_t)0, (size_t)0, (size_t)0 ); // Pointing to the unique normal

// Check
if( lNewTriangles[ iIdx ].HasRepeatedVers() )
    LzLogException("", "Found a triangle with repeated vertices! Tri= "<<lNewTriangles[ iIdx ].ToString())

					// Next
					iIdx++;
				}
			}

			// Add split tris from ver-to-edge collapse
			BrowseList( iST, lSplitTris )
			{
				// Create split triangle
				const Triangle & lSpT = lSplitTris.GetAt( iST );
				lNewTriangles[ iIdx ] = Triangle( lOldToNewVerIdx[ lSpT.mIdxV[0] ],
												  lOldToNewVerIdx[ lSpT.mIdxV[1] ],
                                                  lOldToNewVerIdx[ lSpT.mIdxV[2] ], (size_t)0, (size_t)0, (size_t)0 ); // Pointing to the unique normal
// Check
if( lNewTriangles[ iIdx ].HasRepeatedVers() )
    LzLogException("", "Found a triangle with repeated vertices! Tri= "<<lNewTriangles[ iIdx ].ToString())

				// Next
				iIdx++;
			}

			// Move movers
			BrowseList( iM, lMovers )
			{
				// Get mover
				const Mover & lMov = lMovers.GetAt( iM );

				// Move vertex position
				lNewVertices[ lOldToNewVerIdx[ lMov.mIdx ] ] = lMov.mPos;
//***** DEBUG
// check that none of the movers is collapsing onto another ver...
//***** DEBUG

			}
		}

		//--------------------------------
		// Commit
		//--------------------------------

		pMesh.mVertices = lNewVertices;
		pMesh.mTriangles = lNewTriangles;
#pragma endregion


//***** DEBUG
if( lDEBUG.length() )
{
	pMesh.SmoothNormalize();
    pMesh.Save(LzServices::StartUpPath()+"/_tmp_TopoClean_DEBUG_"+lDEBUG+".obj");
}
//{
//	pMesh.SmoothNormalize();
//	pMesh.Save(LzServices::StartUpPath()+"/_tmp_TopoClean_it_"+std::to_string(iStep)+".obj");
//}
//***** DEBUG


		// Log
        LzLogM("", "Mesh::TopoClean: iteration "<<iStep<<": vertices= "<<pMesh.mVertices.size()<<", triangles= "<<pMesh.mTriangles.size()<<".")
		iStep++;
	}

    // Need to recompute DL
	pMesh.ResetDisplayList();
}
#else
/// OLD..
void TopoClean( Mesh & pMesh, double pVerTol, double /*pEdgTol*/, bool pPreserveCracks )
{
	// Log
    LzLogN("", "LzTriModel::TopoClean: vertex tolerance= "<<pVerTol<<", edge tolerance= "<<"UNUSED"<<".");

	//===============
	// Infinite loop
	//===============

    size_t iStep = 0;
	while( true )
	{
LzLogN("", "Iteration "<<iStep);

		//--------------------------------
		// Reset iteration
		//--------------------------------

		// Reset topology
		MeshTopology lTopo;
        LzServices::Log::Disable();
		{
			lTopo.Set( pMesh );
		}
        LzServices::Log::Enable();
		//
		vector<bool> lLostTopo;
		lLostTopo.assign( pMesh.mVertices.size(), false );

		// Reset vertex merge table
        vector<size_t> lMergesWith;
		lMergesWith.assign( pMesh.mVertices.size(), pMesh.mVertices.size() ); // If idx = mVertices.size() then vertex is not merging

		// Reset rejected triangles
size_t lLostTris = 0;
		vector<bool> lRejectTri;
		lRejectTri.assign( pMesh.mTriangles.size(), false );

		// Build cracked vers table
		vector<bool> lCrackedVer;
		if( pPreserveCracks )
		{
			// Create table of cracked vertices
			lCrackedVer.assign( pMesh.mVertices.size(), false );

			// Fill table
            for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
			{
				// Read edges stemming from v
                const List<size_t> & lEds = lTopo.TopVers()[ v ].mE;
				BrowseList( iE, lEds )
				{
					if( lTopo.TopEdges( lEds.GetAt(iE) ).mT.Count() == 1 )
					{
						lCrackedVer[ v ] = true;
						break;
					}
				}
			}
		}


		//-------------------------------------
#pragma region "Collapse vertices on vertices"
		//-------------------------------------

		// Reset mergers
        size_t lMergers = 0;

		// Try to merge all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
//LzLogN("", "---- v= "<<v);

			// Skip vertices for which the topology has been lost in this iteration
			if( lLostTopo[v] )
				continue;

			// Determine if vertex needs to be merged, its merger and lost triangles

			// Read edges stemming from v
            const List<size_t> & lEds = lTopo.TopVers()[ v ].mE;

			// Skip floating vertices
			if( lEds.Count() == 0 )
				continue;

			// Skip if crack and need to skip cracks
			if( pPreserveCracks && lCrackedVer[v] )
				continue;

			// Find smallest edge-connected neighbor
			double lMinLen = std::numeric_limits<double>::max();
size_t w; //****** a renommer en W et descendre
size_t lMinEdg; //****** inutile, a virer

//List< Sortable<struct{ size_t x, size_t y }, double> > de candidats lOpp potentiels ==> appeles pot_w
// tries par longueur d'edge du plus petit au plus grand

			BrowseList( iE, lEds )
			{
				// Edge index
                const size_t lEdg = lEds.GetAt(iE);

				// Opposite ver idx
                const size_t lOpp = lTopo.TopEdges( lEdg ).OtherVer( v );

				// Skip if lost topology at lOpp
				if( lLostTopo[lOpp] )
					continue;

				// Skip if crack and need to skip cracks
				if( pPreserveCracks && lCrackedVer[lOpp] )
					continue;

				// Skip if lOpp < v
				if( lOpp < v )
					continue;

				// Compute distance and update min
				const double lLen = pMesh.mVertices[v].DistanceTo( pMesh.mVertices[lOpp] );
				if( lMinLen > lLen )
				{
					lMinLen = lLen;
					w = lOpp;
					lMinEdg = lEdg;
				}
			}



// Check tolerance ***** to be moved inside the loop
if( lMinLen >= pVerTol )
	continue;



//**** DEBUG
// Check: rejected tri must share V or W
auto checkRejTri = [&]( size_t pTri, int pXX )
{
	const TopTri & lTri = lTopo.TopTris()[ pTri ];
	if( !lTri.HasVertex(v) && !lTri.HasVertex(w) )
        LzLogErr("", "******** CANNOT REJECT THIS TRI "<<pXX<<" !")
};
//**** DEBUG


//***********************************************************************
			// List of rejected triangles
            List<size_t> lRej;

// BrowseList( iPW, l... )
//{


			// List of triangles sharing collapsed vertices, but not necessarily the whole collapsed edge
            const List<size_t> & lTrisSharing_V = lTopo.TopVers()[ v ].mT;
            const List<size_t> & lTrisSharing_W = lTopo.TopVers()[ w ].mT;

// Update list of rejected triangles for current candidate
//lRej.DelAll();
int x = 666; //*** a virer
			{
				//---------------------------------------------------------
				// Find triangles sharing both V AND W
				//---------------------------------------------------------

				// ---> REJECT these triangles as they collapse
				//
				BrowseList( iT, lTrisSharing_V )
				{
					// Read triangle index
                    const size_t lTri = lTrisSharing_V.GetAt( iT );

					// Dismiss triangle if shares both V AND W
					if( lTopo.TopTris()[ lTri ].HasVertex( w ) )
					{
//*** DEBUG
if( lRejectTri[ lTri ] )
    LzLogException("", "Triangle already rejected 1 = "<<lTri<<" while processing vertex v= "<<v<<", w= "<<w<<" !")

//*** DEBUG
checkRejTri(lTri, 1);
						// Reject triangle
						lRej.AddTail( lTri );
						lLostTris++;
					}
				}

				//---------------------------------------------------------
				// Find triangles sharing V AND NOT W (and vice versa)
				//---------------------------------------------------------

				// Among tese triangles, we will look for annihilating triangles: X Y V and X Y W
				//
				// ---> REJECT these triangles as they collapse one onto another
				{
					// V but NOT W
					List<int> lTrisSharing_V_NOT_W;
					BrowseList( iT, lTrisSharing_V )
					{
                        const size_t lTri = lTrisSharing_V.GetAt(iT);
						if( lRej.FindElt(lTri) == nullptr )
							lTrisSharing_V_NOT_W.AddTail( lTri );
					}

					// W but NOT V
					List<int> lTrisSharing_W_NOT_V;
					BrowseList( iT, lTrisSharing_W )
					{
                        const size_t lTri = lTrisSharing_W.GetAt(iT);
						if( lRej.FindElt(lTri) == nullptr )
							lTrisSharing_W_NOT_V.AddTail( lTri );
					}

					// Look for X-Y-V and X-Y-W (or Y-X-W) triangles
					BrowseList( iV, lTrisSharing_V_NOT_W )
					{
                        const size_t lIdx_V = lTrisSharing_V_NOT_W.GetAt( iV );
						const TopTri & lTri_V = lTopo.TopTris()[ lIdx_V ];

						BrowseList( iW, lTrisSharing_W_NOT_V )
						{
                            const size_t lIdx_W = lTrisSharing_W_NOT_V.GetAt( iW );
							const TopTri & lTri_W = lTopo.TopTris()[ lIdx_W ];

							// If the remaining 2 vertices are the same, reject both triangles
							if( lTri_V.EdgeWithoutVertex(v) == lTri_W.EdgeWithoutVertex(w) )
							{
//*** DEBUG
if( lRejectTri[ lIdx_V ] || lRejectTri[ lIdx_W ] )
    LzLogException("", "Triangle already rejected 2!")

//*** DEBUG
checkRejTri(lIdx_V, 2);
checkRejTri(lIdx_W, 3);

								// Reject triangles
								lRej.AddTail( lIdx_V );
								lRej.AddTail( lIdx_W );
								lLostTris += 2;
							}
						}
					}
				}
			}

			//---------------------------------------------------------
			// Check for possible resulting multiedges
			//---------------------------------------------------------

			// Find all vertices edge-connected with both V AND W
            List<size_t> lVersCnx_V_AND_W;
			{
				// Vertices connected to V
                List<size_t> lVersCnx_V;
				{
					// Edges from V
                    const List<size_t> & lEds = lTopo.TopVers()[ v ].mE;
					BrowseList( iE, lEds )
						lVersCnx_V.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( v ), /*pUnique= */true );
				}

				// Vertices connected to W
                List<size_t> lVersCnx_W;
				{
					// Edges from w
                    const List<size_t> & lEds = lTopo.TopVers()[ w ].mE;
					BrowseList( iE, lEds )
						lVersCnx_W.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( w ), /*pUnique= */true );
				}

				// Vertices connected to V AND W
				lVersCnx_V_AND_W = lVersCnx_V.InterIncrLists( lVersCnx_W, /*pUniqueInInter= */true );
			}

//LzLogM("", "lVersCnx_V_AND_W= "<<lVersCnx_V_AND_W.ListToString())

			// Each vertex X in VersCnx_V_AND_W will produce a single edge X-V (== X-W) once V and W are merged
			// Check that the number of triangles sharing this new edge is not greater than 2
			// to avoid multi-edges
			//
			bool lCanCollapse = true;
			BrowseList( iX, lVersCnx_V_AND_W )
			{
				// X and edges from X
                const size_t lX = lVersCnx_V_AND_W.GetAt(iX);
                const List<size_t> & lEds = lTopo.TopVers()[ lX ].mE;

				// Find edge towards V
                size_t lEd_X_to_V;
				BrowseList( iE, lEds )
				{
                    size_t lEd = lEds.GetAt( iE );
					if( lTopo.TopEdges()[ lEd ].OtherVer( lX ) == v )
					{
						lEd_X_to_V = lEd;
						goto Found_Ed_X_to_V;
					}
				}

				// Error
                LzLogException("", "Could not find Ed_X_to_V!")

			Found_Ed_X_to_V:

//LzLogM("", "Ed_X_to_V= "<<lTopo.TopEdges()[ lEd_X_to_V ].mV[0]<<" - "<<lTopo.TopEdges()[ lEd_X_to_V ].mV[1])
//LzLogM("", "Ed_X_to_V: triangles= "<<lTopo.TopEdges()[ lEd_X_to_V ].mT.ListToString())

				// Find edge towards W
                size_t lEd_X_to_W;
				BrowseList( iE, lEds )
				{
                    size_t lEd = lEds.GetAt( iE );
					if( lTopo.TopEdges()[ lEd ].OtherVer( lX ) == w )
					{
						lEd_X_to_W = lEd;
						goto Found_Ed_X_to_W;
					}
				}

				// Error
                LzLogException("", "Could not find Ed_X_to_W!")

			Found_Ed_X_to_W:
				;

//LzLogM("", "Ed_X_to_W= "<<lTopo.TopEdges()[ lEd_X_to_W ].mV[0]<<" - "<<lTopo.TopEdges()[ lEd_X_to_W ].mV[1])
//LzLogM("", "Ed_X_to_W: triangles= "<<lTopo.TopEdges()[ lEd_X_to_W ].mT.ListToString())
//LzLogM("", "lRej: triangles= "<<lRej.ListToString())

				// Now merge not discarded triangles from these two edges
                List<size_t> lEdTris;
				{
					// Lambda
                    auto addTris = [&]( size_t e )
					{
						const List<int> & lTris = lTopo.TopEdges()[ e ].mT;
						BrowseList( iT, lTris )
						{
							const int lTri = lTris.GetAt(iT);
							if( lRej.FindElt(lTri) == nullptr )
								lEdTris.AddIntoIncrList( lTri, /*pUnique= */true );
						}
					};

					// Add not discarded triangles from edges X to W, and W to V
					addTris( lEd_X_to_V );
					addTris( lEd_X_to_W );
				}

				// And count them
				if( lEdTris.Count() > 2 )
				{
//LzLogM("", "Collapsed edge tris= "<<lEdTris.Count()<<" if X= "<<lX)
//LzLogM("", "Collapsed edge tris= "<<lEdTris.ListToString())
					lCanCollapse = false;
					break;
				}
			}

			// Skip if edge not collapsable
			if( !lCanCollapse )
			{
//LzLogM("", "Cannot collapse edge...")
				// Next candidate for W
				continue;
			}



			// Current candidate has passed all tests

			//W = ...
//}
//***********************************************************************

//LzLogM("", "Edge "<<v<<" - "<<w<<" is collapsable")


			//--------------------------------
			// Collapse edge
			//--------------------------------

			// Mark triangles as rejected
			BrowseList( iT, lRej )
				lRejectTri[ lRej.GetAt(iT) ] = true;

//LzLogM("", "Rejecting triangles "<<lRej.ListToString())

			// Mark vertices as merged
			lMergesWith[v] = w;
			lMergesWith[w] = v;

#if 1
			// Find best merged position to preserve original surface
			{
				// Sum square volumes of local tets depending on merged vertex v'
				//
				// Minimize a L^2 + b L + c, for L in [0, 1]
				//
				// Then v' = v + L vw, i.e. vv' = L vw
				//
				double a = 0;
				double b = 0;

				// V, W, and VW
				const Point3D & V = pMesh.mVertices[v];
				const Point3D & W = pMesh.mVertices[w];
				const Vector3D VW = W - V;

				// Contribution from V
				BrowseList( iT, lTrisSharing_V )
				{
                    const size_t lTri = lTrisSharing_V.GetAt(iT);

					// Skip if rejected
					if( lRejectTri[ lTri ] )
						continue;

					// Edges
					const TopTri & lTT = lTopo.TopTris()[ lTri ];
					const int lE0 = lTT.EdgeWithVertex( v );
					const int lE1 = lTT.EdgeWithVertex( v, lE0 );

					// Points
					const Point3D & A = pMesh.mVertices[ lTopo.TopEdges()[lE0].OtherVer(v) ];
					const Point3D & B = pMesh.mVertices[ lTopo.TopEdges()[lE1].OtherVer(v) ];

					// Vol
					const double lVol = ((A - V)^(B - V)) * VW;

					// Accum
					a += lVol*lVol;
				}

				// Contribution from W
				BrowseList( iT, lTrisSharing_W )
				{
                    const size_t lTri = lTrisSharing_W.GetAt(iT);

					// Skip if rejected
					if( lRejectTri[ lTri ] )
						continue;

					// Edges
					const TopTri & lTT = lTopo.TopTris()[ lTri ];
					const int lE0 = lTT.EdgeWithVertex( w );
					const int lE1 = lTT.EdgeWithVertex( w, lE0 );

					// Points
					const Point3D & A = pMesh.mVertices[ lTopo.TopEdges()[lE0].OtherVer(w) ];
					const Point3D & B = pMesh.mVertices[ lTopo.TopEdges()[lE1].OtherVer(w) ];

					// Vol
					const double lVol = ((A - W)^(B - W)) * VW;

					// Accum
					a += lVol*lVol;
					b += -2 * lVol * lVol;
				}

				// Minimize
				double L;
                if( LzMath::ToolBox::IsZero(2*a) )
					L = 1.0;
				else
					L = -b / (2*a);

				// Clamp
				if( L < 0 ) L = 0;
				if( L > 1 ) L = 1;

				// Set optimal position
				pMesh.mVertices[v] = V + L * VW;
				pMesh.mVertices[w] = pMesh.mVertices[v];
			}
#else
			{
				const Point3D & V = pMesh.mVertices[v];
				const Point3D & W = pMesh.mVertices[w];

				pMesh.mVertices[v] = V.MidPointTo( W );
				pMesh.mVertices[w] = pMesh.mVertices[v];
			}
#endif

			// Mark vertices that lost their topo information
			{
				// List of all contaminated vertices
                List<size_t> lLost;

				// Lambda
                auto addVers = [&]( List<size_t> & pTo, size_t i )
				{
                    const List<size_t> & lEds = lTopo.TopVers()[ i ].mE;
					BrowseList( iE, lEds )
						pTo.AddIntoIncrList( lTopo.TopEdges()[ lEds.GetAt(iE) ].OtherVer( i ), /*pUnique= */true );
				};

				// Add vertices that are edge-connected to either v or w
				// v will add w, and w will add v
				addVers( lLost, v );
				addVers( lLost, w );

				// Flag all
				BrowseList( iL, lLost )
					lLostTopo[ lLost.GetAt(iL) ] = true;
			}

			// Update mergers count
			lMergers++;
		}
#pragma endregion


		//-------------------------------------
#pragma region "Collapse vertices on edges"
		//-------------------------------------

		// Reset projetctor
        size_t lProjectors = 0;

		// Try to project all vertices
        //for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		//{
		//*********************** Verifier que lLostTopo marche aussi pour les projectors
		//}
		//
		//
		//*** TODO
		//
		//
		//*** TODO
		//
		//
#pragma endregion


		//--------------------------------
		// Finish if not merged anything
		//--------------------------------

		if( lMergers==0 && lProjectors==0 )
			break;

		//--------------------------------
		// Create new vers and mapping
		//--------------------------------

//**************** adapter le mapping pour prendre en compte les projectors
//**************** adapter le mapping pour prendre en compte les projectors
//**************** adapter le mapping pour prendre en compte les projectors
//**************** adapter le mapping pour prendre en compte les projectors
//**************** adapter le mapping pour prendre en compte les projectors
//**************** adapter le mapping pour prendre en compte les projectors

		vector<Point3D> lNewVertices( pMesh.mVertices.size() - lMergers );
        vector<size_t> lOldToNewVerIdx( pMesh.mVertices.size() );
		{
            size_t iNewIdx = 0;
            for( size_t old_idx=0 ; old_idx<lMergesWith.size() ; old_idx++ )
			{
				// Read merger
                const size_t w = lMergesWith[ old_idx ];

				// Vertex is merging?
				if( w < pMesh.mVertices.size() )
				{
					// Yes, merging

					// Do once
					if( old_idx < w )
					{
						lNewVertices[ iNewIdx ] = pMesh.mVertices[ old_idx ];
						lOldToNewVerIdx[ old_idx ] = iNewIdx;
						lOldToNewVerIdx[       w ] = iNewIdx;
						iNewIdx++;
					}
				}
				else
				{
					// Nope, not merging

					// Copy previous position
					lNewVertices[ iNewIdx ] = pMesh.mVertices[ old_idx ];
					lOldToNewVerIdx[ old_idx ] = iNewIdx;
					iNewIdx++;
				}
			}
		}

		//--------------------------------
		// Normals are lost
		//--------------------------------

		pMesh.mNormals = { Vector3D(1, 0, 0) };

		//--------------------------------
		// Create new triangles
		//--------------------------------

//************* DEBUG
//************* DEBUG
{
LzLogM("", "**** DEBUG: before, lLostTris= "<<lLostTris)
	lLostTris = 0;
    for( size_t t=0 ; t<pMesh.mTriangles.size() ; t++ )
		if( lRejectTri[ t ] )
			lLostTris++;
LzLogM("", "**** DEBUG: after, lLostTris= "<<lLostTris)
}
//************* DEBUG
//************* DEBUG


		vector<Triangle> lNewTriangles( pMesh.mTriangles.size() - lLostTris );
		{
			// Create new array
            size_t iIdx = 0;
            for( size_t old=0 ; old<pMesh.mTriangles.size() ; old++ )
			{
				// Keeping this triangle?
				if( !lRejectTri[ old ] )
				{
					// Create new triangle
					const Triangle & lOld = pMesh.mTriangles[ old ];
					lNewTriangles[ iIdx ] = Triangle( lOldToNewVerIdx[ lOld.mIdxV[0] ],
													  lOldToNewVerIdx[ lOld.mIdxV[1] ],
													  lOldToNewVerIdx[ lOld.mIdxV[2] ], 0, 0, 0 ); // Pointing to the unique normal

// Check
if( lNewTriangles[ iIdx ].HasRepeatedVers() )
    LzLogException("", "Found a triangle with repeated vertices! Tri= "<<lNewTriangles[ iIdx ].ToString())

					// Next
					iIdx++;
				}
			}
		}

		//--------------------------------
		// Commit
		//--------------------------------

		pMesh.mVertices = lNewVertices;
		pMesh.mTriangles = lNewTriangles;


////***** DEBUG
//{
//	pMesh.SmoothNormalize();
//	pMesh.Save(LzServices::StartUpPath()+"/_tmp_TopoClean_it_"+std::to_string(iStep)+".obj");
//}
////***** DEBUG


		// Log
        LzLogM("", "Mesh::TopoClean: iteration "<<iStep<<": vertices= "<<pMesh.mVertices.size()<<", triangles= "<<pMesh.mTriangles.size()<<".")
		iStep++;
	}

    // Need to recompute DL
	pMesh.ResetDisplayList();
}
#endif
#pragma endregion


#pragma region "Mesh relaxation"
//================================================================================
void SlideVerticesToMeanPos( Mesh & pMesh, const FastDistComp & pOnto, size_t pSteps,
                             /*size_t pVerNeighs= 1,*/ const MeshTopology * ppTopo/*=nullptr*/, const Vector<bool> * pIsVerFixed/*=nullptr*/ )
{
    // Check
    if( pIsVerFixed && pMesh.mVertices.size() != pIsVerFixed->Size() )
        LzLogException("",  "Inconsistent vertex count! Mesh has " << pMesh.mVertices.size() << " vertices, while vector has " << pIsVerFixed->Size() << "." )

//// Check
//if( pVerNeighs == 0 )
//	LzLogException("", "Cannot slide vertices to mean pos using a 0-vertices neighborhood!")

//// Build topo info
//Vector< List<size_t> > lNeighVers( pMesh.mVertices.size() );
//{
		// Pointer to topology: provided by client or my own personal one.
		const MeshTopology * lpTopo;

		MeshTopology lMyOwnTopo;
		if( ppTopo )
			lpTopo = ppTopo;
		else
		{
			lMyOwnTopo.Set( pMesh );
			lpTopo = &lMyOwnTopo;
		}


//**********************************************************************
// Links to topology
const vector<TopVer> & lTopVers = lpTopo->TopVers();
const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();
//**********************************************************************


//		// Temporary done vector
//		Vector<bool> lDone;
//		lDone.Assign( pMesh.mVertices.size(), false );
//
//		// Find neighbors for all vertices
//        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
//		{
//			// Mark this vertex as done
//			lDone[v] = true;
//
//			//--------------------------------
//			// Init first level neighbors
//			//--------------------------------
//			{
//                // Edges for this vertex
//                const List<size_t> & lEdges = lTopVers[v].mE;
//                BrowseList( iE, lEdges )
//                {
//					const size_t lOther = lTopEdges[lEdges.GetAt( iE )].OtherVer( v );
//					if( !lDone[lOther] )
//					{
//						lNeighVers[v].AddTail( lOther );
//						lDone[lOther] = true;
//					}
//				}
//			}
//
//			//--------------------------------
//			// Process higher level neighbors
//			//--------------------------------
//			for( size_t n=1 ; n<pVerNeighs ; n++ )
//			{
//				// Find new neighbors
//				List<size_t> lNewNeighs;
//				BrowseList( iN, lNeighVers[v] )
//				{
//					// Previously found neighbor
//					const size_t lNeigh = lNeighVers[v].GetAt( iN );
//
//					// Add neighbors of this neighbor
//					const List<size_t> & lNeighEdges = lTopVers[lNeigh].mE;
//					BrowseList( iNE, lNeighEdges )
//					{
//						const size_t lOther = lTopEdges[lNeighEdges.GetAt( iNE )].OtherVer( lNeigh );
//						if( !lDone[lOther] )
//						{
//							lNewNeighs.AddTail( lOther );
//							lDone[lOther] = true;
//						}
//					}
//				}
//
//				// Commit new neighbors
//				lNeighVers[v].Append( lNewNeighs );
//			}
//
//			// Reset flags
//			lDone[v] = false;
//			BrowseList( iN, lNeighVers[v] )
//				lDone[ lNeighVers[v].GetAt(iN) ] = false;
//		}
//
////**************************** DEBUG
//// Check
//{
//    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
//	{
//		if( lDone[v] )
//			LzLogException("", "Oopsie!")
//
//		List<size_t> lCheckList;
//		BrowseList( iN, lNeighVers[v] )
//			lCheckList.AddIntoIncrList( lNeighVers[v].GetAt(iN), true );
//
//		if( lCheckList.Count() != lNeighVers[v].Count() )
//			LzLogException("", "Daisie!")
//
//		if( lCheckList.FindElt(v) )
//			LzLogException("", "Well, fuck!")
//	}
//}
////**************************** DEBUG
//	}

    // Create new table of vertices
    vector<Point3D> lNewVertices( pMesh.mVertices.size() );
//**********************
//********************** initialiser new vertices une fois et skipper les noeuds fixes dans la boucle de relax
//**********************
//**********************

    // Loop all steps
    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
        // Compute mean positions
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
        {
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
            {
                // Don't move!
                lNewVertices[v] = pMesh.mVertices[v];
//**********************
//********************** initialiser new vertices une fois et skipper les noeuds fixes dans la boucle de relax
//**********************
//**********************
            }
            else
            {
#if 0
				// New vertex position
Point3D lNewVer(0, 0, 0); // Pure Laplacian
//Point3D lNewVer = pMesh.mVertices[v];

                // New mean point and neighbors count
                size_t lNeighborsCount = lNeighVers[v].Count();
				if( lNeighborsCount )
				{
					// Add
					BrowseList( iN, lNeighVers[v] )
					{
						// Neighboring vertex
						const Point3D & lOtherVer = pMesh.mVertices[ lNeighVers[v].GetAt(iN) ];

						// Acc
						lNewVer.X() += lOtherVer.X();
						lNewVer.Y() += lOtherVer.Y();
						lNewVer.Z() += lOtherVer.Z();
					}

					// Normalize
lNewVer /= lNeighborsCount; // Pure Laplacian
//lNewVer /= lNeighborsCount + 1;
				}
				else
				{
					lNewVer = pMesh.mVertices[v];
				}
#else
                // Edges for this vertex
                const List<size_t> & lEdges = lTopVers[v].mE;

                Point3D lNewVer = pMesh.mVertices[v];
                size_t lNeighborsCount = 1 + lEdges.Count();

                BrowseList( iE, lEdges )
                {
                    // Get edge connected vertex
                    const Point3D & lOtherVer = pMesh.mVertices[ lTopEdges[lEdges.GetAt( iE )].OtherVer( v ) ];

                    // Acc
                    lNewVer.X() += lOtherVer.X();
                    lNewVer.Y() += lOtherVer.Y();
                    lNewVer.Z() += lOtherVer.Z();
                }

                // Normalize
                lNewVer /= lNeighborsCount;
#endif
                // Stash
                lNewVertices[v] = lNewVer;
            }
        }

		// Resolve constraints
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
				continue;

			// Apply constraint to free vertices only
			Point3D lProj;
			Vector3D lTmp;
			if( pOnto.SignedDistance(Tree3D::NearestMode::Accurate, lNewVertices[v], lProj, lTmp) < 0 )
				lNewVertices[v] = lProj;
		}

        // Commit
        pMesh.mVertices = lNewVertices;
    }

    // Need to recompute DL
    pMesh.ResetDisplayList();
}

//================================================================================
void SmoothSmoothing( Mesh & pMesh, size_t pSteps, /*size_t pTriNeighs=1, */const MeshTopology * ppTopo/*=nullptr*/, const Vector<bool> * pIsVerFixed/*=nullptr*/ )
{
    // Check
    if( pIsVerFixed && pMesh.mVertices.size() != pIsVerFixed->Size() )
        LzLogException("",  "Inconsistent vertex count! Mesh has " << pMesh.mVertices.size() << " vertices, while vector has " << pIsVerFixed->Size() << "." )

//// Check
//if( pVerNeighs == 0 )
//	LzLogException("", "Cannot slide vertices to mean pos using a 0-vertices neighborhood!")

	// Build topo info
    const MeshTopology * lpTopo;
    Vector< List<size_t> > lNeighVers( pMesh.mVertices.size() );
	{
		// Pointer to topology: provided by client or my own personal one.

		MeshTopology lMyOwnTopo;
		if( ppTopo )
			lpTopo = ppTopo;
		else
		{
			lMyOwnTopo.Set( pMesh );
			lpTopo = &lMyOwnTopo;
		}


//**********************************************************************
// Links to topology
const vector<TopVer> & lTopVers = lpTopo->TopVers();
const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();
//**********************************************************************


		// Temporary done vector
		Vector<bool> lDone;
		lDone.Assign( pMesh.mVertices.size(), false );

		// Find neighbors for all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Mark this vertex as done
			lDone[v] = true;

			//--------------------------------
			// Init first level neighbors
			//--------------------------------
			{
                // Edges for this vertex
                const List<size_t> & lEdges = lTopVers[v].mE;
                BrowseList( iE, lEdges )
                {
                    const size_t lOther = lTopEdges[lEdges.GetAt( iE )].OtherVer( v );
					if( !lDone[lOther] )
					{
						lNeighVers[v].AddTail( lOther );
						lDone[lOther] = true;
					}
				}
			}

			//--------------------------------
			// Process higher level neighbors
			//--------------------------------
            //for( size_t n=1 ; n<pVerNeighs ; n++ )
			//{
			//	// Find new neighbors
            //	List<size_t> lNewNeighs;
			//	BrowseList( iN, lNeighVers[v] )
			//	{
			//		// Previously found neighbor
            //		const size_t lNeigh = lNeighVers[v].GetAt( iN );
			//
			//		// Add neighbors of this neighbor
            //		const List<size_t> & lNeighEdges = lTopVers[lNeigh].mE;
			//		BrowseList( iNE, lNeighEdges )
			//		{
            //			const size_t lOther = lTopEdges[lNeighEdges.GetAt( iNE )].OtherVer( lNeigh );
			//			if( !lDone[lOther] )
			//			{
			//				lNewNeighs.AddTail( lOther );
			//				lDone[lOther] = true;
			//			}
			//		}
			//	}
			//
			//	// Commit new neighbors
			//	lNeighVers[v].Append( lNewNeighs );
			//}

			// Reset flags
			lDone[v] = false;
			BrowseList( iN, lNeighVers[v] )
				lDone[ lNeighVers[v].GetAt(iN) ] = false;
		}

//**************************** DEBUG
// Check
if( 0 )
{
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		if( lDone[v] )
            LzLogException("", "Oopsie!")

        List<size_t> lCheckList;
		BrowseList( iN, lNeighVers[v] )
			lCheckList.AddIntoIncrList( lNeighVers[v].GetAt(iN), true );

		if( lCheckList.Count() != lNeighVers[v].Count() )
            LzLogException("", "Daisie!")

		if( lCheckList.FindElt(v) )
            LzLogException("", "Well, fuck!")
	}
}
//**************************** DEBUG
	}


//******************** create FDC
FastDistComp lFDC;
lFDC.Set( pMesh, *lpTopo, LzTriModel::Patching::Split );
//******************** create FDC


    // Create new table of vertices, and table of retro-projections
    vector<Point3D> lNewVertices( pMesh.mVertices.size() );
//**********************
//********************** initialiser new vertices une fois et skipper les noeuds fixes dans la boucle de relax
//**********************
//**********************

	vector<Vector3D> lProjVertices( pMesh.mVertices.size() );

    // Loop all steps
    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
		//------------------------------------
        // Compute mean positions
		//------------------------------------

        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
        {
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
            {
                // Don't move!
                lNewVertices[v] = pMesh.mVertices[v];
            }
            else
            {
				// New vertex position
				Point3D lNewVer(0, 0, 0); // Pure Laplacian

                // New mean point and neighbors count
                size_t lNeighborsCount = lNeighVers[v].Count();
				if( lNeighborsCount )
				{
					// Add
					BrowseList( iN, lNeighVers[v] )
					{
						// Neighboring vertex
						const Point3D & lOtherVer = pMesh.mVertices[ lNeighVers[v].GetAt(iN) ];

						// Acc
						lNewVer.X() += lOtherVer.X();
						lNewVer.Y() += lOtherVer.Y();
						lNewVer.Z() += lOtherVer.Z();
					}

					// Normalize
					lNewVer /= lNeighborsCount; // Pure Laplacian
				}
				else
				{
					lNewVer = pMesh.mVertices[v];
				}

				// Stash
                lNewVertices[v] = lNewVer;
            }
        }

		//------------------------------------
        // Compute projections
		//------------------------------------

        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
            // Skip fixed vertices
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
				continue;

			// Compute projection for moving vertices
			Point3D lProj;
			Vector3D tmp;
            lFDC.SignedDistance( LzGeom::Tree3D::NearestMode::Accurate, lNewVertices[v], lProj, tmp);
			//
			lProjVertices[v] = lProj - lNewVertices[v];
		}

		//------------------------------------
        // Smooth projections
		//------------------------------------

		{
			// Buffer for smoothing
			vector<Vector3D> lNewProjVertices = lProjVertices;

//****************************************************
            const size_t SMOOTH = 200;
//            const size_t SMOOTH = 20;
//****************************************
//****************************************
//**
//**
//**
//** CA MARCHE !!!!!!
//**
//**
//**
//****************************************
//****************************************
//const size_t SMOOTH = 1;
//****************************************************
            for( size_t sm=0 ; sm<SMOOTH ; sm++ )
			{
                for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
				{
					// Skip fixed vertices
					if( pIsVerFixed && ( *pIsVerFixed )[v] )
						continue;

					// New mean point and neighbors count
                    size_t lNeighborsCount = lNeighVers[v].Count();
					if( lNeighborsCount )
					{
						// Add
						Vector3D lNewProj(0, 0, 0); // Pure Laplacian
						BrowseList( iN, lNeighVers[v] )
							lNewProj += lProjVertices[ lNeighVers[v].GetAt(iN) ];

						// Normalize
						lNewProj /= lNeighborsCount; // Pure Laplacian

						// Stash
						lNewProjVertices[v] = lNewProj;
					}
				}

				// Commit
				lProjVertices = lNewProjVertices;
			}
		}

		//------------------------------------
        // Apply projections
		//------------------------------------

        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			lNewVertices[v] += lProjVertices[v];
		}

        // Commit
        pMesh.mVertices = lNewVertices;
    }

    // Need to recompute DL
    pMesh.ResetDisplayList();
}


#pragma region "....failures"
//================================================================================
void TangentMeshSmoothing( Mesh & pMesh, size_t pSteps, size_t pTriNeighs/*=1*/, const MeshTopology * ppTopo/*=nullptr*/, const Vector<bool> * pIsVerFixed/*=nullptr*/ )
{
    // Check
    if( pIsVerFixed && pMesh.mVertices.size() != pIsVerFixed->Size() )
        LzLogException("",  "Inconsistent vertex count! Mesh has " << pMesh.mVertices.size() << " vertices, while vector has " << pIsVerFixed->Size() << "." )

	// Check
	if( pTriNeighs == 0 )
        LzLogException("", "Cannot tangent mesh smooth using a 0-triangles neighborhood!")

	// Build topo info
    Vector< List<size_t> > lNeigTris( pMesh.mVertices.size() );
	{
		// Pointer to topology: provided by client or my own personal one.
		const MeshTopology * lpTopo;

		MeshTopology lMyOwnTopo;
		if( ppTopo )
			lpTopo = ppTopo;
		else
		{
			lMyOwnTopo.Set( pMesh );
			lpTopo = &lMyOwnTopo;
		}


//**********************************************************************
// Links to topology
const vector<TopVer> & lTopVers = lpTopo->TopVers();
const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();
const vector<TopTri> & lTopTris = lpTopo->TopTris();
//**********************************************************************


		// Temporary done vector
		Vector<bool> lDone;
		lDone.Assign( pMesh.mTriangles.size(), false );

		// Find neighbor tris for all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			//--------------------------------
			// Init first level neighbors
			//--------------------------------
			{
                // Edges for this vertex
                const List<size_t> & lTris = lTopVers[v].mT;
                BrowseList( iT, lTris )
                {
                    const size_t lTri = lTris.GetAt( iT );
					if( !lDone[lTri] )
					{
						lNeigTris[v].AddTail( lTri );
						lDone[lTri] = true;
					}
				}
			}

			//--------------------------------
			// Process higher level neighbors
			//--------------------------------
            for( size_t n=1 ; n<pTriNeighs ; n++ )
			{
				// Find new neighbors
                List<size_t> lNewTris;
				BrowseList( iT, lNeigTris[v] )
				{
					// Previously found neighbor
                    const size_t lNeigh = lNeigTris[v].GetAt( iT );

					// Add neighbors of this neighbor
					for( int e=0 ; e<3 ; e++ )
					{
                        const size_t lOther = lTopEdges[ lTopTris[lNeigh].mE[e] ].OtherTri( lNeigh );
						if( !lDone[lOther] )
						{
							lNewTris.AddTail( lOther );
							lDone[lOther] = true;
						}
					}
				}

				// Commit new neighbors
				lNeigTris[v].Append( lNewTris );
			}

			// Reset flags
			BrowseList( iT, lNeigTris[v] )
				lDone[ lNeigTris[v].GetAt(iT) ] = false;
		}

//**************************** DEBUG
// Check
{
    LzLogN("", "** CHECK");
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		if( lDone[v] )
            LzLogException("", "Oopsie!")

        List<size_t> lCheckList;
		BrowseList( iN, lNeigTris[v] )
			lCheckList.AddIntoIncrList( lNeigTris[v].GetAt(iN), true );

		if( lCheckList.Count() != lNeigTris[v].Count() )
            LzLogException("", "Daisie!")

        LzLogM("", "Triangles for ver "<<v<<": "<<lNeigTris[v].Count())
	}
}
//**************************** DEBUG
	}





	// Precompute data
	Vector<Point3D> lAllCenters;
	Vector<double> lAllAreas;
	pMesh.GetFacesCentroids( lAllCenters, lAllAreas );

    // Create new table of vertices
    vector<Point3D> lNewVertices( pMesh.mVertices.size() );

    // Loop all steps
    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
        // Compute mean positions
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
        {
            //// Fixed vertex?
            //if( pIsVerFixed && ( *pIsVerFixed )[v] )
            //{
            //    // Don't move!
            //    lNewVertices[v] = pMesh.mVertices[v];
            //}
            //else
            {
                //const LzTriModel::TopVer & lTV = lTopVers[v];

				// Collect data
				Vector<Point3D> lCenters;
				Vector<double> lAreas;
				BrowseList( iT, lNeigTris[v] )
				{
                    const size_t lT = lNeigTris[v].GetAt( iT );
					lCenters.PushBack( lAllCenters[lT] );
					lAreas.PushBack( lAllAreas[lT] );
				}

				// Compute PCA
				Point3D lMean;
				Vector3D lEigV[3];
                LzMath::ToolBox::WeightedPCA( lCenters, lAreas, lMean, lEigV );

				// Project
const Plane3D lPi( lMean, lEigV[2] );
//const Plane3D lPi( lMean, lEigV[0], lEigV[1] );
				lNewVertices[v] = lPi.Projection( pMesh.mVertices[v] );
			}
		}

        // Commit
        pMesh.mVertices = lNewVertices;
    }

    // Need to recompute DL
    pMesh.ResetDisplayList();
}

//================================================================================
void TanVerSmooth( Mesh & pMesh/*, const FastDistComp & pOnto*/, size_t pSteps, size_t pVerNeighs/*=1*/, const MeshTopology * ppTopo/*=nullptr*/, const Vector<bool> * pIsVerFixed/*=nullptr*/ )
{
    // Check
    if( pIsVerFixed && pMesh.mVertices.size() != pIsVerFixed->Size() )
        LzLogException("",  "Inconsistent vertex count! Mesh has " << pMesh.mVertices.size() << " vertices, while vector has " << pIsVerFixed->Size() << "." )

	// Check
	if( pVerNeighs == 0 )
        LzLogException("", "Cannot slide vertices to mean pos using a 0-vertices neighborhood!")

	// Pointer to topology: provided by client or my own personal one.
	const MeshTopology * lpTopo;
	MeshTopology lMyOwnTopo;
	if( ppTopo )
		lpTopo = ppTopo;
	else
	{
		lMyOwnTopo.Set( pMesh );
		lpTopo = &lMyOwnTopo;
	}

	// Build topo info
    Vector< List<size_t> > lNeighVers( pMesh.mVertices.size() );
	{
//**********************************************************************
// Links to topology
const vector<TopVer> & lTopVers = lpTopo->TopVers();
const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();
//**********************************************************************


		// Temporary done vector
		Vector<bool> lDone;
		lDone.Assign( pMesh.mVertices.size(), false );

		// Find neighbors for all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Mark this vertex as done
			lDone[v] = true;

			//--------------------------------
			// Init first level neighbors
			//--------------------------------
			{
                // Edges for this vertex
                const List<size_t> & lEdges = lTopVers[v].mE;
                BrowseList( iE, lEdges )
                {
                    const size_t lOther = lTopEdges[lEdges.GetAt( iE )].OtherVer( v );
					if( !lDone[lOther] )
					{
						lNeighVers[v].AddTail( lOther );
						lDone[lOther] = true;
					}
				}
			}

			//--------------------------------
			// Process higher level neighbors
			//--------------------------------
            for( size_t n=1 ; n<pVerNeighs ; n++ )
			{
				// Find new neighbors
                List<size_t> lNewNeighs;
				BrowseList( iN, lNeighVers[v] )
				{
					// Previously found neighbor
                    const size_t lNeigh = lNeighVers[v].GetAt( iN );

					// Add neighbors of this neighbor
                    const List<size_t> & lNeighEdges = lTopVers[lNeigh].mE;
					BrowseList( iNE, lNeighEdges )
					{
                        const size_t lOther = lTopEdges[lNeighEdges.GetAt( iNE )].OtherVer( lNeigh );
						if( !lDone[lOther] )
						{
							lNewNeighs.AddTail( lOther );
							lDone[lOther] = true;
						}
					}
				}

				// Commit new neighbors
				lNeighVers[v].Append( lNewNeighs );
			}

			// Reset flags
			lDone[v] = false;
			BrowseList( iN, lNeighVers[v] )
				lDone[ lNeighVers[v].GetAt(iN) ] = false;
		}

//**************************** DEBUG
// Check
{
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		if( lDone[v] )
            LzLogException("", "Oopsie!")

        List<size_t> lCheckList;
		BrowseList( iN, lNeighVers[v] )
			lCheckList.AddIntoIncrList( lNeighVers[v].GetAt(iN), true );

		if( lCheckList.Count() != lNeighVers[v].Count() )
            LzLogException("", "Daisie!")

		if( lCheckList.FindElt(v) )
            LzLogException("", "Well, fuck!")
	}
}
//**************************** DEBUG
	}

    // Create new table of vertices
    vector<Point3D> lNewVertices( pMesh.mVertices.size() );

////******************* SAME PLAYER SHOOTS AGAIN
//// Smooth normalize
//pMesh.SmoothNormalize( lpTopo );
////******************* SAME PLAYER SHOOTS AGAIN

    // Loop all steps
    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
//******************* SAME PLAYER SHOOTS AGAIN
// Smooth normalize
pMesh.SmoothNormalize( lpTopo );
//******************* SAME PLAYER SHOOTS AGAIN

        // Compute mean positions
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
        {
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
            {
                // Don't move!
                lNewVertices[v] = pMesh.mVertices[v];
            }
            else
            {
				// New vertex position
				Point3D lNewVer;

                // New mean point and neighbors count
                size_t lNeighborsCount = lNeighVers[v].Count();
				if( lNeighborsCount )
				{
#if 1
					const Point3D & lVer = pMesh.mVertices[v];
					const Vector3D & lNor = pMesh.mNormals[v];

					double L = 0;
					BrowseList( iN, lNeighVers[v] )
						L += (pMesh.mVertices[ lNeighVers[v].GetAt(iN) ] - lVer)*lNor;
					L /= lNeighborsCount;

					if( L > 0 )
						lNewVer = lVer + L*lNor;
					else
						lNewVer = lVer;
					//lNewVer = lVer + L*lNor;
#else
					// Gather all neighbors
					Vector<Point3D> lNeighs( lNeighborsCount );
                    size_t n = 0;
					BrowseList( iN, lNeighVers[v] )
						lNeighs[n++] = pMesh.mVertices[ lNeighVers[v].GetAt(iN) ];

					// Compute LSQ plane and project
					Plane3D lLSQ( lNeighs );

					//lNewVer = lLSQ.Projection( pMesh.mVertices[v] );
//******************* SAME PLAYER SHOOTS AGAIN
					const Point3D & lVer = pMesh.mVertices[v];
					if( lLSQ.Normal() * pMesh.mNormals[v] * lLSQ.SignedDistanceTo(lVer) < 0 )
						lNewVer = lLSQ.Projection( lVer );
					else
						lNewVer = lVer;
#endif
//******************* SAME PLAYER SHOOTS AGAIN
				}
				else
				{
					lNewVer = pMesh.mVertices[v];
				}

				// Stash
                lNewVertices[v] = lNewVer;
            }
        }

//		// Resolve constraints
//if( 1 )
//		for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
//		{
//            // Fixed vertex?
//            if( pIsVerFixed && ( *pIsVerFixed )[v] )
//				continue;
//
//			// Apply constraint to free vertices only
//			Point3D lProj;
//			if( pOnto.SignedDistance(Tree3D::NearestMode::Accurate, lNewVertices[v], lProj, Vector3D()) < 0 )
//				lNewVertices[v] = lProj;
//		}

        // Commit
        pMesh.mVertices = lNewVertices;
    }

    // Need to recompute DL
    pMesh.ResetDisplayList();
}

//================================================================================
void TanVerSmooth_2( Mesh & pMesh, const FastDistComp & pOnto, size_t pSteps, size_t pVerNeighs/*=1*/, const MeshTopology * ppTopo/*=nullptr*/, const Vector<bool> * pIsVerFixed/*=nullptr*/ )
{
    // Check
    if( pIsVerFixed && pMesh.mVertices.size() != pIsVerFixed->Size() )
        LzLogException("",  "Inconsistent vertex count! Mesh has " << pMesh.mVertices.size() << " vertices, while vector has " << pIsVerFixed->Size() << "." )

	// Check
	if( pVerNeighs == 0 )
        LzLogException("", "Cannot slide vertices to mean pos using a 0-vertices neighborhood!")

	// Build topo info
    Vector< List<size_t> > lNeighVers( pMesh.mVertices.size() );
	{
		// Pointer to topology: provided by client or my own personal one.
		const MeshTopology * lpTopo;

		MeshTopology lMyOwnTopo;
		if( ppTopo )
			lpTopo = ppTopo;
		else
		{
			lMyOwnTopo.Set( pMesh );
			lpTopo = &lMyOwnTopo;
		}


//**********************************************************************
// Links to topology
const vector<TopVer> & lTopVers = lpTopo->TopVers();
const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();
//**********************************************************************


		// Temporary done vector
		Vector<bool> lDone;
		lDone.Assign( pMesh.mVertices.size(), false );

		// Find neighbors for all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Mark this vertex as done
			lDone[v] = true;

			//--------------------------------
			// Init first level neighbors
			//--------------------------------
			{
                // Edges for this vertex
                const List<size_t> & lEdges = lTopVers[v].mE;
                BrowseList( iE, lEdges )
                {
                    const size_t lOther = lTopEdges[lEdges.GetAt( iE )].OtherVer( v );
					if( !lDone[lOther] )
					{
						lNeighVers[v].AddTail( lOther );
						lDone[lOther] = true;
					}
				}
			}

			//--------------------------------
			// Process higher level neighbors
			//--------------------------------
            for( size_t n=1 ; n<pVerNeighs ; n++ )
			{
				// Find new neighbors
                List<size_t> lNewNeighs;
				BrowseList( iN, lNeighVers[v] )
				{
					// Previously found neighbor
                    const size_t lNeigh = lNeighVers[v].GetAt( iN );

					// Add neighbors of this neighbor
                    const List<size_t> & lNeighEdges = lTopVers[lNeigh].mE;
					BrowseList( iNE, lNeighEdges )
					{
                        const size_t lOther = lTopEdges[lNeighEdges.GetAt( iNE )].OtherVer( lNeigh );
						if( !lDone[lOther] )
						{
							lNewNeighs.AddTail( lOther );
							lDone[lOther] = true;
						}
					}
				}

				// Commit new neighbors
				lNeighVers[v].Append( lNewNeighs );
			}

			// Reset flags
			lDone[v] = false;
			BrowseList( iN, lNeighVers[v] )
				lDone[ lNeighVers[v].GetAt(iN) ] = false;
		}

//**************************** DEBUG
// Check
{
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		if( lDone[v] )
            LzLogException("", "Oopsie!")

        List<size_t> lCheckList;
		BrowseList( iN, lNeighVers[v] )
			lCheckList.AddIntoIncrList( lNeighVers[v].GetAt(iN), true );

		if( lCheckList.Count() != lNeighVers[v].Count() )
            LzLogException("", "Daisie!")

		if( lCheckList.FindElt(v) )
            LzLogException("", "Well, fuck!")
	}
}
//**************************** DEBUG
	}



//// Create new table of vertices
//vector<Point3D> lNewVertices( pMesh.mVertices.size() );




    // Loop all steps
    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
//**** BUFFER
vector< List<Point3D> > lProjLists( pMesh.mVertices.size() );
//**** BUFFER



        // Compute mean positions
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
        {
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
            {
//// Don't move!
//lNewVertices[v] = pMesh.mVertices[v];
            }
            else
            {
				// New vertex position
				Point3D lNewVer;

                // New mean point and neighbors count
                size_t lNeighborsCount = lNeighVers[v].Count();
				if( lNeighborsCount )
				{
					// Gather all neighbors
					Vector<Point3D> lNeighs( lNeighborsCount );
                    size_t n = 0;
					BrowseList( iN, lNeighVers[v] )
						lNeighs[n++] = pMesh.mVertices[ lNeighVers[v].GetAt(iN) ];

					// Compute LSQ plane and project
					Plane3D lLSQ( lNeighs );

					// Project self
					lProjLists[ v ].AddTail( lLSQ.Projection( pMesh.mVertices[v] ) );

					// Project neighbors
					BrowseList( iN, lNeighVers[v] )
					{
						unsigned lIdx = lNeighVers[v].GetAt(iN);
						lProjLists[ lIdx ].AddTail( lLSQ.Projection( pMesh.mVertices[lIdx] ) );
					}
				}
				else
				{
//lNewVer = pMesh.mVertices[v];
				}

//// Stash
//lNewVertices[v] = lNewVer;
            }
        }


		// Create new table of vertices
		vector<Point3D> lNewVertices( pMesh.mVertices.size() );
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Fixed?
			if( lProjLists[v].Count() == 0 )
			{
				lNewVertices[v] = pMesh.mVertices[v];
			}
			else
			{
				Point3D lMean(0, 0, 0);
				BrowseList( iP, lProjLists[v] )
				{
					const Point3D & lProj = lProjLists[v].GetAt(iP);
					for( int i=0 ; i<3 ; i++ )
						lMean.mV[i] += lProj.mV[i];
				}
				lMean /= lProjLists[v].Count();

				// Stash
				lNewVertices[v] = lMean;
			}
		}


		// Resolve constraints
//if( 1 ) //********** il faut resoudre les contraintes sinon on squizze la forme
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
            // Fixed vertex?
            if( pIsVerFixed && ( *pIsVerFixed )[v] )
				continue;

			// Apply constraint to free vertices only
			Point3D lProj;
			Vector3D tmp;
			if( pOnto.SignedDistance(Tree3D::NearestMode::Accurate, lNewVertices[v], lProj, tmp) < 0 )
				lNewVertices[v] = lProj;
		}

        // Commit
        pMesh.mVertices = lNewVertices;
    }

    // Need to recompute DL
    pMesh.ResetDisplayList();
}

//================================================================================
void SmoothPerVerNormals( Mesh & pMesh, size_t pSteps, size_t pVerNeighs/*=1*/, const MeshTopology * ppTopo/*=nullptr*/ )
{
	// Check
	if( pMesh.mVertices.size() != pMesh.mNormals.size() )
        LzLogException("", "Cannot smooth per vertex normals! Found "<<pMesh.mVertices.size()<<" vertices vs "<<pMesh.mNormals.size()<<" normals.")

	// Check
	if( pVerNeighs == 0 )
        LzLogException("", "Cannot smooth per vertex normals using a 0-vertices neighborhood!")

	// Build topo info
    Vector< List<size_t> > lNeighVers( pMesh.mVertices.size() );
	{
		// Pointer to topology: provided by client or my own personal one.
		const MeshTopology * lpTopo;

		MeshTopology lMyOwnTopo;
		if( ppTopo )
			lpTopo = ppTopo;
		else
		{
			lMyOwnTopo.Set( pMesh );
			lpTopo = &lMyOwnTopo;
		}


//**********************************************************************
// Links to topology
const vector<TopVer> & lTopVers = lpTopo->TopVers();
const vector<TopEdge> & lTopEdges = lpTopo->TopEdges();
//**********************************************************************


		// Temporary done vector
		Vector<bool> lDone;
		lDone.Assign( pMesh.mVertices.size(), false );

		// Find neighbors for all vertices
        for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
		{
			// Mark this vertex as done
			lDone[v] = true;

			//--------------------------------
			// Init first level neighbors
			//--------------------------------
			{
                // Edges for this vertex
                const List<size_t> & lEdges = lTopVers[v].mE;
                BrowseList( iE, lEdges )
                {
                    const size_t lOther = lTopEdges[lEdges.GetAt( iE )].OtherVer( v );
					if( !lDone[lOther] )
					{
						lNeighVers[v].AddTail( lOther );
						lDone[lOther] = true;
					}
				}
			}

			//--------------------------------
			// Process higher level neighbors
			//--------------------------------
            for( size_t n=1 ; n<pVerNeighs ; n++ )
			{
				// Find new neighbors
                List<size_t> lNewNeighs;
				BrowseList( iN, lNeighVers[v] )
				{
					// Previously found neighbor
                    const size_t lNeigh = lNeighVers[v].GetAt( iN );

					// Add neighbors of this neighbor
                    const List<size_t> & lNeighEdges = lTopVers[lNeigh].mE;
					BrowseList( iNE, lNeighEdges )
					{
                        const size_t lOther = lTopEdges[lNeighEdges.GetAt( iNE )].OtherVer( lNeigh );
						if( !lDone[lOther] )
						{
							lNewNeighs.AddTail( lOther );
							lDone[lOther] = true;
						}
					}
				}

				// Commit new neighbors
				lNeighVers[v].Append( lNewNeighs );
			}

			// Reset flags
			lDone[v] = false;
			BrowseList( iN, lNeighVers[v] )
				lDone[ lNeighVers[v].GetAt(iN) ] = false;
		}

//**************************** DEBUG
// Check
{
    for( size_t v=0 ; v<pMesh.mVertices.size() ; v++ )
	{
		if( lDone[v] )
            LzLogException("", "Oopsie!")

        List<size_t> lCheckList;
		BrowseList( iN, lNeighVers[v] )
			lCheckList.AddIntoIncrList( lNeighVers[v].GetAt(iN), true );

		if( lCheckList.Count() != lNeighVers[v].Count() )
            LzLogException("", "Daisie!")

		if( lCheckList.FindElt(v) )
            LzLogException("", "Well, fuck!")
	}
}
//**************************** DEBUG
	}

    // Loop all steps
    for( size_t iStep=0 ; iStep<pSteps ; iStep++ )
    {
		// Create new table of normals
		vector<Vector3D> lNewNormals( pMesh.mNormals.size() );

        for( size_t n=0 ; n<pMesh.mNormals.size() ; n++ )
		{
			lNewNormals[n] = pMesh.mNormals[n];

			BrowseList( iN, lNeighVers[n] )
				lNewNormals[n] += pMesh.mNormals[ lNeighVers[n].GetAt(iN) ];

			try
			{
				lNewNormals[n].Normalize();
			}
			catch(...)
			{
				// Reset to previous normal
				lNewNormals[n] = pMesh.mNormals[n];
			}
		}

        // Commit
        pMesh.mNormals = lNewNormals;
	}

    // Need to recompute DL
    pMesh.ResetDisplayList();
}
#pragma endregion
#pragma endregion


}
