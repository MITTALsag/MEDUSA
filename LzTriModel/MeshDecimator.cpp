#include "MeshDecimator.h"
#include "Triangle.h"
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>
//#include <LzOpenGL/gl.h>
//#include <float.h>


namespace LzTriModel
{
#pragma region "EnhancedVertex"
//================================================================================
EnhancedVertex::EnhancedVertex()
{
    mIdx = -1;
    mCollapseTgt = -1;
    mCollapseCost = std::numeric_limits<double>::max();
}

//================================================================================
void EnhancedVertex::Log() const
{
//#ifdef USE_CLI
    // Triangles
	std::string lTrianglesLog;
	{
		std::stringstream lStrStr;
		lStrStr << "{ ";
//String ^lTrianglesLog = "{ ";
//{
        BrowseList(iT, mTriangles)
			lStrStr << mTriangles.GetAt(iT) << " ";
//{
//    lTrianglesLog += mTriangles.GetAt(iT);
//    lTrianglesLog += " ";
//}
		lStrStr << "}";
//lTrianglesLog += "}";

		lTrianglesLog = lStrStr.str();
    }

    // Neighbors
	std::string lNeighborsLog;
	{
		std::stringstream lStrStr;
		lStrStr << "{ ";
//String ^lNeighborsLog = "{ ";
//{
        BrowseList(iN, mNeighbors)
			lStrStr << mNeighbors.GetAt(iN) << " ";
//{
//    lNeighborsLog += mNeighbors.GetAt(iN);
//    lNeighborsLog += " ";
//}
		lStrStr << "}";
//lNeighborsLog += "}";

		lNeighborsLog = lStrStr.str();
    }

    // Log all
    {
        LzLogN("", "Enhanced Vertex "<<mIdx)
        LzLogM("", "Valid: "<<IsValid())
        LzLogM("", "Triangles: "<<lTrianglesLog)
        LzLogM("", "Neighbors: "<<lNeighborsLog)
        LzLogM("", "Best collapse: ["<<mIdx<<" -> "<<mCollapseTgt<<"] cost: "<<mCollapseCost)
    }
        
//#else
//        TxLogException("BOUYYYYYA");
//#endif
}
#pragma endregion


#pragma region "Construction & destruction"
//================================================================================
MeshDecimator::~MeshDecimator()
{
	Free();
}

//================================================================================
void MeshDecimator::Free()
{
    mEV.clear();
    mET.clear();
    mCostMin = -1;
    mCostMax = -1;
}
#pragma endregion


#pragma region "Decimate - GameDeveloper algo 1998"
//================================================================================
void MeshDecimator::Decimate_GameDeveloper1998( Mesh & pMesh, const MeshTopology & pTopo, unsigned int pTgtTrianglesNo, bool pPreserveCracks/*=false*/, double pMaxNormalAngles_Deg/*=180*/ )
{
    _LzLogTry

    LzLogNode("", "Decimating Mesh")

    // First of all
    Free();

    // Log decimation parameters
    if( pPreserveCracks)
        LzLogM("", "Preserving cracks.")

    if( pMaxNormalAngles_Deg < 180 )
        LzLogM("", "Preserving topology (normal angles >= " << pMaxNormalAngles_Deg << " degrees).")

    // Watch decimation performance
    auto lStartDecim = LzServices::Chrono_Now();

    // Check targeted triangles number < initial mesh's triangles.
    if( pTgtTrianglesNo > pMesh.mTriangles.size() )
    {
        LzLogException("", "Can't decimate to more triangles than the mesh's initial number. Aborting.")
        return;
    }


#pragma region "A. Initialize structures"
    // Topology info
    {
        // Update 
        mEV.resize( pMesh.mVertices.size() );
        mET.resize( pMesh.mTriangles.size() );

        for ( unsigned int iT=0; iT<pMesh.mTriangles.size(); iT++ )
        {
            // Enhanced triangles 
            {
                mET[iT].mIdx = iT;
                mET[iT].mNormal = ( pMesh.mNormals[ pMesh.mTriangles[iT].mIdxN[0] ]
                                    +  pMesh.mNormals[ pMesh.mTriangles[iT].mIdxN[1] ]
                                    +  pMesh.mNormals[ pMesh.mTriangles[iT].mIdxN[2] ] ) / 3;
            }

            // Enhanced vertices
            {
                // Triangles' vertices
                const int & lIdxV0 = pMesh.mTriangles[iT].mIdxV[0];
                const int & lIdxV1 = pMesh.mTriangles[iT].mIdxV[1];
                const int & lIdxV2 = pMesh.mTriangles[iT].mIdxV[2];

                // Idx
                mEV[ lIdxV0 ].mIdx = lIdxV0;
                mEV[ lIdxV1 ].mIdx = lIdxV1;
                mEV[ lIdxV2 ].mIdx = lIdxV2;

                // Store triangles (unique)
                mEV[ lIdxV0 ].mTriangles.AddIntoIncrList( iT, true );
                mEV[ lIdxV1 ].mTriangles.AddIntoIncrList( iT, true );
                mEV[ lIdxV2 ].mTriangles.AddIntoIncrList( iT, true );

                // Store neighbor vertices (unique)
                mEV[ lIdxV0 ].mNeighbors.AddIntoIncrList( lIdxV1, true );
                mEV[ lIdxV0 ].mNeighbors.AddIntoIncrList( lIdxV2, true );
                mEV[ lIdxV1 ].mNeighbors.AddIntoIncrList( lIdxV0, true );
                mEV[ lIdxV1 ].mNeighbors.AddIntoIncrList( lIdxV2, true );
                mEV[ lIdxV2 ].mNeighbors.AddIntoIncrList( lIdxV0, true );
                mEV[ lIdxV2 ].mNeighbors.AddIntoIncrList( lIdxV1, true );
            }
        }
    }
            
    // Collapse info
    {
        // Get minimum collapse info for each vertex
        for( unsigned int iIdxVOrg=0 ; iIdxVOrg<mEV.size() ; iIdxVOrg++ )
            ComputeBestCost( pMesh, pTopo, mEV[iIdxVOrg], pPreserveCracks, pMaxNormalAngles_Deg );
    }
#pragma endregion 


#pragma region "B. Decimation loop"
    // Number of Triangles to remove
    unsigned int lTNoToRemove = (unsigned int)(pMesh.mTriangles.size() - pTgtTrianglesNo);

    // Store triangles to remove, remove them once lTNoToRemove has been reached
    List<int> lTtoRemove;

    // Check triangles flip, i.e. before & after decimation normals' angle difference > lFlipTolerance 
    double lFlipTolerance_deg = 120.0; 
    double lFlipTolerance_rad = LzMath::DEG_2_RAD * lFlipTolerance_deg ;

    // Main loop: perform edge collapse until targeted number of triangles has been reached
    while( lTtoRemove.Count() < lTNoToRemove )
    {
        // Find best edge to collapse
        EnhancedVertex lBestEV;
        for( unsigned int iEV=0 ; iEV<mEV.size() ; iEV++)
        {
            if ( mEV[ iEV ].IsValid() && lBestEV.mCollapseCost >= mEV[ iEV ].mCollapseCost )
                lBestEV = mEV[ iEV ];
        }

        // Performing edge collapse lBestEV.mIdx -> lBestEV.mCollapseTgt
        {
            // 1. Neighborhood update
            BrowseList(iN, lBestEV.mNeighbors)
            {
                const unsigned int & lIdxN = lBestEV.mNeighbors.GetAt(iN);

                // a. Tell my neighbors not to link with me anymore
                mEV[ lIdxN ].mNeighbors.DelElts( lBestEV.mIdx );

                // b. Tell my neighbors to link with mCollapseTgt
                if( lIdxN != lBestEV.mCollapseTgt)
                    mEV[ lIdxN ].mNeighbors.AddIntoIncrList( lBestEV.mCollapseTgt, true );

                // c. Tell mCollapseTgt to link with all my neighbors
                if( lIdxN != lBestEV.mCollapseTgt)
                    mEV[ lBestEV.mCollapseTgt ].mNeighbors.AddIntoIncrList(lIdxN, true);
            }
         
            // 2. Triangle updates
            BrowseList(iT, lBestEV.mTriangles)
            {
                const int & lIdxT = lBestEV.mTriangles.GetAt(iT);
                Triangle & lT = pMesh.mTriangles[ lIdxT ];

                // a. Remove triangles sharing mIdx & mCollapseTgt
                {
                    if ( lT.mIdxV[0] == lBestEV.mCollapseTgt || 
                            lT.mIdxV[1] == lBestEV.mCollapseTgt || 
                            lT.mIdxV[2] == lBestEV.mCollapseTgt )
                    {
                        // Remove lT reference from my EV neighbors
                        BrowseList(iN, lBestEV.mNeighbors)
                        {
                            const unsigned int & lIdxN = lBestEV.mNeighbors.GetAt(iN);
                            mEV[lIdxN].mTriangles.DelElts(lIdxT);
                        }

                        // Tell T has to be removed from mesh structure
                        lTtoRemove.AddIntoDecrList( lIdxT, false );

                        // Invalidate T in the enhanced triangles structure
                        mET[lIdxT].mIdx = -1;
                    }
                }

                // b. Update triangles' vertices in the neighborhood
                {
                    if ( lT.mIdxV[0] != lBestEV.mCollapseTgt && 
                            lT.mIdxV[1] != lBestEV.mCollapseTgt && 
                            lT.mIdxV[2] != lBestEV.mCollapseTgt )
                    {
                        // Vertex update mIdx -> mCollapseTgt
                        if ( lT.mIdxV[0] == lBestEV.mIdx)
                                lT.mIdxV[0]  = lBestEV.mCollapseTgt;
                        if ( lT.mIdxV[1] == lBestEV.mIdx)
                                lT.mIdxV[1]  = lBestEV.mCollapseTgt;
                        if ( lT.mIdxV[2] == lBestEV.mIdx)
                                lT.mIdxV[2]  = lBestEV.mCollapseTgt;

                        // Add this triangle to mCollapseTgt's ones
                        mEV[ lBestEV.mCollapseTgt ].mTriangles.AddIntoIncrList(lIdxT, true);
                    }
                }
            }
                
            // 3. Recompute all neighbors cost
            {
                List<int> lAllNeighbors;
                lAllNeighbors.AddIntoIncrList( lBestEV.mNeighbors, true );
                lAllNeighbors.AddIntoIncrList( mEV[ lBestEV.mCollapseTgt ].mNeighbors, true );
                BrowseList(iN, lAllNeighbors)
                {
                    const unsigned int & lIdxN = lAllNeighbors.GetAt(iN);

                    // Recompute cost
                    ComputeBestCost( pMesh, pTopo, mEV[ lIdxN ], pPreserveCracks, pMaxNormalAngles_Deg );
                }
            }

            // 4. Update all neighbor triangles' normals
            {
                List<int> lAllTriangles;
                lAllTriangles.AddIntoIncrList( lBestEV.mTriangles, true );
                lAllTriangles.AddIntoIncrList( mEV[ lBestEV.mCollapseTgt ].mTriangles, true );
                BrowseList(iT, lAllTriangles)
                {
                    const unsigned int & lIdxT = lAllTriangles.GetAt(iT);

                    // Save previous normal
                    mET[ lIdxT ].mLastNormal = mET[ lIdxT ].mNormal;

                    // Compute current normals
                    const Triangle & lT = pMesh.mTriangles[ lIdxT ];
                    const Point3D & lT_Pt0 = pMesh.mVertices[ lT.mIdxV[0] ];
                    const Point3D & lT_Pt1 = pMesh.mVertices[ lT.mIdxV[1] ];
                    const Point3D & lT_Pt2 = pMesh.mVertices[ lT.mIdxV[2] ];
//************ throws exception on _MDL Samples_/CubeBig.mdl
//************ throws exception on _MDL Samples_/CubeBig.mdl
Plane3D lTPlane;
try
{
	lTPlane = Plane3D( lT_Pt0, lT_Pt1, lT_Pt2 );
}
catch( ... )
{
    LzLogErr("", "Could not build plane from points")
    LzLogErr("", "Pt0= "<<lT_Pt0.ToString())
    LzLogErr("", "Pt1= "<<lT_Pt1.ToString())
    LzLogErr("", "Pt2= "<<lT_Pt2.ToString())

	Vector3D lNormal = ( lT_Pt1 - lT_Pt0 ) ^ ( lT_Pt2 - lT_Pt0 );
    LzLogErr("", "Norm= "<<lNormal.Norm())

	lTPlane = Plane3D();
}

//const Plane3D lTPlane( lT_Pt0, lT_Pt1, lT_Pt2 );
//************ throws exception on _MDL Samples_/CubeBig.mdl
//************ throws exception on _MDL Samples_/CubeBig.mdl
                    mET[ lIdxT ].mNormal = lTPlane.Normal();

                    // Check this triangle did not flip too much
                    {
                        double lAngle = mET[ lIdxT ].mLastNormal.ShortPositiveAngleTo( mET[ lIdxT ].mNormal );
                        if( lAngle > lFlipTolerance_rad )
                            mET[ lIdxT ].mHasFlipped = true;
                    }
                }
            }

            // 5. Collapse to the middle of the edge [VOrg, VTgt]
            if( !pTopo.HasCrackVer(lBestEV.mIdx) && !pTopo.HasCrackVer(lBestEV.mCollapseTgt) )
            {
                Point3D & lTgtPt = pMesh.mVertices[ lBestEV.mCollapseTgt ];
                {
                    // vertex position
                    const Point3D & lOrgPt = pMesh.mVertices[ lBestEV.mIdx ];
                    lTgtPt = lTgtPt.MidPointTo(lOrgPt);
                }

            }
         
            // 6. Invalidate lBestEV: do not consider this vertex anymore
            mEV[ lBestEV.mIdx ].mIdx = -1;
        }
    }
#pragma endregion 


#pragma region "C. Clean-up structures"
    // 1. Delete triangles of the Mesh data structure
    BrowseList(iT, lTtoRemove)
    {
        const int & lIdxT = lTtoRemove.GetAt(iT);
        pMesh.mTriangles.erase( pMesh.mTriangles.begin() + lIdxT );
    }

    // 2. Update all normals, advice user if some T flipped during decimation
    {
        pMesh.mNormals.clear();
        pMesh.mNormals.resize( pMesh.mTriangles.size() );
        unsigned int lNormalCpt = 0;
        unsigned int lNoTFlipped = 0;
        for( unsigned int iET=0; iET<mET.size(); iET++ )
        {
            if( !mET[ iET ].IsValid() )
                continue;

            // Save the actual normal in the Mesh datastructure
            pMesh.mNormals[ lNormalCpt ] = mET[ iET ].mNormal;
            pMesh.mTriangles[ lNormalCpt ].mIdxN[0] = lNormalCpt;
            pMesh.mTriangles[ lNormalCpt ].mIdxN[1] = lNormalCpt;
            pMesh.mTriangles[ lNormalCpt ].mIdxN[2] = lNormalCpt;
            lNormalCpt++;

            // Associated valid T has flipped?
            if( mET[ iET ].HasFlipped() )
                lNoTFlipped;
        }

        // Advice user if some triangles have flipped during decimation
        if( lNoTFlipped > 0 )
            LzLogException("", lNoTFlipped << " triangle(s) have flipped during decimation! \n" <<
                                    "Inverted triangles appear in green on the mesh surface.\n\n" <<
                                    "You might want to run a laplacian filter to clean them.")
    }

    // 3. Delete vertices of the Mesh data structure (normals are up to date)
    pMesh.RemoveUnusedVerticesAndNormals();
#pragma endregion


    // Log
    auto lDecimTime = LzServices::Chrono_MSecs(lStartDecim);
    LzLogMsg("", "Removed "<<lTtoRemove.Count()<<" triangles in "<<lDecimTime<<" msec.")

    _LzLogCatchAndThrow("Error while decimating the mesh.")
}

//================================================================================
double MeshDecimator::ComputeCost( const Mesh & pMesh, const MeshTopology & pTopo,
                                   const EnhancedVertex & pVOrg, const EnhancedVertex & pVTgt,
                                   bool pPreserveCracks, double pMaxNormalAngles_Deg )
{
    // Preserve cracks?
    if( pPreserveCracks &&
      (  (pTopo.HasCrackVer( pVOrg.mIdx ) && !pTopo.HasCrackVer( pVTgt.mIdx ))
      || (pTopo.HasCrackVer( pVTgt.mIdx ) && !pTopo.HasCrackVer( pVOrg.mIdx ))
       ) ) // This rule allows collapsing if pVOrg & pVTgt are BOTH cracks.
    {
        return std::numeric_limits<double>::max();
    }

    // Avoid creating multi-edges
    {
        // Using connectivity criterion: 
        // http://stackoverflow.com/questions/27049163/mesh-simplification-edge-collapse-conditions
        // Check card(Nu intersect Nv / {u, v}) == 2 (particular case: 1 if u - v is a crack)
        List<int> lCommonNeighbors;
        BrowseList(iNVOrg, pVOrg.mNeighbors)
        {
            const int & lNVorg = pVOrg.mNeighbors.GetAt(iNVOrg);
            if( pVTgt.mNeighbors.FindElt(lNVorg) )
                lCommonNeighbors.AddTail(lNVorg);
        }

        // Crack collapse
        if ( pTopo.HasCrackVer( pVOrg.mIdx ) && pTopo.HasCrackVer( pVTgt.mIdx ) &&
                lCommonNeighbors.Count() != 1 ) // expected connexity for collapsing a crack
            return std::numeric_limits<double>::max();
            
        // Other collapses
        if ( (!pTopo.HasCrackVer( pVOrg.mIdx ) || !pTopo.HasCrackVer( pVTgt.mIdx )) && 
                lCommonNeighbors.Count() != 2 ) 
            return std::numeric_limits<double>::max();
    } 

    // Edge [pVOrg, pVTgt]'s length
    double lEdgeLgth = pMesh.mVertices[ pVOrg.mIdx ].DistanceTo( pMesh.mVertices[ pVTgt.mIdx ] );

    // Find triangles having pVorg & pVTgt as vertices
    List<int> lSharedTriangles;
    {
        BrowseList(iT, pVOrg.mTriangles)
        {
            const int & lTIdx = pVOrg.mTriangles.GetAt(iT);
            if( pVTgt.mTriangles.FindElt(lTIdx) )
                lSharedTriangles.AddTail(lTIdx);
        }
    }

    // Find maximum cost over all triangles sharing pVOrg
    double lFinalCost = 0.0;
    {
        BrowseList(iT, pVOrg.mTriangles)
        {
            const Triangle & lT =  pMesh.mTriangles[ pVOrg.mTriangles.GetAt(iT) ];
            Vector3D lTN = ( pMesh.mNormals[ lT.mIdxN[0] ] + 
                                pMesh.mNormals[ lT.mIdxN[1] ] + 
                                pMesh.mNormals[ lT.mIdxN[2] ] ) / 3.0;

            // Sometime normals are not normalized ...
            lTN.Normalize();

            // Find minimum cost over triangles having pVorg & pVTgt as vertices
            double lCostMin = std::numeric_limits<double>::max();
            BrowseList(iT_shared, lSharedTriangles)
            {
                const Triangle & lT_shared =  pMesh.mTriangles[ lSharedTriangles.GetAt(iT_shared) ];
                Vector3D lTN_shared = ( pMesh.mNormals[ lT_shared.mIdxN[0] ] + 
                                        pMesh.mNormals[ lT_shared.mIdxN[1] ] + 
                                        pMesh.mNormals[ lT_shared.mIdxN[2] ] ) / 3.0;

                // Sometime normals are not normalized ...
                lTN_shared.Normalize();

// NICO Tweak cost function such as coplanary triangles => lCost is still > 0
// double lCost = ( 1 - (lTN * lTN_shared) ) / 2.0;
// double lCost = exp( (1 - (lTN * lTN_shared)) / 2.0 );
double lCost = (2 - (lTN * lTN_shared)) / 2.0;

// Preserve topology A FIGNOLER parce que plus compliqué que prévu (dépend du collapse au milieu ou pas)
#if 0
                {
                    double lNormalsAngle = TxMath::ToolBox::Rad2Deg( lTN.PositiveAngleTo( lTN_shared, lTN^lTN_shared ) );
                
                    if ( lNormalsAngle >= pMaxNormalAngles_Deg)
                    {
                        TxLogM("Normal angle = "+lNormalsAngle);
                        lCost = std::numeric_limits<double>::max();
                    }
                }
#endif

                if (lCostMin > lCost)
                    lCostMin = lCost;
            }

            // Take the maximum amoung all triangles sharing pVOrg
            if( lFinalCost < lCostMin )
                lFinalCost = lCostMin;
        }
    }

    return lEdgeLgth * lFinalCost;
}

//================================================================================
void MeshDecimator::ComputeBestCost(const Mesh & pMesh, const MeshTopology & pTopo, EnhancedVertex & pV, bool pPreserseCracks, double pMaxNormalAngles_Deg)
{
    pV.mCollapseCost = std::numeric_limits<double>::max();
    pV.mCollapseTgt = -1;

    // Get minimum collapse cost over vertex's neighbors
    BrowseList(iVTgt, pV.mNeighbors)
    {
        const int & lIdxVTgt = pV.mNeighbors.GetAt(iVTgt);
        double lCost = ComputeCost( pMesh, pTopo, pV, mEV[lIdxVTgt], pPreserseCracks, pMaxNormalAngles_Deg );
        if( lCost <= pV.mCollapseCost )
        {
            pV.mCollapseCost = lCost;
            pV.mCollapseTgt  = lIdxVTgt;
        }
    }
}
#pragma endregion


#pragma region "Edge Collapse Cost"
//================================================================================
void MeshDecimator::ComputeEnhancedVertexStructure(const Mesh & pMesh, const MeshTopology & pTopo)
{
    _LzLogTry

    LzLogNode("", "Computing Enhanced Vertice structure.")

    // First of all
    Free();

    // By default, preserve cracks for cost computation
    const bool lPreserveCracks = true;

    // Topology info
    {
        // Update 
        mEV.resize( pMesh.mVertices.size() );

        for ( unsigned int iT=0; iT<pMesh.mTriangles.size(); iT++ )
        {
            // Triangles vertices
            const int & lIdxV0 = pMesh.mTriangles[iT].mIdxV[0];
            const int & lIdxV1 = pMesh.mTriangles[iT].mIdxV[1];
            const int & lIdxV2 = pMesh.mTriangles[iT].mIdxV[2];

            // Idx
            mEV[ lIdxV0 ].mIdx = lIdxV0;
            mEV[ lIdxV1 ].mIdx = lIdxV1;
            mEV[ lIdxV2 ].mIdx = lIdxV2;

            // Store triangles (unique)
            mEV[ lIdxV0 ].mTriangles.AddIntoIncrList( iT, true );
            mEV[ lIdxV1 ].mTriangles.AddIntoIncrList( iT, true );
            mEV[ lIdxV2 ].mTriangles.AddIntoIncrList( iT, true );

            // Store neighbor vertices (unique)
            mEV[ lIdxV0 ].mNeighbors.AddIntoIncrList( lIdxV1, true );
            mEV[ lIdxV0 ].mNeighbors.AddIntoIncrList( lIdxV2, true );
            mEV[ lIdxV1 ].mNeighbors.AddIntoIncrList( lIdxV0, true );
            mEV[ lIdxV1 ].mNeighbors.AddIntoIncrList( lIdxV2, true );
            mEV[ lIdxV2 ].mNeighbors.AddIntoIncrList( lIdxV0, true );
            mEV[ lIdxV2 ].mNeighbors.AddIntoIncrList( lIdxV1, true );
        }
    }

    // Collapse info
    {
        // Get minimum collapse info for each vertex
        for( unsigned int iIdxVOrg=0 ; iIdxVOrg<mEV.size() ; iIdxVOrg++ )
            ComputeBestCost( pMesh, pTopo, mEV[iIdxVOrg], lPreserveCracks, 180 );
    }

    // Get min and max costs
    {
        mCostMin = std::numeric_limits<double>::max();
        for( unsigned int iEV=0; iEV<mEV.size(); iEV++)
        {
            if( mCostMin > mEV[iEV].mCollapseCost )
                mCostMin = mEV[iEV].mCollapseCost;
            if( mCostMax < mEV[iEV].mCollapseCost && mEV[iEV].mCollapseCost != std::numeric_limits<double>::max() )
                mCostMax = mEV[iEV].mCollapseCost;
        }
        LzLogM("", "Cost min =" << mCostMin)
        LzLogM("", "Cost max =" << mCostMax)
    }
    
    _LzLogCatchAndThrow("Error while computing Enhanced vertex structure.")
}

#ifdef USE_CLI
//================================================================================
void MeshDecimator::DrawCollapseCostEdge( const Mesh & pMesh )
{
    glBegin( GL_LINES );
	for( unsigned int iEV=0 ; iEV < mEV.size() ; iEV++ )
	{
		const EnhancedVertex & lEV = mEV[iEV];

        if ( lEV.IsValid() )
        {
            TxMath::ToolBox::SetScaleColor( mCostMin, mCostMax, lEV.mCollapseCost );
			glVertex3dv( pMesh.mVertices[ lEV.mIdx ].mV );
			glVertex3dv( pMesh.mVertices[ lEV.mCollapseTgt ].mV );
        }
	}
	glEnd();
}
#endif
#pragma endregion


}
