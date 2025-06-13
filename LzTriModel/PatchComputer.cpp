#include "PatchComputer.h"
#include "Mesh.h"
#include <LzGeom/Line3D.h>
#include <LzServices/List.h>
#include <LzServices/Vector.h>


namespace LzTriModel
{
using LzGeom::Line3D;


//================================================================================
void PatchComputer::SplitMesh( const Mesh & pMesh,
                               List<Mesh> & pPatches,
                               bool pRemoveUnusedversNors,
                               const List<size_t> * ppHermeticEdges/*=nullptr*/ )
{
	// Compute topology
	MeshTopology lTopo;
	lTopo.Set( pMesh );

	// Split mesh
    SplitMesh( pMesh, lTopo, pPatches, pRemoveUnusedversNors, ppHermeticEdges );
}

//================================================================================
void PatchComputer::SplitMesh( const Mesh & pMesh,
                               const MeshTopology & pTopo,
                               List<Mesh> & pPatches,
                               bool pRemoveUnusedVersNors,
                               const List<size_t> * ppHermeticEdges/*=nullptr*/ )
{
    // Log
    LzLogNode("", "Splitting mesh into patches...")

	// Clear previous patch list
	pPatches.DelAll();

	// Mark all triangles as not processed
    vector<bool> lDone;
    lDone.assign( pTopo.TopTris().size(), false );

    // Extract patches from all triangles
    for( size_t iTopTri=0 ; iTopTri<lDone.size() ; iTopTri++ )
	{
		// Already processed?
        if( lDone[iTopTri] )
			continue;

		// Extract patch from triangle index
        List<size_t> lPatchIdx;
        ExtractFrom( lDone, pTopo, iTopTri, &lPatchIdx, ppHermeticEdges );

        // Stash empty mesh
        pPatches.AddTail( Mesh() );

        // In-place form mesh
        Mesh & lPatch = pPatches.GetTail();
		lPatch.mVertices = pMesh.mVertices;
		lPatch.mTexCoords2D = pMesh.mTexCoords2D;
		lPatch.mNormals = pMesh.mNormals;
		BrowseList( iIdx, lPatchIdx )
			lPatch.mTriangles.push_back( pMesh.mTriangles[ lPatchIdx.GetAt(iIdx) ] );

        // Clean unused?
        if( pRemoveUnusedVersNors )
            lPatch.RemoveUnusedVerticesAndNormals();

		// Log
        LzLogM("", "Found patch: "<<lPatch.mTriangles.size()<<" triangles, "<<lPatch.mVertices.size()<<" vertices, "<<lPatch.mNormals.size()<<" normals.")
	}

	// Log
    LzLogM("", "Found "<<pPatches.Count()<<" patch(es).")
}

//================================================================================
void * PatchComputer::FindPosOfMaxArea( const List<Mesh> & pPatches )
{
    void * lMaxPos = nullptr;
    double lMaxArea = 0.0;

    BrowseList( iP, pPatches )
    {
        // Compute area of patch
        const double lArea = pPatches.GetAt(iP).Area();

        // Update max?
        if( lMaxArea < lArea )
        {
            lMaxArea = lArea;
            lMaxPos = iP;
        }
    }

    return lMaxPos;
}

//===========================================================   =====================
size_t PatchComputer::CountPatches( const Mesh & pMesh,
                                    const List<size_t> * ppHermeticEdges/*=nullptr*/ )
{
    // Compute topology
    MeshTopology lTopo;
    lTopo.Set( pMesh );

    // Split mesh
    return CountPatches( lTopo, ppHermeticEdges );
}

//================================================================================
size_t PatchComputer::CountPatches( const MeshTopology & pTopo,
                                    const List<size_t> * ppHermeticEdges/*=nullptr*/ )
{
    // Log
    LzLogNode("", "Counting patches in mesh...")

    // Reset counter
    size_t lNbPatches = 0;

    // Mark all triangles as not processed
    vector<bool> lDone;
    lDone.assign( pTopo.TopTris().size(), false );

    // Extract patches from all triangles
    for( size_t iTopTri=0 ; iTopTri<lDone.size() ; iTopTri++ )
    {
        // Already processed?
        if( lDone[iTopTri] )
            continue;

        // Extract patch from triangle index
        ExtractFrom( lDone, pTopo, iTopTri, nullptr/* no need to collect tri idx */, ppHermeticEdges );

        // Update number of patches
        lNbPatches++;
    }

    // Log
    LzLogM("", "Counted "<<lNbPatches<<" patch(es).")

    // Finish
    return lNbPatches;
}

/*
//================================================================================
bool MeshTopology::CheckPatchOrientation( const List<size_t> & pPatch ) const
{
	// Empty patch?
	if( !pPatch.Count() )
		TxLogException("Cannot check orientation in an empty patch!");

	// To avoid pushing same triangle twice
	vector<bool> lPushed;
	lPushed.assign( mTopTris.size(), false );

	// Init
    List<size_t> lStack;
	lStack.AddTail( pPatch.GetHead() );
	lPushed[pPatch.GetHead()] = true;

	do
	{
		// Pop
        size_t lTriIdx = lStack.GetHead();
		lStack.DelHead();

		// Check
		if( lTriIdx >= mTopTris.size() )
			TxLogException("Triangle index "+lTriIdx+" is out of range! Have only "+mTopTris.size()+".");

		// Get current triangle
		const TopTri & lCurrTopT = mTopTris[lTriIdx];

		// Check all 3 neighbors
		for( int e=0 ; e<3 ; e++ )
		{
			int lEdgeIdx = lCurrTopT.mE[e];
			const TopEdge & lTopE = mTopEdges[lEdgeIdx];

//// Check if not multiedge
//if( lTopE.mT.Count() > 2 )
//	TxLogException("Found a multiedge in this patch! (index= "+lEdgeIdx+")");

			// Regular edge
			if( lTopE.mT.Count() == 2 )
			{
				int lOtherTriIdx = lTopE.OtherTri(lTriIdx);

				// Already pushed?
				if( lPushed[lOtherTriIdx] )
					continue;

				// Check orientation
				if( !lCurrTopT.HasSameOrientationAs( mTopTris[lOtherTriIdx] ) )
					return false;
				else
				{
					lStack.AddTail( lOtherTriIdx );
					lPushed[lOtherTriIdx] = true;
				}
			}
		}
	}
	while( lStack.Count() );

	// Everything went pretty well
	return true;
}
*/

//================================================================================
void PatchComputer::SplitMesh_HermeticME( const Mesh & pMesh,
                                          List< List<size_t> > & pPatchIdx )
{
	// Compute topology
	MeshTopology lTopo;
	lTopo.Set( pMesh );

	// Split mesh
    SplitMesh_HermeticME( lTopo, pPatchIdx );
}

//================================================================================
void PatchComputer::SplitMesh_HermeticME( const MeshTopology & pTopo,
                                          List< List<size_t> > & pPatchIdx )
{
    // Log
    LzLogNode("", "Splitting mesh into patches considering hermetic multi-edges...")

	// Clear previous patch list
	pPatchIdx.DelAll();

    // Mark all triangles as not processed
    vector<bool> lDone;
    lDone.assign( pTopo.TopTris().size(), false );

    for( size_t iTopTri=0 ; iTopTri<lDone.size() ; iTopTri++ )
	{
		// Already processed?
        if( lDone[iTopTri] )
			continue;

        // Stash empty list
        pPatchIdx.AddTail( List<size_t>() );

        // In-place extract patch from triangle index
        ExtractFrom_HermeticME( lDone, pTopo, iTopTri, pPatchIdx.GetTail() );

		// Log
        LzLogM("", "Found patch: "<<pPatchIdx.GetTail().Count()<<" triangle(s).")
	}

	// Log
    LzLogM("", "Found "<<pPatchIdx.Count()<<" patch(es).")
}

//================================================================================
void PatchComputer::ExtractPatch( const Mesh & pFrom,
                                  const List<size_t> & pTris,
                                  Mesh & pTo,
                                  bool pRemoveUnusedversNors )
{
	// Copy vertices, tex coords, and normals 
	pTo.mVertices = pFrom.mVertices;
	pTo.mTexCoords2D = pFrom.mTexCoords2D;
	pTo.mNormals = pFrom.mNormals;

	// Copy patch triangles only
	pTo.mTriangles.resize( pTris.Count() );
	{
        size_t t = 0;
		BrowseList( iT, pTris )
			pTo.mTriangles[ t++ ] = pFrom.mTriangles[ pTris.GetAt(iT) ];
	}

    // Clean unused?
    if( pRemoveUnusedversNors )
        pTo.RemoveUnusedVerticesAndNormals();

    // Display list
	pTo.ResetDisplayList();
}

//================================================================================
void PatchComputer::ExtractFrom( vector<bool> & pDone,
                                 const MeshTopology & pTopo,
                                 size_t pIdx,
                                 List<size_t> * ppPatchIdx,
                                 const List<size_t> * ppHermeticEdges )
{
	// Push
    List<size_t> lStack;
    lStack.AddTail( pIdx );
    pDone[pIdx] = true;

    while( lStack.Count() )
	{
		// Pop triangle index
        const size_t lSrcTri = lStack.GetHead();
        lStack.DelHead();

		// Store in patch
        if( ppPatchIdx )
            ppPatchIdx->AddTail( lSrcTri );

		// Current triangle
		const TopTri & lTopTri = pTopo.TopTris()[ lSrcTri ];

		// Check all edges
		for( int e=0 ; e<3 ; e++ )
		{
			// The index
			const int lEdgeIdx = lTopTri.mE[e];

			// Check if forbidden edge
			if( ppHermeticEdges )
			{
                // By default, yes, consider this edge
                bool lConsiderThisEdge = true;
				BrowseList( iHE, *ppHermeticEdges )
				{
                    if( lEdgeIdx == static_cast<int>(ppHermeticEdges->GetAt(iHE)) )
					{
                        lConsiderThisEdge = false;
						break;
					}
				}

				// Skip?
                if( !lConsiderThisEdge )
					continue;
			}

			// The edge
			const TopEdge & lEdge = pTopo.TopEdges()[ lEdgeIdx ];

			// Push triangles
			BrowseList( iT, lEdge.mT )
			{
                const size_t lOtherTri = lEdge.mT.GetAt( iT );

                // Skip already processed triangles (including lSrcTri, already marked as pDone)
                if( pDone[lOtherTri] )
                    continue;

				// Push
                lStack.AddTail( lOtherTri );
                pDone[ lOtherTri ] = true;
			}
		}
	}
}

//================================================================================
void PatchComputer::ExtractFrom_HermeticME( vector<bool> & pDone,
                                            const MeshTopology & pTopo,
                                            size_t pIdx,
                                            List<size_t> & pPatchIdx )
{
    // Push
    List<size_t> lStack;
    lStack.AddTail( pIdx );
    pDone[pIdx] = true;

    while( lStack.Count() )
	{
		// Pop triangle index
        const size_t lSrcTri = lStack.GetHead();
        lStack.DelHead();

		// Store in patch
		pPatchIdx.AddTail( lSrcTri );

		// Current triangle
		const TopTri & lTopTri = pTopo.TopTris()[ lSrcTri ];

		// Check all edges
		for( int e=0 ; e<3 ; e++ )
		{
			// The index
			const int lEdgeIdx = lTopTri.mE[e];

			// Check if forbidden edge
			if( pTopo.TopEdges()[ lEdgeIdx ].mT.Count() > 2 )
				continue;

			// The edge
			const TopEdge & lEdge = pTopo.TopEdges()[ lEdgeIdx ];

			// Push triangles
			BrowseList( iT, lEdge.mT )
			{
                size_t lOtherTri = lEdge.mT.GetAt( iT );

                // Skip already processed triangles (including lSrcTri, already marked as pDone)
                if( pDone[lOtherTri] )
                    continue;

				// Push
                lStack.AddTail( lOtherTri );
                pDone[ lOtherTri ] = true;
			}
		}
	}
}
}
