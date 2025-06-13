#include "MeshTopology.h"
#include <LzServices/LzLog.h>
#include <LzServices/HashTable.h>


namespace LzTriModel
{
using LzServices::HashTable;


//********************************************************************************
//********************************************************************************
#pragma region "TopVer"
//================================================================================
TopVer::TopVer()
{
}

//================================================================================
TopVer::~TopVer()
{
}
#pragma endregion


//********************************************************************************
//********************************************************************************
#pragma region "TopTri"
//================================================================================
TopTri::TopTri( int pV0/*=-1*/, int pV1/*=-1*/, int pV2/*=-1*/ )
{
	mV[0] = pV0;
	mV[1] = pV1;
	mV[2] = pV2;

	mE[0] = mE[1] = mE[2] = -1;
}

//================================================================================
TopTri::~TopTri()
{
}

//================================================================================
int TopTri::EdgeWithVertex( int pIdxV, int pExcludeEdge/*=-1*/ ) const
{
	for( int v=0 ; v<3 ; v++ )
	{
		if( pIdxV == mV[v] )
		{
			int lE1 = mE[ (v+1) % 3 ];
			return lE1!=pExcludeEdge ? lE1 : mE[ (v+2) % 3 ] ;
		}
	}

    LzLogException("", "Vertex #"<<pIdxV<<" not found in this triangle!");
}

//================================================================================
int TopTri::EdgeWithoutVertex( int pIdxV ) const
{
	for( int v=0 ; v<3 ; v++ )
	{
		if( pIdxV == mV[v] )
			return mE[v];
	}

    LzLogException("", "Vertex #"<<pIdxV<<" not found in this triangle!");
}

//================================================================================
int TopTri::VertexOppositeToEdge( int pIdxE, int * pTriIdx/*=nullptr*/ ) const
{
	for( int v=0 ; v<3 ; v++ )
	{
		if( pIdxE == mE[v] )
		{
			if( pTriIdx )
				*pTriIdx = v;

			return mV[v];
		}
	}

    LzLogException("", "Edge #"<<pIdxE<<" not found in this triangle!");
}

//================================================================================
int TopTri::OtherEdge( int pExcludeEdge1, int pExcludeEdge2/*=-1*/ ) const
{
	for( int v=0 ; v<3 ; v++ )
	{
		if( mE[v]!=pExcludeEdge1 && mE[v]!=pExcludeEdge2 )
			return mE[v];
	}

    LzLogException("", "Edges are kinda mixed-up in this triangle!");
}

//================================================================================
bool TopTri::HasSameOrientationAs( const TopTri & pOther ) const
{
	int lOppV[3] = { mV[2], mV[1], mV[0] };

	// Check all 3 shifts
	for( int i=0 ; i<3 ; i++ )
	{
		int lMatch = 0;
		for( int j=0 ; j<3 ; j++ )
		{
			if( pOther.mV[j] == lOppV[(j+i)%3] )
				lMatch++;
		}

		if( lMatch == 2 )
			return true;
	}

	return false;
}

//================================================================================
bool TopTri::HasVertex( int pIdxV ) const
{
	if( mV[0]==pIdxV || mV[1]==pIdxV || mV[2]==pIdxV )
		return true;

	return false;
}
#pragma endregion


//********************************************************************************
//********************************************************************************
#pragma region "TopEdge"
//================================================================================
TopEdge::TopEdge( int pV0/*=-1*/, int pV1/*=-1*/ )
{
	mV[0] = pV0;
	mV[1] = pV1;
}

//================================================================================
TopEdge::~TopEdge()
{
}

//================================================================================
bool TopEdge::HasVertex( int pIdxV ) const
{
	return mV[0]==pIdxV || mV[1]==pIdxV;
}

//================================================================================
int TopEdge::OtherVer( int pExcludeVer ) const
{
	return mV[0]!=pExcludeVer ? mV[0] : mV[1] ;
}

//================================================================================
int TopEdge::OtherTri( int pExcludeTri ) const
{
	if( mT.Count() != 2 )
        LzLogException("", "Edge must have exactly 2 triangles in order to use this method!");

	if( mT.GetHead() != pExcludeTri )
		return mT.GetHead();
	else
		return mT.GetTail();
}

//================================================================================
bool TopEdge::IsConnected( const TopEdge & pEdge ) const
{
	return pEdge.HasVertex(mV[0]) || pEdge.HasVertex(mV[1]);
}
#pragma endregion


//********************************************************************************
//********************************************************************************
#pragma region "MeshTopology"

#pragma region "Construction & destruction"
//================================================================================
MeshTopology::MeshTopology()
{
}

//================================================================================
MeshTopology::MeshTopology( std::function<void(MeshTopology * pThis)> pInit )
{
    // Initialize
    pInit(this);
}

//================================================================================
MeshTopology::~MeshTopology()
{
	Free();
}

//================================================================================
void MeshTopology::Free()
{
//mpVertices = nullptr;
//mpTexCoords2D = nullptr;
//mpNormals = nullptr;
//mpTriangles = nullptr;

	mTopVers.clear();
	mTopTris.clear();
	mTopEdges.clear();
	mCracks.clear();
	mMultiEdges.clear();
}
#pragma endregion


#pragma region "Find elements"
//================================================================================
bool MeshTopology::HasCrackVer( size_t pVer ) const
{
    for( size_t i=0 ; i<mCracks.size() ; i++ )
    {
        if( mTopEdges[ mCracks[i] ].mV[0] == (int)pVer
         || mTopEdges[ mCracks[i] ].mV[1] == (int)pVer )
			return true;
    }

	return false;
}
#pragma endregion


#pragma region "Analysis"
//================================================================================
void MeshTopology::Set( const Mesh & pMesh )
{
try
{
	// Free previous topology
	Free();

	//--------------------------------------
    // Vertices
	//--------------------------------------

	mTopVers.resize( pMesh.mVertices.size() );

	//--------------------------------------
    // Triangles
	//--------------------------------------

	mTopTris.resize( pMesh.mTriangles.size() );

	for( vector<Triangle>::const_iterator iTri=pMesh.mTriangles.begin() ; iTri!=pMesh.mTriangles.end() ; iTri++ )
	{
		// Index of current triangle
		int lTriIdx = (int)(iTri - pMesh.mTriangles.begin());

		// Create new triangle
		TopTri lTopTri( iTri->mIdxV[0], iTri->mIdxV[1], iTri->mIdxV[2] );

		// Check triangle semantics
		for( int i=0 ; i<3 ; i++ )
		{
			if( lTopTri.mV[i] == lTopTri.mV[(i+1)%3] )
                LzLogException("", "Found degenerate triangle #"<<lTriIdx<<"! Duplicated index= "<<lTopTri.mV[i]<<".");
		}

// Store triangle
//mTopTris.push_back( lTopTri );
mTopTris[lTriIdx] = lTopTri;

		// Register triangle for all vertices
		for( int i=0 ; i<3 ; i++ )
			mTopVers[ lTopTri.mV[i] ].mT.AddTail( lTriIdx );
	}

	//--------------------------------------
	// Edges
	//--------------------------------------

	// Fast edge analysis
	HashTable<TopEdge> lHashEdges(10000/*carefully chosen magic number*/);

#ifdef _DEBUG
    size_t lMinKey = std::numeric_limits<size_t>::max();
    size_t lMaxKey = 0;
#endif

    // Create edges
	for( vector<TopTri>::const_iterator iTopTri=mTopTris.begin() ; iTopTri!=mTopTris.end() ; iTopTri++ )
	{
		// Index of current triangle
		int lTopTriIdx = (int)(iTopTri - mTopTris.begin());

		// For each new edge
		for( int i=0 ; i<3 ; i++ )
		{
			// Edge vertices
			int lV0 = iTopTri->mV[(i+0)%3];
			int lV1 = iTopTri->mV[(i+1)%3];

			// Try to find edge
            size_t lEdgeCode = HashTable<TopEdge>::HashKey( lV0 + lV1 );
#ifdef _DEBUG
//    LzLogM("", "Code= "<<lEdgeCode)
//    LzLogM("", "lMinKey= "<<lMinKey)
//    LzLogM("", "lMaxKey= "<<lMaxKey)

    if( lMinKey > lEdgeCode )
    {
//        LzLogM("", "Code= "<<lEdgeCode<<", min key= "<<lMinKey)
        lMinKey = lEdgeCode;
//        LzLogM("", "New min key= "<<lMinKey)
    }
    if( lMaxKey < lEdgeCode )
    {
//        LzLogM("", "Code= "<<lEdgeCode<<", max key= "<<lMaxKey)
        lMaxKey = lEdgeCode;
//        LzLogM("", "New max key= "<<lMaxKey)
    }
//    LzLogM("", "")
#endif
            //
            TopEdge * lpFoundEdge = nullptr;
			List<TopEdge> & lEdgeList = lHashEdges.ListOf(lEdgeCode);
			BrowseList( iEdge, lEdgeList )
			{
				TopEdge & lEdge = lEdgeList.GetAt(iEdge);
				if( lEdge.HasVertex(lV0) && lEdge.HasVertex(lV1) )
				{
					lpFoundEdge = &lEdge;
					break;
				}
			}

			// Edge already in list
			if( lpFoundEdge )
			{
				// Add triangle
				lpFoundEdge->mT.AddTail( lTopTriIdx );
			}
			else
			{
				// Create new edge, add triangle & store edge
				TopEdge lNewEdge( lV0, lV1 );
				lNewEdge.mT.AddTail( lTopTriIdx );
				lEdgeList.AddTail( lNewEdge );
			}
		}
    }

#ifdef _DEBUG
    LzLogM("", "*** DEBUG HASH TABLE:                       Min key= "<<lMinKey)
    LzLogM("", "*** DEBUG HASH TABLE:                       Max key= "<<lMaxKey)
    LzLogM("", "*** DEBUG HASH TABLE: AMPLITUDE = Max key - Min key= "<<lMaxKey-lMinKey)
    lHashEdges.Log("HashTable for MeshTopology edges");
#endif

	// Flush hash table to vector
	lHashEdges.FlushToVector( mTopEdges );

	//--------------------------------------
	// Register edges
	//--------------------------------------

	for( vector<TopEdge>::const_iterator iTopEdge=mTopEdges.begin() ; iTopEdge!=mTopEdges.end() ; iTopEdge++ )
	{
		// Index of current edge
		int lTopEdgeIdx = (int)(iTopEdge - mTopEdges.begin());

		// To vertices
		mTopVers[ iTopEdge->mV[0] ].mE.AddTail( lTopEdgeIdx );
		mTopVers[ iTopEdge->mV[1] ].mE.AddTail( lTopEdgeIdx );

		// To triangles
		BrowseList( iT, iTopEdge->mT )
		{
			// Index and current triangle
			int lTopTriIdx = iTopEdge->mT.GetAt(iT);
			TopTri & lTopTri = mTopTris[ lTopTriIdx ];

			// TopTri::mE[0] is opposed to TopTri::mV[0], etc
			for( int v=0 ; v<3 ; v++ )
			{
				// Found triangle vertex opposite to edge
				if( !iTopEdge->HasVertex(lTopTri.mV[v]) )
				{
					// Set edge
					lTopTri.mE[v] = lTopEdgeIdx;
					break;
				}
			}
		}
	}

#ifdef _DEBUG
	//--------------------------------------
	// Check data coherency in *this
	//--------------------------------------

    _LzLogTry

	// Vertices
	int lUnusedVertices = 0;
	for( vector<TopVer>::const_iterator iTopVer=mTopVers.begin() ; iTopVer!=mTopVers.end() ; iTopVer++ )
	{
		// Check unused vertex: both lists must be empty
		if( iTopVer->mE.Count()==0 && iTopVer->mT.Count()==0 )
		{
			lUnusedVertices++;
			continue;
		}

		// Index of current vertex
		int lTopVerIdx = (int)(iTopVer - mTopVers.begin());

		// Vertex is used
		if( iTopVer->mE.Count()<2 || iTopVer->mT.Count()<1 )
            LzLogException("", "Edge inconsistency in vertex #"<<lTopVerIdx<<"!");

		// E
		BrowseList( iE, iTopVer->mE )
		{
			if( !mTopEdges[iTopVer->mE.GetAt(iE)].HasVertex(lTopVerIdx) )
                LzLogException("", "Edge inconsistency in vertex #"<<lTopVerIdx<<"!");
		}

		// T
		BrowseList( iT, iTopVer->mT )
		{
			// Throws exception if vertex not found
			mTopTris[iTopVer->mT.GetAt(iT)].EdgeWithVertex(lTopVerIdx);
		}
	}

	// Gentle reminder
	if( lUnusedVertices )
        LzLogMsg("", "Found "<<lUnusedVertices<<" unused vertices in mesh. Consider checking its consistency...")

	// Triangles
	for( vector<TopTri>::const_iterator iTopTri=mTopTris.begin() ; iTopTri!=mTopTris.end() ; iTopTri++ )
	{
		// Index of current triangle
		int lTopTriIdx = (int)(iTopTri - mTopTris.begin());

		// V
		if( iTopTri->mV[0]==-1 || iTopTri->mV[1]==-1 || iTopTri->mV[2]==-1
			|| iTopTri->mV[0]==iTopTri->mV[1] || iTopTri->mV[1]==iTopTri->mV[2] || iTopTri->mV[2]==iTopTri->mV[0] )
            LzLogException("", "Vertex inconsistency in triangle #"<<lTopTriIdx<<"!");

		// E
		if( iTopTri->mE[0]==-1 || iTopTri->mE[1]==-1 || iTopTri->mE[2]==-1
			|| iTopTri->mE[0]==iTopTri->mE[1] || iTopTri->mE[1]==iTopTri->mE[2] || iTopTri->mE[2]==iTopTri->mE[0] )
            LzLogException("", "Edge inconsistency in triangle #"<<lTopTriIdx<<"!");

		// Assumption
		for( int i=0 ; i<3 ; i++ )
		{
			if( mTopEdges[iTopTri->mE[i]].HasVertex(iTopTri->mV[i]) )
                LzLogException("", "Edge/vertex inconsistency in triangle #"<<lTopTriIdx<<"!");
		}
	}

	// Edges
	for( vector<TopEdge>::const_iterator iTopEdge=mTopEdges.begin() ; iTopEdge!=mTopEdges.end() ; iTopEdge++ )
	{
		// Index of current edge
		int lTopEdgeIdx = (int)(iTopEdge - mTopEdges.begin());

		// V
		if( iTopEdge->mV[0]==-1 || iTopEdge->mV[1]==-1 
			|| iTopEdge->mV[0]==iTopEdge->mV[1] )
            LzLogException("", "Vertex inconsistency in edge #"<<lTopEdgeIdx<<"!");

		// T
		if( iTopEdge->mT.Count() < 1 )
            LzLogException("", "Triangle inconsistency in edge #"<<lTopEdgeIdx<<"!");

		BrowseList( iT, iTopEdge->mT )
		{
			// Throws exception if edge not found
			mTopTris[iTopEdge->mT.GetAt(iT)].VertexOppositeToEdge(lTopEdgeIdx);
		}
	}

    _LzLogCatchAndThrow("Topology check failed!")
#endif

	//--------------------------------------
	// Find cracks and multi-edges
	//--------------------------------------

	for( vector<TopEdge>::const_iterator iTopEdge=mTopEdges.begin() ; iTopEdge!=mTopEdges.end() ; iTopEdge++ )
	{
		// Index of current edge
		int lTopEdgeIdx = (int)(iTopEdge - mTopEdges.begin());

		// Crack ?
		if( iTopEdge->mT.Count() == 1 )
			mCracks.push_back(lTopEdgeIdx);
		else
		// Multi-edge ?
		if( iTopEdge->mT.Count() > 2 )
			mMultiEdges.push_back(lTopEdgeIdx);
	}

	// That's all folks
    LzLogMsg("", "Results:")
    LzLogMsg("", "\t"<<mTopVers.size()<<" vertices,")
    LzLogMsg("", "\t"<<mTopTris.size()<<" triangles,")
    LzLogMsg("", "\t"<<mTopEdges.size()<<" edges,")
    LzLogMsg("", "\t"<<mCracks.size()<<" cracks,")
    LzLogMsg("", "\t"<<mMultiEdges.size()<<" multi-edges.")
}
catch( const std::exception & e )
{
    LzLogErr("", e.what());
    Free();
    LzLogException("", "Could not analyze mesh!");
}
}

//================================================================================
void MeshTopology::GetConnectedCracks( List<List<int>> & pTo ) const
{
    LzLogN("", "Sorting cracks by connexity.");
    LzLogM("", "Nb of cracks to be sorted: " << Cracks().size());

    // Free
    pTo.DelAll();

    // All the unsorted cracks
    const vector<int> & lAllCracksIdx = Cracks();
    const vector<TopEdge> & lAllCracks = TopEdges();

    // Sort cracks by connectivity
    for (vector<int>::const_iterator iCrack = lAllCracksIdx.begin(); iCrack != lAllCracksIdx.end(); iCrack++)
    {
        const TopEdge & lCrack = lAllCracks[*iCrack];

        // All the cracks connected to lCrack
        List<int> lConnectedCracks;

        // Browse already connected crack lists
        BrowseList(iListIdx, pTo)
        {
            const List<int> & lCrackList = pTo.GetAt(iListIdx);
            BrowseList(iCrackIdx, lCrackList)
            {
                // Check lCrack is connected to at least one crack of lCrackList
                const LzTriModel::TopEdge & lOtherCrack = lAllCracks[ lCrackList.GetAt(iCrackIdx) ];
                if (lCrack.IsConnected(lOtherCrack))
                {
                    // Retrieve all cracks in a separate list and remove it from mConnectedCracks
                    lConnectedCracks.Append(lCrackList);
                    pTo.DelAt(iListIdx);
                    break;
                }
            }
        }

        // lConnectedCracks contains all the cracks connected to lCrack
        lConnectedCracks.AddTail(*iCrack);
        pTo.AddTail(lConnectedCracks);
    }

    // Log
    LzLogM("", "Found " << pTo.Count() << " connected crack groups")
}

//*********** DEPRECATED
//================================================================================
//void MeshTopology::FindPatches( /*bool pMultiEdgesAreBoundaries,*/ List< List<size_t> > & pPatches ) const
//{
//	// Patch index for each triangle
//	vector<int> lTriToPatch;
//	lTriToPatch.assign( mTopTris.size(), -1 );
//
//	// Current patch id
//	int lCurrPatchId = 0;
//
//	// Triangles stack
//	List<size_t> lTriStack;
//
//	// Seed patches from all triangles until exhaustion
//	for( size_t t=0 ; t<mTopTris.size() ; t++ )
//	{
//		// Already in a patch?
//		if( lTriToPatch[t] > -1 )
//			continue;
//
//		// Start a new patch
//		lTriStack.AddTail( t );
//		lTriToPatch[ t ] = lCurrPatchId;
//
//		do
//		{
//			// Pop
//			size_t lTriIdx = lTriStack.GetHead();
//			lTriStack.DelHead();
//
//			// Push all 3 neighbors
//			for( int e=0 ; e<3 ; e++ )
//			{
//				const TopEdge & lTopE = mTopEdges[ mTopTris[lTriIdx].mE[e] ];
//
//				// Regular edge?
//				if( lTopE.mT.Count() == 2 )
//				{
//					// Get other triangle index
//					size_t lOtherTriIdx = lTopE.OtherTri( lTriIdx );
//
//					// Not yet in the patch?
//					if( lTriToPatch[lOtherTriIdx] == -1 )
//					{
//						lTriStack.AddTail( lOtherTriIdx );
//						lTriToPatch[ lOtherTriIdx ] = lCurrPatchId;
//					}
//				}
//			}
//		}
//		while( lTriStack.Count() );
//
//		// Next patch id
//		lCurrPatchId++;
//	}
//
//	// Log
//	TxLogNode("Found "+lCurrPatchId+" patches.");
//
//	// Remove previous
//	pPatches.DelAll();
//
//	// Convert to patch list
//	for( int p=0 ; p<lCurrPatchId ; p++ )
//	{
//		// Retrieve patch list
//		List<size_t> lPatchList;
//		for( size_t t=0 ; t<lTriToPatch.size() ; t++ )
//		{
//			if( lTriToPatch[t] == p )
//				lPatchList.AddTail( t );
//		}
//
//		// Store list
//		pPatches.AddTail( lPatchList );
//
//		// Log
//		TxLogM("Patch #"+p+": "+lPatchList.Count()+" triangles.");
//	}
//}
//*********** DEPRECATED
//
//================================================================================
//	bool MeshTopology::CheckPatchOrientation( const List<size_t> & pPatch ) const
//	{
//		// Empty patch?
//		if( !pPatch.Count() )
//			LzLogException("", "Cannot check orientation in an empty patch!");
//
//		// To avoid pushing same triangle twice
//		vector<bool> lPushed;
//		lPushed.assign( mTopTris.size(), false );
//
//		// Init
//		List<size_t> lStack;
//		lStack.AddTail( pPatch.GetHead() );
//		lPushed[pPatch.GetHead()] = true;
//
//		do
//		{
//			// Pop
//			size_t lTriIdx = lStack.GetHead();
//			lStack.DelHead();
//
//			// Check
//			if( lTriIdx >= mTopTris.size() )
//				LzLogException("", "Triangle index "+lTriIdx+" is out of range! Have only "+mTopTris.size()+".");
//
//			// Get current triangle
//			const TopTri & lCurrTopT = mTopTris[lTriIdx];
//
//			// Check all 3 neighbors
//			for( int e=0 ; e<3 ; e++ )
//			{
//				int lEdgeIdx = lCurrTopT.mE[e];
//				const TopEdge & lTopE = mTopEdges[lEdgeIdx];
//
////// Check if not multiedge
////if( lTopE.mT.Count() > 2 )
////	LzLogException("", "Found a multiedge in this patch! (index= "+lEdgeIdx+")");
//
//				// Regular edge
//				if( lTopE.mT.Count() == 2 )
//				{
//					int lOtherTriIdx = lTopE.OtherTri(lTriIdx);
//
//					// Already pushed?
//					if( lPushed[lOtherTriIdx] )
//						continue;
//
//					// Check orientation
//					if( !lCurrTopT.HasSameOrientationAs( mTopTris[lOtherTriIdx] ) )
//						return false;
//					else
//					{
//						lStack.AddTail( lOtherTriIdx );
//						lPushed[lOtherTriIdx] = true;
//					}
//				}
//			}
//		}
//		while( lStack.Count() );
//
//		// Everything went pretty well
//		return true;
//	}
#pragma endregion


#pragma region "Data"
//================================================================================
const TopVer & MeshTopology::TopVers( size_t pIdx ) const
{
	// Check
	if( pIdx >= mTopVers.size() )
        LzLogException("", "Cannot access TopVer "<<pIdx<<"! Size= "<<mTopVers.size()<<".")

	return mTopVers[ pIdx ];
}

//================================================================================
const TopTri & MeshTopology::TopTris( size_t pIdx ) const
{
	// Check
	if( pIdx >= mTopTris.size() )
        LzLogException("", "Cannot access TopTri "<<pIdx<<"! Size= "<<mTopTris.size()<<".")

	return mTopTris[ pIdx ];
}

//================================================================================
const TopEdge & MeshTopology::TopEdges( size_t pIdx ) const
{
	// Check
	if( pIdx >= mTopEdges.size() )
        LzLogException("", "Cannot access TopEdge "<<pIdx<<"! Size= "<<mTopEdges.size()<<".")

	return mTopEdges[ pIdx ];
}
#pragma endregion

#pragma endregion
}
