#include "Outliner.h"
#include "Mesh.h"
//#include "MeshTopology.h"
//
// Shewchuk's Triangle mesher
//
#ifdef COMPILE_SHEWCHUK
    #ifdef VOID
    #undef VOID
    #endif
    #define REAL double
    #define VOID int
    #include "TriangleMesher.h"
    #undef VOID
#endif
//
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>
#include <LzGeom/Line3D.h>


//*************** conflicting with std::numeric_limits<double>::max()
#undef max
//*************** conflicting with std::numeric_limits<double>::max()


namespace LzTriModel
{
using LzGeom::Line3D;


#pragma region "Construction & destruction"
//================================================================================
Outliner::Outliner() : mpMesh(nullptr), mpTopo(nullptr)
{
}

//================================================================================
Outliner::~Outliner()
{
	Free();
}

//================================================================================
void Outliner::Free()
{
	mpMesh = nullptr;
	mpTopo = nullptr;
    mTopo.Free();

	mOutVers.clear();
	mOutTrisChecked.clear();
	mOutEdges.clear();
}
#pragma endregion


#pragma region "Outlining"
//================================================================================
void Outliner::Set( const Mesh & pMesh, const MeshTopology & pTopo )
{
    // Free previous data
    Free();

    // Set containers
    mOutVers.resize( pTopo.TopVers().size() );
    mOutTrisChecked.resize( pTopo.TopTris().size() );
    mOutEdges.resize( pTopo.TopEdges().size() );

    // Link to mesh and topology
    mpMesh = &pMesh;
    mpTopo = &pTopo;
}

//================================================================================
void Outliner::Set( const Mesh & pMesh )
{
    // Free previous data
    Free();

    // Set my local topology
    mTopo.Set( pMesh );

    // Set containers
    mOutVers.resize( mTopo.TopVers().size() );
    mOutTrisChecked.resize( mTopo.TopTris().size() );
    mOutEdges.resize( mTopo.TopEdges().size() );

    // Link to mesh and topology
    mpMesh = &pMesh;
    mpTopo = &mTopo;
}

//================================================================================
#ifdef _DEBUG
// _DEBUG only statistics
static int sVertexCuts;
static int sEdgeCuts;
// _DEBUG only statistics
#endif

//================================================================================
void Outliner::Cut( const Plane3D & pCutPlane, List< List<Point3D> > & pOutlines, double pTolerance/*=1e-6*/ )
{
	//--------------------------------------
	// Init
	//--------------------------------------

	// Clear previous
	pOutlines.DelAll();

	// For all vertices compute distance to cut and status
    List<size_t> lStartingVers;
	for( vector<OutVer>::iterator iOutVer=mOutVers.begin() ; iOutVer!=mOutVers.end() ; iOutVer++ )
	{
		// Index of current vertex
        size_t lVerIdx = iOutVer - mOutVers.begin();

		// Skip unused vertices
		if( mpTopo->TopVers()[lVerIdx].mE.Count() == 0 )
			continue;

		// Compute distance and status
		iOutVer->mD = pCutPlane.SignedDistanceTo( mpMesh->mVertices[lVerIdx] );
		if( iOutVer->mD < -pTolerance )
			iOutVer->mStatus = -1;
		else
		if( iOutVer->mD > +pTolerance )
			iOutVer->mStatus = +1;
		else
		{
			iOutVer->mStatus = 0;

			// Store as potential starting point
			lStartingVers.AddTail( lVerIdx );
		}
	}

	// For all triangles reset checked flags
//	for( vector<bool>::iterator iOutTri=mOutTrisChecked.begin() ; iOutTri!=mOutTrisChecked.end() ; iOutTri++ )
//	{
//		// Reset checked flag
//		*iOutTri = false;
//	}
    std::fill(mOutTrisChecked.begin(), mOutTrisChecked.end(), false);

	// For all edges reset checked flags and compute strict cuts
	List<int> lStartingEdges;
	for( vector<OutEdge>::iterator iOutEdge=mOutEdges.begin() ; iOutEdge!=mOutEdges.end() ; iOutEdge++ ) 
	{
		// Reset checked flag
		iOutEdge->mAlreadyChecked = false;

		// Index of current edge
		int lEdgeIdx = (int)(iOutEdge - mOutEdges.begin());

		// Check for strict cut
		int lV0 = mpTopo->TopEdges()[lEdgeIdx].mV[0];
		int lV1 = mpTopo->TopEdges()[lEdgeIdx].mV[1];
		if( mOutVers[lV0].mStatus * mOutVers[lV1].mStatus == -1 )
		{
			// Mark edge as having a strict cut
			iOutEdge->mHasStrictCut = true;

			// Compute cut point
			double k = mOutVers[lV0].mD / (mOutVers[lV0].mD - mOutVers[lV1].mD);
			const Point3D & lPt0 = mpMesh->mVertices[lV0];
			const Point3D & lPt1 = mpMesh->mVertices[lV1];
			iOutEdge->mCutPoint = lPt0 + k*(lPt1 - lPt0);

			// Store as potential starting point
			lStartingEdges.AddTail( lEdgeIdx );
		}
		else
		{
			// Mark edge as not having a strict cut
			iOutEdge->mHasStrictCut = false;
		}
	}

	//--------------------------------------
    // Explore all starting entities
	//--------------------------------------
#ifdef _DEBUG
    LzLogMsg("", "Starting vertices: "<<lStartingVers.Count()<<", starting edges: "<<lStartingEdges.Count()<<".");

	// Reset statistics
	sVertexCuts = 0;
	sEdgeCuts = 0;
#endif
	// Edges
	BrowseList( iStartEdge, lStartingEdges )
	{
		// Current starting edge
		int lOutEdgeIdx = lStartingEdges.GetAt(iStartEdge);
		const OutEdge & lOutEdge = mOutEdges[ lOutEdgeIdx ];

		// Forward outline
		List<Point3D> lFwdOutline;
		CutFromEdge( lOutEdgeIdx, lFwdOutline );

		// No outline found starting from this edge: edge is exhausted ==> skip to next one
		if( !lFwdOutline.Count() )
			continue;

		// Complete outline
		lFwdOutline.AddHead( lOutEdge.mCutPoint );
#ifdef _DEBUG
		sEdgeCuts++;
#endif
		// Try a backwards outline
		List<Point3D> lBwdOutline;
		CutFromEdge( lOutEdgeIdx, lBwdOutline );

		// Merge outlines
		if( lBwdOutline.Count() )
		{
			BrowseList( iPt, lBwdOutline )
				lFwdOutline.AddHead( lBwdOutline.GetAt(iPt) );
		}

		// Store new outline (forward + backward)
		pOutlines.AddTail( lFwdOutline );

		// Starting edge may not be exhausted: stay on it
	}

	// Vertices
	BrowseList( iStartVer, lStartingVers )
	{
		// Current starting vertex
        const size_t lOutVerIdx = lStartingVers.GetAt(iStartVer);
//		const OutVer & lOutVer = mOutVers[ lOutVerIdx ];

		// Forward outline
		List<Point3D> lFwdOutline;
		CutFromVertex( lOutVerIdx, lFwdOutline );

		// No outline found starting from this vertex: vertex is exhausted ==> skip to next one
		if( !lFwdOutline.Count() )
			continue;

		// Complete outline
		lFwdOutline.AddHead( mpMesh->mVertices[lOutVerIdx] );
#ifdef _DEBUG
		sVertexCuts++;
#endif
		// Try a backwards outline
		List<Point3D> lBwdOutline;
		CutFromVertex( lOutVerIdx, lBwdOutline );

		// Merge outlines
		if( lBwdOutline.Count() )
		{
			BrowseList( iPt, lBwdOutline )
				lFwdOutline.AddHead( lBwdOutline.GetAt(iPt) );
		}

		// Store new outline (forward + backward)
		pOutlines.AddTail( lFwdOutline );

		// Starting vertex may not be exhausted: stay on it
	}
#ifdef _DEBUG
    LzLogMsg("", "Found "<<sVertexCuts<<" 'vertex cuts' and "<<sEdgeCuts<<" 'edge cuts'.");
#endif
}

//================================================================================
void Outliner::CutFromVertex( int pVerIdx, List<Point3D> & pOutline )
{
	// Search strict cut on (non checked-triangle)-connected opposite edges
    const List<size_t> & lT = mpTopo->TopVers()[pVerIdx].mT;
	BrowseList( iT, lT )
	{
		// Tri index
		int lTIdx = lT.GetAt(iT);

		// Skip checked triangles
		if( mOutTrisChecked[lTIdx] )
			continue;

		// Mark triangle as checked
		mOutTrisChecked[lTIdx] = true;

		// Index of opposite edge
		int lOppEIdx = mpTopo->TopTris()[lTIdx].EdgeWithoutVertex(pVerIdx);

		// Strict cut on edge?
		if( mOutEdges[lOppEIdx].mHasStrictCut )
		{
			// Add cut point
			pOutline.AddTail( mOutEdges[lOppEIdx].mCutPoint );
#ifdef _DEBUG
			sEdgeCuts++;
#endif
			// Continue from opposite edge
			CutFromEdge( lOppEIdx, pOutline );

			// Linear recursivity only
			return;
		}
	}

	// Search cut on (non checked-edge)-connected vertices
    const List<size_t> & lE = mpTopo->TopVers()[pVerIdx].mE;
	BrowseList( iE, lE )
	{
		// Edge index
		int lEIdx = lE.GetAt(iE);

		// Skip checked edges
		if( mOutEdges[lEIdx].mAlreadyChecked )
			continue;

		// Mark edge as checked
		mOutEdges[lEIdx].mAlreadyChecked = true;

		// Index of opposite vertex
		int lOppVIdx = mpTopo->TopEdges()[lEIdx].OtherVer(pVerIdx);

		// Vertex on cut?
		if( mOutVers[lOppVIdx].mStatus == 0 )
		{
			// Add cut point
			pOutline.AddTail( mpMesh->mVertices[lOppVIdx] );
#ifdef _DEBUG
			sVertexCuts++;
#endif
			// Continue from opposite vertex
			CutFromVertex( lOppVIdx, pOutline );

			// Linear recursivity only
			return;
		}
	}
}

//================================================================================
void Outliner::CutFromEdge( int pEdgeIdx, List<Point3D> & pOutline )
{
	// Search cut on (non checked-triangle)-connected opposite edges and vertices
	const List<int> & lT = mpTopo->TopEdges()[pEdgeIdx].mT;
	BrowseList( iT, lT )
	{
		// Tri index
		int lTIdx = lT.GetAt(iT);

		// Skip checked triangles
		if( mOutTrisChecked[lTIdx] )
			continue;

		// Mark triangle as checked
		mOutTrisChecked[lTIdx] = true;

		// Index of first opposite edge
		const TopTri & lTopTri = mpTopo->TopTris()[lTIdx];
		int lOppEIdx1 = lTopTri.OtherEdge(pEdgeIdx);

		// Strict cut on edge 1?
		if( mOutEdges[lOppEIdx1].mHasStrictCut )
		{
			// Add cut point
			pOutline.AddTail( mOutEdges[lOppEIdx1].mCutPoint );
#ifdef _DEBUG
			sEdgeCuts++;
#endif
			// Continue from opposite edge
			CutFromEdge( lOppEIdx1, pOutline );

			// Linear recursivity only
			return;
		}

		// Index of second opposite edge
		int lOppEIdx2 = lTopTri.OtherEdge(pEdgeIdx,lOppEIdx1);

		// Strict cut on edge 2?
		if( mOutEdges[lOppEIdx2].mHasStrictCut )
		{
			// Add cut point
			pOutline.AddTail( mOutEdges[lOppEIdx2].mCutPoint );
#ifdef _DEBUG
			sEdgeCuts++;
#endif
			// Continue from opposite edge
			CutFromEdge( lOppEIdx2, pOutline );

			// Linear recursivity only
			return;
		}

		// Index of opposite vertex
		int lOppVIdx = lTopTri.VertexOppositeToEdge(pEdgeIdx);

		// Vertex on cut?
		if( mOutVers[lOppVIdx].mStatus == 0 )
		{
			// Add cut point
			pOutline.AddTail( mpMesh->mVertices[lOppVIdx] );
#ifdef _DEBUG
			sVertexCuts++;
#endif
			// Continue from opposite vertex
			CutFromVertex( lOppVIdx, pOutline );

			// Linear recursivity only
			return;
		}
	}
}

//================================================================================
static void RecursiveMeshOutline( const List<int> & pOutIdx, Mesh & pMesh )
{
	if( pOutIdx.Count() == 3 )
	{
		//-------------------------------------------
		// Outline = triangle: finished splitting
		//-------------------------------------------

        void * iIdx = pOutIdx.HeadPos();
		int lV0 = pOutIdx.GetAtAndNext(iIdx);
		int lV1 = pOutIdx.GetAtAndNext(iIdx);
		int lV2 = pOutIdx.GetAt(iIdx);

		pMesh.mTriangles.push_back( Triangle(lV0,lV1,lV2,0,0,0) );
	}
	else
	{
		//-------------------------------------------
		// Need for split
		//-------------------------------------------

		// Best split
		void * lBestSplitA = nullptr;
		void * lBestSplitB = nullptr;

		// Read outline orientation
		const Vector3D & lCCW = pMesh.mNormals[0];

		// Maximal min outline vertex dist to admissible split
		double lMaxMinDistToSplit = -std::numeric_limits<double>::max();

		// Check all possible splits starting from A...
		void * lUpperBoundA = pOutIdx.PrevPos( pOutIdx.TailPos() );
		for( void * iSplitA=pOutIdx.HeadPos() ; iSplitA!=lUpperBoundA ; pOutIdx.Next(iSplitA) )
		{
			// Split vertex A
			const Point3D & lVerA = pMesh.mVertices[ pOutIdx.GetAt(iSplitA) ];

			// A prev and next vertices
			const Point3D & lVerA_prev = pMesh.mVertices[ iSplitA==pOutIdx.HeadPos() ? pOutIdx.GetTail() : pOutIdx.GetAt(pOutIdx.PrevPos(iSplitA)) ];
			const Point3D & lVerA_next = pMesh.mVertices[ pOutIdx.GetAt(pOutIdx.NextPos(iSplitA)) ];

			// Angle at A
			Vector3D lForwardA = lVerA_next - lVerA;
			double lMaxAngleA = lForwardA.PositiveAngleTo(lVerA_prev-lVerA,lCCW);

			// ... and arriving in B
			void * lUpperBoundB = iSplitA==pOutIdx.HeadPos() ? pOutIdx.TailPos() : nullptr ;
			for( void * iSplitB=pOutIdx.NextPos(pOutIdx.NextPos(iSplitA)) ; iSplitB!=lUpperBoundB ; pOutIdx.Next(iSplitB) )
			{
				// Split vertex B
				const Point3D & lVerB = pMesh.mVertices[ pOutIdx.GetAt(iSplitB) ];

				// Split vector
				Vector3D lSplitAB = lVerB - lVerA;

				// Found best split ever?
                if( LzMath::ToolBox::IsZero(lSplitAB.Norm(),1e-6) )
				{
					lBestSplitA = iSplitA;
					lBestSplitB = iSplitB;

					// Jump out of both nested loops
					goto FoundBestSplitEver;
				}

				//--------------------------------------------------
				// Check angle
				{
					// B prev and next vertices
					const Point3D & lVerB_prev = pMesh.mVertices[ pOutIdx.GetAt(pOutIdx.PrevPos(iSplitB)) ];
					const Point3D & lVerB_next = pMesh.mVertices[ iSplitB==pOutIdx.TailPos() ? pOutIdx.GetHead() : pOutIdx.GetAt(pOutIdx.NextPos(iSplitB)) ];

					// Check angle at A
					double lAngleA = lForwardA.PositiveAngleTo(lSplitAB,lCCW);
					if( lAngleA >= lMaxAngleA )
                        continue;

					// Angle at B
					Vector3D lForwardB = lVerB_next - lVerB;
					double lMaxAngleB = lForwardB.PositiveAngleTo(lVerB_prev-lVerB,lCCW);

					// Check angle at B
					double lAngleB = lForwardB.PositiveAngleTo(-lSplitAB,lCCW);
					if( lAngleB >= lMaxAngleB )
                        continue;
				}
				//--------------------------------------------------

				//--------------------------------------------------
				// Split line
				Line3D lSplitLine( lVerA, lVerB );

				// Check self-intersect
				{
					bool lFoundIntersection = false;

					for( void * iVer1=pOutIdx.HeadPos() ; iVer1 ; pOutIdx.Next(iVer1) )
					{
						// Skip current split
						if( iVer1==iSplitA || iVer1==iSplitB )
							continue;

						// Next vertex
						void * iVer2 = iVer1==pOutIdx.TailPos() ? pOutIdx.HeadPos() : pOutIdx.NextPos(iVer1);

						// Skip current split
						if( iVer2==iSplitA || iVer2==iSplitB )
							continue;
			
						// Check edge intersection with split
						const Point3D & lVer1 = pMesh.mVertices[ pOutIdx.GetAt(iVer1) ];
						const Point3D & lVer2 = pMesh.mVertices[ pOutIdx.GetAt(iVer2) ];

						// Compute lines intersection
						//Point3D lInter;
						Line3D lEdgeLine( lVer1, lVer2 );
						if( lSplitLine.IsParallelTo(lEdgeLine) )
						{
							// Lines are parallel

							// Compute distance between the two parallel lines
							if( lSplitLine.DistanceTo(lEdgeLine) < 1e-6 )
							{
                                // Exclude splits that "might be" overlapping an edge - /!\ other tests are required to actualy assess overlapping
								lFoundIntersection = true;
								break;
							}
						}
						else
						{
							// Lines are secant

							// Compute intersection
							Point3D lInter = lSplitLine.PseudoIntersection( lEdgeLine );

							// Lines are not parallel

							// Check intersection
							if( (lVerA-lInter)*(lVerB-lInter)<0 && (lVer1-lInter)*(lVer2-lInter)<0 )
							{
								lFoundIntersection = true;
								break;
							}
						}
#if 0 // OLD
					//	try
					//	{
					//		// Will throw exception if lines are parallel
					//		lInter = lSplitLine.PseudoIntersection( lEdgeLine );

					//		// Lines are not parallel

					//		// Check intersection
					//		if( (lVerA-lInter)*(lVerB-lInter)<0 && (lVer1-lInter)*(lVer2-lInter)<0 )
					//		{
					//			lFoundIntersection = true;
					//			break;
					//		}
					//	}
					//	catch( ... )
					//	{
					//		// Lines are parallel

					//		// Compute distance between the two parallel lines
					//		if( lSplitLine.DistanceTo(lEdgeLine) < 1e-6 )
					//		{
					//			// Exclude splits that "might be" overlapping an edge - /!\ other tests are required to actualy assess overlapping /!\ 
					//			lFoundIntersection = true;
					//			break;
					//		}
					//	}
#endif
					}

					// Found an intersection?
					if( lFoundIntersection )
						continue;
				}
				//--------------------------------------------------

				//--------------------------------------------------
				// Check min distance
				double lMinDistToSplit = +std::numeric_limits<double>::max();
				{
					// Check distance for all outline vertices other than those in split
					for( void * iVer=pOutIdx.HeadPos() ; iVer ; pOutIdx.Next(iVer) )
					{
						// Skip split vertices
						if( iVer==iSplitA || iVer==iSplitB )
							continue;

						// Compute distance from vertex to split
						const Point3D & lVer = pMesh.mVertices[ pOutIdx.GetAt(iVer) ];
						Point3D lProjVer = lSplitLine.Projection( lVer );

						// Localize projection and compute distance
						double lDistToSplit;
						if( (lVerA-lProjVer)*(lVerB-lProjVer) <= 0 )
							lDistToSplit = lVer.DistanceTo( lProjVer );
						else
						if( lSplitAB*(lProjVer-lVerA) > 0 )
							lDistToSplit = lVer.DistanceTo( lVerB );
						else
							lDistToSplit = lVer.DistanceTo( lVerA );

						// Found shorter dist?
						if( lMinDistToSplit > lDistToSplit )
							lMinDistToSplit = lDistToSplit;
					}
				}
				//--------------------------------------------------

				// Found better split?
				if( lMaxMinDistToSplit < lMinDistToSplit )
				{
					lBestSplitA = iSplitA;
					lBestSplitB = iSplitB;

					lMaxMinDistToSplit = lMinDistToSplit;
				}
			}
		}

//**************************************************** DEBUG
//**************************************************** DEBUG
if( !lBestSplitA || !lBestSplitB )
{
LzLogN("", "Outline that issued the error");
BrowseList( iIdx, pOutIdx )
    LzLogM("", ""+pMesh.mVertices[ pOutIdx.GetAt(iIdx)].ToString());
LzLogM("", ""+pMesh.mVertices[ pOutIdx.GetHead()].ToString());

LzLogException("", "Could not mesh this outline! ... because !lBestSplitA || !lBestSplitB");
}
//**************************************************** DEBUG
//**************************************************** DEBUG

FoundBestSplitEver:

		// Split
		List<int> lSplitOutIdx1;
		List<int> lSplitOutIdx2;
		List<int> * lpSplitOutIdx = &lSplitOutIdx1;

		for( void * iIdx=pOutIdx.HeadPos() ; iIdx ; pOutIdx.Next(iIdx) )
		{
			int lIdx = pOutIdx.GetAt(iIdx);

			// Hit a landmark?
			if( iIdx == lBestSplitA )
			{
				lSplitOutIdx1.AddTail( lIdx );
				lSplitOutIdx2.AddTail( lIdx );

				lpSplitOutIdx = &lSplitOutIdx2;
			}
			else
			if( iIdx == lBestSplitB )
			{
				lSplitOutIdx1.AddTail( lIdx );
				lSplitOutIdx2.AddTail( lIdx );

				lpSplitOutIdx = &lSplitOutIdx1;
			}
			else
				lpSplitOutIdx->AddTail( lIdx );
		}

		RecursiveMeshOutline( lSplitOutIdx1, pMesh );
		RecursiveMeshOutline( lSplitOutIdx2, pMesh );
	}
}

//================================================================================
void Outliner::MeshOutline( const List<Point3D> & pOutline,
                            const Vector3D & pCutNormal,
                            Mesh & pMesh )
{
	// Check that we have a closed outline with at least 4 points (triangle)
    if( pOutline.Count()<4 || !IsClosed(pOutline) )
        LzLogException("", "Cannot mesh an open outline!");

	// Set mesh normal
	pMesh.mNormals.clear();
	pMesh.mNormals.push_back( pCutNormal );
	pMesh.mNormals[0].Normalize();

	// Check outline orientation w.r.t. the specified normal
    bool lIsCCW = GetOrientation( pOutline, pCutNormal ) == Orientation::CCW;

	// Store vertices & build index list in CCW order w.r.t. the specified normal
	pMesh.mVertices.clear();
	List<int> lOutIdx;
	{
		int i = lIsCCW ? 0 : pOutline.Count()-2 ;

		for( void * iVer=pOutline.HeadPos() ; iVer!=pOutline.TailPos() ; pOutline.Next(iVer) )
		{
			pMesh.mVertices.push_back( pOutline.GetAt(iVer) );
			lOutIdx.AddTail( i );
			lIsCCW ? i++ : i-- ;
		}
	}

	// Compute mesh
	pMesh.mTriangles.clear();
	RecursiveMeshOutline( lOutIdx, pMesh );

    LzLogMsg("", "Generated mesh: "<<pMesh.mVertices.size()<<" vertices, "<<pMesh.mNormals.size()<<" normal(s), "<<pMesh.mTriangles.size()<<" triangle(s).");
}

#ifdef COMPILE_SHEWCHUK
//================================================================================
void Outliner::MeshOutline_Shewchuk( const List<Point3D> & pOutline, const Vector3D & pCutNormal, double pMaxTriArea, Mesh & pMesh )
{
	// Check that we have a closed outline with at least 4 points (triangle)
    if( pOutline.Count()<4 || !LzMath::ToolBox::IsZero( pOutline.GetHead().DistanceTo(pOutline.GetTail()) ) )
        LzLogException("", "Cannot mesh an open outline!");

	// Create local referential
	Point3D lEigO;
	Vector3D lEigVs[3];
    LzMath::ToolBox::PCA( pOutline, lEigO, lEigVs );

	// Orient PCA referential according to pCutNormal
	if( (lEigVs[0]^lEigVs[1]) * pCutNormal < 0 )
		lEigVs[1] *= -1;

	// Transforms
	const RigidTr3D lAlias_PCA_to_W( lEigVs[0], lEigVs[1], lEigO );
	const RigidTr3D lAlias_W_to_PCA = lAlias_PCA_to_W.Inverse();

	// Create data structures
	struct triangulateio in, mid, out;
	Reset_triangulateio( in );
	Reset_triangulateio( mid );
	Reset_triangulateio( out );

	// Fill points
	in.numberofpoints = pOutline.Count() - 1;
	in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
	{
		int i = 0;
		BrowseList( iP, pOutline )
		{
			// Convert point to PCA coordinates
			Point3D lPt = pOutline.GetAt( iP );
			lPt *= lAlias_W_to_PCA;

			in.pointlist[i++] = lPt.X();
			in.pointlist[i++] = lPt.Y();

			// Skip last, repeated, point in pOutline
			if( i == in.numberofpoints * 2 )
				break;
		}
	}

	// Fill segments
	in.numberofsegments = in.numberofpoints;
	in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
	for( int s=0 ; s<in.numberofsegments ; s++ )
	{
		in.segmentlist[2*s + 0] = s;
		in.segmentlist[2*s + 1] = (s + 1) % in.numberofpoints;
	}

	/* Triangulate the points.  Switches are chosen to read and write a  */
	/*   PSLG (p), preserve the convex hull (c), number everything from  */
	/*   zero (z), assign a regional attribute to each element (A), and  */
	/*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
	/*   neighbor list (n).                                              */


//*** DEBUG
//{
//    LzLogN("", "Points")
//    for( int p=0 ; p<in.numberofpoints ; p++ )
//        LzLogM("", "Point "<<p<<"= "<<in.pointlist[2*p]<<", "<<in.pointlist[2*p+1]);
//}
//{
//    LzLogN("", "Segments")
//    for( int s=0 ; s<in.numberofsegments ; s++ )
//        LzLogM("", "Segment "<<s<<"= "<<in.segmentlist[2*s]<<", "<<in.segmentlist[2*s+1]);
//}
//*** DEBUG


	//triangulate("pczAevn", &in, &mid, &vorout);
	//triangulate("pcz", &in, &mid, nullptr);
    triangulate((char*)"pzY", &in, &mid, nullptr);
	//*** Y prohibits the insertion of Steiner points on external boundaries


//*** DEBUG
//    LzLogMsg("", "triangulate((char*)'pzY', &in, &mid, nullptr) : OK")
//*** DEBUG


	// Needed only if -r and -a switches used
	mid.trianglearealist = (REAL *) malloc(mid.numberoftriangles * sizeof(REAL));
	for( int t=0 ; t<mid.numberoftriangles ; t++ )
		mid.trianglearealist[t] = pMaxTriArea;

    triangulate((char*)"prazY", &mid, &out, nullptr);
	//triangulate("prazBP", &mid, &out, nullptr);
	//triangulate("prazP", &mid, &out, nullptr); *** crash
	//*** Y prohibits the insertion of Steiner points on external boundaries


//*** DEBUG
//    LzLogMsg("", "triangulate((char*)'prazY', &mid, &out, nullptr) : OK")
//*** DEBUG


    // Generate mesh
	pMesh.Free();
	pMesh.mVertices.resize( out.numberofpoints );
	for( int v=0 ; v<out.numberofpoints ; v++ )
		pMesh.mVertices[v] = lAlias_PCA_to_W * Point3D( out.pointlist[2*v], out.pointlist[2*v + 1], 0.0 );
	//
	pMesh.mNormals = { pCutNormal };
	//
	pMesh.mTriangles.resize( out.numberoftriangles );
	for( int t=0 ; t<out.numberoftriangles ; t++ )
        pMesh.mTriangles[t] = LzTriModel::Triangle(out.trianglelist[3*t], out.trianglelist[3*t + 1], out.trianglelist[3*t + 2], 0, 0, 0);

	// Clean-up mess
	Delete_triangulateio( in );
	Delete_triangulateio( mid );
	Delete_triangulateio( out );
}
#endif
#pragma endregion


#pragma region "Tools"
//================================================================================
bool Outliner::IsClosed( const List<Point3D> & pOutline )
{
	// Check
	if( pOutline.Count() < 2 )
        LzLogException("", "Cannot determine closedness of the outline! Too few points (count= "<<pOutline.Count()<<").");
#if 0
    // Test that identical points undergoing identical transformations retain IDENTICAL
    // internal representation, making the operator == suitable for comparing doubles
    // IN THIS SPECIFIC CONTEXT
    //
    std::uniform_real_distribution<double> lDist(-1000, +1000 );
    Point3D lP( lDist(LzServices::RandomEngine()), lDist(LzServices::RandomEngine()), lDist(LzServices::RandomEngine()) );

    Point3D lQ = lP;

    LzLogM("", "P= "<<lP.ToString())
    LzLogM("", "Q= "<<lQ.ToString())

    for( int i=0 ; i<100 ; i++ )
    {
        RigidTr3D lTr( lDist(LzServices::RandomEngine()), lDist(LzServices::RandomEngine()), lDist(LzServices::RandomEngine()), // Vector3D() );
                       Vector3D( lDist(LzServices::RandomEngine()), lDist(LzServices::RandomEngine()), lDist(LzServices::RandomEngine()) ) );
        lP *= lTr;
        lQ *= lTr;
    }

    LzLogM("", "P= "<<lP.ToString())
    LzLogM("", "Q= "<<lQ.ToString())

    if( lP == lQ )
        LzLogInfoMessage("OK", "")
    else
        LzLogInfoMessage("*** KO", "")


//    //**** WARNING ****
//    static bool sShowOnce = true;
//    if( sShowOnce )
//    {
//        sShowOnce = false;
//        for( int i=0 ; i<10 ; i++ )
//            LzLogMsg("", "Outliner::IsClosed: Using == comparison between tail and head points in outline...")
//    }
#endif

    // Outline is CLOSED iff head and tail points are IDENTICAL
    //
    // Which means that they were generated by the Outliner class and
    // underwent the same transformations, see tests above
    return pOutline.GetHead() == pOutline.GetTail();    
}

//================================================================================
Outliner::Orientation Outliner::GetOrientation( const List<Point3D> & pOutline, const Vector3D & pRef )
{
	// Check
	if( !IsClosed(pOutline) )
        LzLogException("", "Cannot compute orientation of an open outline!")

	// Compute orientation
    double lAngle = 0;
    void * iPos = pOutline.HeadPos();
    do
    {
        const Point3D & lA = pOutline.GetAtAndNext(iPos);
        const Point3D & lB = pOutline.GetAt(iPos);

        void * lNextPos;
        if( iPos == pOutline.TailPos() )
            lNextPos = pOutline.NextPos( pOutline.HeadPos() );
        else
            lNextPos = pOutline.NextPos( iPos );

        const Point3D & lC = pOutline.GetAt(lNextPos);

        Vector3D U = lB - lA;
        Vector3D V = lC - lB;

        // Ignore overlappling points
        double lDeltaAngle;
        if( U.SignedAngleTo_NoExc( V, pRef, lDeltaAngle ) )
            lAngle += lDeltaAngle;
    }
    while( iPos != pOutline.TailPos() );


//LzLogMsg("", "Count= "<<pOutline.Count()<<", Angle= "<<lAngle)
//while(true);


    // Angle = +2Pi if outline is ordered CCW w.r.t. pCutNormal, -2Pi otherwise
//return lAngle > 0 ? Orientation::CCW : Orientation::CW ;
    if( std::abs(lAngle - LzMath::DPI) < 1e-6 ) //*** old version lAngle > +LzMath::PI )
    {
        // Actually, if no self intersection then = +2Pi
        return Orientation::CCW;
    }
    else
    if( std::abs(lAngle + LzMath::DPI) < 1e-6 ) //*** old version lAngle < -LzMath::PI )
    {
        // Actually, if no self intersection then = -2Pi
        return Orientation::CW;
    }
    else
    {
        // If self intersecting, then = 0
        return Orientation::Undefined;
    }
}

//================================================================================
string Outliner::OrientationToString( Orientation pOr )
{
    switch( pOr )
    {
    case Orientation::CCW: return "CCW";
    case Orientation::CW: return "CW";
    case Orientation::Undefined: return "Undefined";
    default:
        LzLogException("", "Unexpected orientation!")
    }
}

//================================================================================
double Outliner::OutlineLength( const List<Point3D> & pOutline )
{
	double lLen = 0;
	BrowseList( iPt, pOutline )
	{
		void * iNextPt = pOutline.NextPos( iPt );
		if( iNextPt )
			lLen += pOutline.GetAt(iPt).DistanceTo( pOutline.GetAt(iNextPt) );
		else
			break;
	}

	return lLen;
}

//================================================================================
double Outliner::SignedOutlineArea( const List<Point3D> & pOutline, const Vector3D * ppNormal/*=nullptr*/ )
{
	// Check
	if( !IsClosed(pOutline) )
        LzLogException("", "Cannot compute area of an open outline!");

	// Compute outline orientation
	Vector3D lMyNormal;
	const Vector3D * lpNor;
	if( ppNormal )
	{
		// Use user provided normal
		lpNor = ppNormal;
	}
	else
	{
		// Compute normal using a least squares plane
		Plane3D lLSQ( pOutline );
		lMyNormal = lLSQ.Normal();
		lpNor = &lMyNormal;
	}

	// Find outline center (for better accuracy)
    Point3D lOutCen = LzMath::ToolBox::Centroid( pOutline );

	// Compute area
    double lSgnOutArea = 0.0;
	BrowseList( iP, pOutline )
	{
		// Get segment endpoints
		const Point3D & lA = pOutline.GetAt( iP );
		void * iNextP = pOutline.NextPos( iP );
		const Point3D & lB = iNextP ? pOutline.GetAt( iNextP ) : pOutline.GetHead() ;

		// Compute area vector
		const Vector3D lAreaVec = (lA - lOutCen)^(lB - lOutCen);
		const double lSgnArea = 0.5 * lAreaVec * (*lpNor);

		// Update
        lSgnOutArea += lSgnArea;
	}

//	// Take absolute value before returning
//	return fabs( lSgnOutArea );
//*** NO, return signed Area, which makes it possible to check for outline orientation
    return lSgnOutArea;
}

//================================================================================
Point3D Outliner::OutlineWeightedCentroid( const List<Point3D> & pOutline )
{
	// Check
	if( pOutline.Count() == 0 )
        LzLogException("", "Empty outline! Undefined centroid.");

	// Check
	if( pOutline.Count() == 1 )
		return pOutline.GetHead();

	// Outline has >= 2 points
	Vector<Point3D> lPoints( pOutline.Count() - 1 );
	Vector<double> lLengths( pOutline.Count() - 1 );
	int i = 0;
	BrowseList( iPt, pOutline )
	{
		void * iNextPt = pOutline.NextPos( iPt );
		if( iNextPt )
		{
			const Point3D & lA = pOutline.GetAt(iPt);
			const Point3D & lB = pOutline.GetAt(iNextPt);

			lPoints[i]  = lA.MidPointTo( lB );
			lLengths[i] = lA.DistanceTo( lB );

			i++;
		}
		else
			break;
	}

    return LzMath::ToolBox::WeightedCentroid( lPoints, lLengths );
}

//================================================================================
void Outliner::LogOutlinesInfo( const List< List<Point3D> > & pOutlines, const std::string & pTag/*="Outlines"*/ )
{
    LzLogMsg("", pTag<<": "<<pOutlines.Count()<<" outline(s) in set.");

	int i = 0;
	BrowseList( iOut, pOutlines )
	{
		const List<Point3D> & lOut = pOutlines.GetAt(iOut);
		if( IsClosed(lOut) )
            LzLogMsg("", pTag<<"["<<i<<"]: "<<(lOut.Count()-1)<<"+1 points in CLOSED outline.")
		else
            LzLogMsg("", pTag+"["<<i<<"]: "<<lOut.Count()<<" points in OPEN outline.")

		i++;
	}
}

//================================================================================
void Outliner::OutlineIntersections( const List<Point3D> & pOutline, const Plane3D & pCutPlane, List<Point3D> & pInters, double pTolerance/*=1e-6*/ )
{
	// Check
	if( pOutline.Count() < 2 )
        LzLogException("", "Invalid outline!");

	// Clear previous
	pInters.DelAll();

	// Each segment is considered as [A; B[
	// - Open outline ==> need to consider the last point of the outline: [A; B[ _and_ B
	// - Closed outline ==> do not consider the last point as it is also the first point
	//
	const bool lTestLastPoint = !Outliner::IsClosed( pOutline );

	// Compute intersection for the current outline
	void * iPos = pOutline.HeadPos();
	void * iNextPos = pOutline.NextPos( iPos );
	do
	{
		// Segment [A; B[
		const Point3D & lA = pOutline.GetAt( iPos );
		const Point3D & lB = pOutline.GetAt( iNextPos );

		// Distances
		double lDstA = pCutPlane.SignedDistanceTo( lA );
		double lDstB = pCutPlane.SignedDistanceTo( lB );

		// Status A
		int lStatA;
		if( lDstA < -pTolerance )
			lStatA = -1;
		else
		if( lDstA > +pTolerance )
			lStatA = +1;
		else
			lStatA = 0;

		// A on cut ?
		if( lStatA == 0 )
		{
			// A is on cut
			pInters.AddTail( lA );
		}
		else
		{
			// Status B
			int lStatB;
			if( lDstB < -pTolerance )
				lStatB = -1;
			else
			if( lDstB > +pTolerance )
				lStatB = +1;
			else
				lStatB = 0;

			// Only consider case where B is not on cut and have a clean cut
			if( lStatA*lStatB < 0 )
			{
				// Compute weights
				double lAbsA = fabs( lDstA );
				double lAbsB = fabs( lDstB );
				double lSum = lAbsA + lAbsB;

				// Interpolate clean cut
				Point3D lCleanCut;
				for( int i=0 ; i<3 ; i++ )
					lCleanCut.mV[i] = (lAbsB*lA.mV[i] + lAbsA*lB.mV[i]) / lSum;

				// Stash
				pInters.AddTail( lCleanCut );
			}
		}

		// Next segment
		iPos = iNextPos;
		iNextPos = pOutline.NextPos( iPos );
	}
	while( iNextPos );

	// Last point?
	if( lTestLastPoint )
	{
		const Point3D & lB = pOutline.GetTail();

		if( fabs(pCutPlane.SignedDistanceTo( lB )) <= pTolerance )
			pInters.AddTail( lB );
	}
}

//================================================================================
void Outliner::OutlineIntersections( const List< List<Point3D> > & pOutlines, const Plane3D & pCutPlane, List< List<Point3D> > & pInters, double pTolerance/*=1e-6*/  )
{
	// Clear previous
	pInters.DelAll();

	// Compute intersection for each outline
	BrowseList( iOuts, pOutlines )
	{
		// Compute intersections for this outline
		List<Point3D> lInters;
		OutlineIntersections( pOutlines.GetAt(iOuts), pCutPlane, lInters, pTolerance );

		// Stash
		pInters.AddTail( lInters );
	}
}

//================================================================================
List<Point3D> Outliner::ConvexifyOutline( /*const*/ List<Point3D> /*&*/ pClosedOut, const Vector3D * ppNormal/*=nullptr*/ )
{
#if 1
    // Check
    if( !IsClosed( pClosedOut ) )
        LzLogException("", "Cannot convexify outline! Outline ("<<pClosedOut.Count()<<" point(s)) is not closed.")

    // Find a normal vector to the polygon
    const Vector3D lNormal = ppNormal ? *ppNormal : Plane3D( pClosedOut ).Normal() ;

    /* Proof that MaxPt is a singular point on the convex hull
     * -------------------------------------------------------
     *
     * Find MaxPt, the outline point furthest away from the coplanar reference point
     * We know that this point will lie on the convex hull since all other points are
     * withing a circle centered on RefPt with a radius from RefPt to MaxPt, and as a
     * consequence, no other outline point lie past the tangent to the circle at MaxPt.
     *
     * Which means that MaxPt is an extremum of the outline along the radius direction.
     *
     * Any point of the outline can be obtained as a linear combination of 2 points from
     * the convex hull. Now since no point is further than the tangent, these two points
     * will have to be on the tangent itself to be able to produce MaxPt. However there
     * are no other points on the tangent except MaxPt (if there were, these points would
     * be dubbed MaxPt instead of MaxPt as their distance to RefPt would be greater).
     *
     * So bottom line, not only is MaxPt on the convex hull, he is there by himself, as
     * a singularity, and the 2 outgoing edges point towards the inside of the circle.
     *
     * In other words MaxPt cannot either be a redundant point lying on a segment of
     * the convex hull. It is a good starting and ending point for our convexification
     * process.
     *
     */
    void * iMaxPos = nullptr;
    {
        // Compute polygon plane (without computing LSQ plane)
        const Plane3D lOutPlane( pClosedOut.GetHead(), lNormal );

        // Pick a reference point
        const Point3D lRefPt = lOutPlane.Projection( Point3D{0,0,0} );

        // Find max points
        double lMaxDist = 0.0;
        BrowseList( iPt, pClosedOut )
        {
            const double lDist = lRefPt.DistanceTo( pClosedOut.GetAt(iPt) );
            if( lMaxDist < lDist )
            {
                lMaxDist = lDist;
                iMaxPos = iPt;
            }
        }
    }



//*************** si pas de modification de l'outline, alors besoin de circular NEXT --> a rajouter dans LzOutlinerTools
//*************** si pas de modification de l'outline, alors besoin de circular NEXT --> a rajouter dans LzOutlinerTools
//*************** si pas de modification de l'outline, alors besoin de circular NEXT --> a rajouter dans LzOutlinerTools
//*************** si pas de modification de l'outline, alors besoin de circular NEXT --> a rajouter dans LzOutlinerTools
//*************** si pas de modification de l'outline, alors besoin de circular NEXT --> a rajouter dans LzOutlinerTools


    //******************************* Rotate outline (deal with the fact that the last point is repeated)
    {
        if( iMaxPos != pClosedOut.HeadPos() )
        {
//LzLogM("", "------ Rotating")

            // Remove duplicate head
            pClosedOut.DelHead();

            // Rotate
            pClosedOut.RotateToHead( iMaxPos );

            // Close outline
            pClosedOut.AddTail( pClosedOut.GetHead() );

            // Update max pos
            iMaxPos = pClosedOut.HeadPos();
        }
        else
        {
            // If iMaxPos == HeadPos, then nothing to do
        }
    }




    //******************************* Reverse outline if not CCW
    {
        // Determine orientation
        const Orientation lOrient = GetOrientation( pClosedOut, lNormal );

        if( lOrient != Orientation::CCW )
        {
//LzLogM("", "***** Reversing")
            pClosedOut.Reverse();

            // Update max pos
            iMaxPos = pClosedOut.HeadPos();
        }
    }


//LzLogM("", "***** MaxPt= "<<pClosedOut.GetHead().ToString())
//LzLogM("", "***** lNormal= "<<lNormal.ToString())


    // Initialize convex outline
    List<Point3D> lRes;

    // Initialize convex hull point pos
    void * iHullPos = iMaxPos;

    // Initialize rotation reference vector
    Vector3D lRotRef = lNormal ^ (pClosedOut.GetAt(iMaxPos) - Point3D{0,0,0}); // Well defined vector

//LzLogM("", "lRotRef= "<<lRotRef.ToString())

    // Iterate
    while( true )
    {
//********* Garder references sur points au lieu d'acceder via GetAt
//********* Garder references sur points au lieu d'acceder via GetAt
//********* Garder references sur points au lieu d'acceder via GetAt


        // Add hull point
        lRes.AddTail( pClosedOut.GetAt(iHullPos) );

        // Get Next pos (NO NEED FOR 'CIRCULAR NEXT' HERE, since the list has been rotated above)
        void * lNextPos = pClosedOut.NextPos(iHullPos);

        // Find next hull pos if any
        iHullPos = nullptr;
        double lMinAngle_rad = 2.0*LzMath::PI + 0.1;
        //
        for( void * iPos=lNextPos ; iPos ; pClosedOut.Next(iPos) )
        {
            // Get outline point
            const Point3D & lPt = pClosedOut.GetAt( iPos );

            try
            {
                // Compute angle
                const double lAngle_rad = lRotRef.ShortPositiveAngleTo( lPt - lRes.GetTail() );

//LzLogM("", lPt.ToString()<<" --> "<<LzMath::RAD_2_DEG*lAngle_rad)

                // New min?
                if( lMinAngle_rad > lAngle_rad )
                {
                    lMinAngle_rad = lAngle_rad;
                    iHullPos = iPos;

//LzLogM("","New min")
                }
            }
            catch( ... )
            {
                // NOP: ignore outline aberrations
            }
        }

        // Break if no next hull pos
        if( !iHullPos )
            break;

        // Update rotation reference, if found a hull pos
        lRotRef = pClosedOut.GetAt(iHullPos) - lRes.GetTail();

//LzLogM("", "---------------------------- New lRotRef= "<<lRotRef.ToString())
    }

    // Close outline if needed
    // (either because the outline collapsed into a single point black hole,
    // or because for some reason only the final point is missing)
    if( lRes.Count()==1 || !Outliner::IsClosed(lRes) )
    {
        LzLogM("", "Closing outline")
        lRes.AddTail( lRes.GetHead() );
    }

    return lRes;
#else
    LzLogException("", "**** CODE INCORRECT")


    // Check
    if( !IsClosed( pClosedOut ) )
        LzLogException("", "Cannot convexify outline! Outline ("<<pClosedOut.Count()<<" point(s)) is not closed.")

// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)
// **** need at least 2 points (i.e. Count >= 3)



    // Find a normal vector to the polygon
    const Vector3D lNormal = ppNormal ? *ppNormal : Plane3D( pClosedOut ).Normal() ;

    /* Proof that MaxPt is a singular point on the convex hull
     * -------------------------------------------------------
     *
     * Find MaxPt, the outline point furthest away from the coplanar reference point
     * We know that this point will lie on the convex hull since all other points are
     * withing a circle centered on RefPt with a radius from RefPt to MaxPt, and as a
     * consequence, no other outline point lie past the tangent to the circle at MaxPt.
     *
     * Which means that MaxPt is an extremum of the outline along the radius direction.
     *
     * Any point of the outline can be obtained as a linear combination of 2 points from
     * the convex hull. Now since no point is further than the tangent, these two points
     * will have to be on the tangent itself to be able to produce MaxPt. However there
     * are no other points on the tangent except MaxPt (if there were, these points would
     * be dubbed MaxPt instead of MaxPt as their distance to RefPt would be greater).
     *
     * So bottom line, not only is MaxPt on the convex hull, he is there by himself, as
     * a singularity, and the 2 outgoing edges point towards the inside of the circle.
     *
     * In other words MaxPt cannot either be a redundant point lying on a segment of
     * the convex hull. It is a good starting and ending point for our convexification
     * process.
     *
     */
    void * iMaxPos = nullptr;
    {
        // Compute polygon plane (without computing LSQ plane)
        const Plane3D lOutPlane( pClosedOut.GetHead(), lNormal );

        // Pick a reference point
        const Point3D lRefPt = lOutPlane.Projection( Point3D{0, 0, 0} );

        // Find max points
        double lMaxDist = 0.0;
        BrowseList( iPt, pClosedOut )
        {
            const double lDist = lRefPt.DistanceTo( pClosedOut.GetAt(iPt) );
            if( lMaxDist < lDist )
            {
                lMaxDist = lDist;
                iMaxPos = iPt;
            }
        }
    }

    // To make sure that the stop criterion works properly, otherwise
    // P0 will never reach HeadPos, while skipping using Next below
    if( iMaxPos == pClosedOut.HeadPos() )
        iMaxPos = pClosedOut.TailPos();




    // Check (pretty useless but ok)
    if( iMaxPos == nullptr )
        LzLogException("", "Could not find MaxPt, for some obscure reason!")

    // Determine orientation
    const Orientation lOrient = GetOrientation( pClosedOut, lNormal );




    // Adjust orientation vector according to Orientation
    const Vector3D lPosOrient = lOrient==Orientation::CCW ? +lNormal : -lNormal ;




    // Next position circular finder
    auto Next = [&pClosedOut]( void * pPos ) -> void *
        {
            // Find regular next pos
            void * lNextPos = pClosedOut.NextPos(pPos);

            // Went overboard?
            if( lNextPos == nullptr )
                lNextPos = pClosedOut.NextPos( pClosedOut.HeadPos() );

            return lNextPos;
        };


LzLogN("", "ConvexifyOutline")

    // Use iMaxPt to remember start pos, and init cursors
    List<Point3D> lResOut;
    void * iP0 = iMaxPos;
    do
    {
LzLogN("", "P0= "<<pClosedOut.GetAt(iP0).ToString())

        // Store current starting point
        lResOut.AddTail( pClosedOut.GetAt(iP0) );

        // P1: look for next point on convex hull; P1 runs in ]P0, iMaxPos[
        void * iP1;
        for( iP1=Next(iP0) ; iP1!=iMaxPos ; iP1=Next(iP1) )
        {
LzLogN("", "Considering P1= "<<pClosedOut.GetAt(iP1).ToString())


            // Compute reference vector P0P1
            const Vector3D P0P1 = pClosedOut.GetAt(iP1) - pClosedOut.GetAt(iP0);

            // Check norm for angle computations
            if( P0P1.Norm() < 1e-6 )
            {
                // Skip this candidate since if we chose it, it would be
                // redundant with the current point P0
                continue;
            }

            // Check if iP1 is a convex hull point

            // No need to check the first convex vertex since everybody is BENEATH the first edge and P1 is BENEATH the edge P-1P0
            // The angle between MaxPos-P0 and MaxPos-P1 can be deduced from the consecutive angles over the convex contour
            // Add proper explanation here !!!
            // Add proper explanation here !!!
            // Add proper explanation here !!!
            for( void * iP2=Next(iP1) ; iP2!=iMaxPos ; iP2=Next(iP2) )
            {
LzLogN("", "Test point P2= "<<pClosedOut.GetAt(iP2).ToString())


                // Compute test vector P1P2
                const Vector3D P1P2 = pClosedOut.GetAt(iP2) - pClosedOut.GetAt(iP1);

                // Check norm for angle computations
                if( P1P2.Norm() < 1e-6 )
                {
                    // Ignore this test point P2 as it won't have an impact on the
                    // rejection of P1
                    continue;
                }

                // Compute angle from P0P1 to P1P2 to see if P0P1 is an edge on the convex hull
                const double lSgnRad = P0P1.SignedAngleTo( P1P2, lPosOrient );
LzLogM("", "Alpha= "<<LzMath::RAD_2_DEG * lSgnRad)

                // Check
                if( lSgnRad <= 0.0 )
                {
LzLogM("", "P1 rejected, negative angle")
                    // Nope: P1 is no good since it creates a concavity or a flat edge
                    goto ConsiderNextP1;
                }
LzLogM("", "Okay for P2= "<<iP2)
            }


            // P1 passed the test
LzLogM("", "P1 passed the test")
            goto FoundNewHullPoint;

            // P1 failed the test
        ConsiderNextP1:
//LzLogM("", "P1 failed the test for P2= "<<pClosedOut.GetAt(iP2).ToString())
                ;
        }

        // Could not find suitable candidate P1 in ]P0, iMaxPos[
        // The only choice left is iMaxPos
        break;

    FoundNewHullPoint:
LzLogM("", "Moving P0= "<<pClosedOut.GetAt(iP0).ToString()<<" to "<<pClosedOut.GetAt(iP1).ToString())
        // Move pointer
        iP0 = iP1;
    }
    while( iP0 != iMaxPos );

    // Close the outline and return
    lResOut.AddTail( pClosedOut.GetAt(iMaxPos) );


//*** DEBUG
{
    LzLogN("", "CONVEX")
    BrowseList( iP, lResOut )
        LzLogM("", lResOut.GetAt(iP).ToString())
}
//*** DEBUG



    return lResOut;
#endif
}

//================================================================================
static double DistanceToSegment( const Point3D & pPt, const Point3D & pA, const Point3D & pB, Point3D & pProj )
{
	double lMinDist;

	// A
	lMinDist = pPt.DistanceTo( pA );
	pProj = pA;

	// B
	double lDistB = pPt.DistanceTo( pB );
	if( lMinDist > lDistB )
	{
		lMinDist = lDistB;
		pProj = pB;
	}

	// [A, B]
	Line3D lLine( pA, pB );
	Point3D lProj = lLine.Projection( pPt );

	// If point projects out of segment return current min
	if( (pA-lProj) * (pB-lProj) > 0 )
		return lMinDist;

	// Orthogonal projection gives the minimum distance to segment
	pProj = lProj;
	return pPt.DistanceTo( lProj );
}

//================================================================================
double Outliner::SgnDistToClosedOut( const List<Point3D> & pOutline, const Point3D & pPt, Point3D & pProjP, const Vector3D * ppNormal/*=nullptr*/)
{
    // Check
    if( !IsClosed(pOutline) )
        LzLogException("", "Cannot compute signed distance to closed outline! Outline is NOT closed.")

    // Find a normal vector to the polygon
    const Vector3D & lNorm = ppNormal ? *ppNormal : Plane3D(pOutline).Normal() ;

    // Compute projection and sum of angles to determine in/out sign
	double lMinPosDist = std::numeric_limits<double>::max();
    double lAngleSum = 0;

	// Check all segments in closed outline: last point repeats first so no need to go around the list
	BrowseList( iPt, pOutline )
	{
		// Last segment?
		void * iNextPt = pOutline.NextPos( iPt );
		if( !iNextPt )
			break;

        // Two points
        const Point3D & lV = pOutline.GetAt( iPt );
        const Point3D & lW = pOutline.GetAt( iNextPt );

		// Find projection on segment and absolute distance
		Point3D lProjP;
		double lDist = DistanceToSegment( pPt, lV, lW, lProjP );
		
		// Update
		if( lMinPosDist > lDist )
		{
			lMinPosDist = lDist;
			pProjP = lProjP;
		}
		
		// Two vectors
        const Vector3D lPtV = lV - pPt;
        const Vector3D lPtW = lW - pPt;

        // Check if point lies near one of the two vertices v or w
//        if( LzMath::ToolBox::IsZero( lPtV.Norm(), 1e-6 )
//         || LzMath::ToolBox::IsZero( lPtW.Norm(), 1e-6 ) )
        if( LzMath::ToolBox::IsZero(lPtV.Norm() * lPtW.Norm()) )
        {
            // Point on boundary
//return PointAndPolygon::Surface;
			return 0.0;
        }

        // Compute angle: it is safe since both vectors are not null
        const double lAngle = lPtV.SignedAngleTo( lPtW, lNorm );

        // Check if point lies on edge [v,w] iff. angle close to +Pi or -Pi
        if( LzMath::ToolBox::IsZero( lAngle - LzMath::PI, 1e-6 )
         || LzMath::ToolBox::IsZero( lAngle + LzMath::PI, 1e-6 ) )
        {
            // Point on boundary
//return PointAndPolygon::Surface;
			return 0.0;
        }

        // Point inside or outside
        lAngleSum += lAngle;
    }

    // Absolute value
    double lAbsAngleSum = std::abs( lAngleSum );

    // Out
    if( lAbsAngleSum < 1e-6 ) // Any tolerance lower than Pi/2 will do
	{
//return PointAndPolygon::Out;
		return +lMinPosDist;
	}
    else
	{
//return PointAndPolygon::In;
 		return -lMinPosDist;
	}
}

//================================================================================
double Outliner::SgnDistToClosedOuts( const List< List<Point3D> > & pOutlines, const Point3D & pPt, Point3D & pProjP/*, Vector3D & pProjN*/, const Vector3D * /*ppNormal=nullptr*/)
{
    double lMinAbsDist = std::numeric_limits<double>::max();
    double lMinSgnDist;
    BrowseList( iOut, pOutlines )
    {
        Point3D lProjP;
        double lSgnDist = SgnDistToClosedOut( pOutlines.GetAt(iOut), pPt, lProjP );

        double lAbsDist = fabs( lSgnDist );

        if( lMinAbsDist > lAbsDist )
        {
            lMinAbsDist = lAbsDist;
            lMinSgnDist = lSgnDist;
            pProjP = lProjP;
        }
    }

    return lMinSgnDist;
}

//================================================================================
//
// ** See: double FastDistComp::SignedDistance( Tree3D::NearestMode pTMode, const Point3D & pFrom, Point3D & pProjP, Vector3D & pProjN, SgnDstInfo * ppInfo/*=nullptr*/ ) const
//
//double Outliner::SgnDistToClosedOuts( const List< List<Point3D> > & pOutlines, const Point3D & pPt, Point3D & pProjP/*, Vector3D & pProjN*/ )
//{
//	double lMinAbsDist = std::numeric_limits<double>::max();
//	double lMinSgnDist;
//
//	BrowseList( iOut, pOutlines )
//	{
//		Point3D lProjP;
//		double lSgnDist = SgnDistToClosedOut( pOutlines.GetAt(iOut), pPt, lProjP );
//
//		double lAbsDist = fabs( lSgnDist );
//
//		if( lMinAbsDist > lAbsDist ) Check semantics !!! (outlines inside outlines ?? CCW vs CW orientation ??)
//		{
//			lMinAbsDist = lAbsDist;
//			lMinSgnDist = lSgnDist;
//			pProjP = lProjP;
//		}
//	}
//
//	return lMinSgnDist;
//}


#pragma region "CSG"
namespace // Hide names in local translation unit
{
//================================================================================
class CSNode
{
public:
    CSNode( const Point3D & pPos=Point3D()/*, size_t pIdx*/, bool pCut=false ) : mIsFinal(true), mIdx(-1/*pIdx*/), mPos(pPos), mCut(pCut) {}

	// Node init and relinking (when 2 outlines intersect at a node)
//void SetIdx( size_t pIdx )   { mFinalIdx = true;  mIdx = pIdx; }
	void RelinkTo( int pIdx )
	{
		// Useless: there should be no self-relink (where final becomes false, and idx = my index in sNods ==> infinite rabbit hole)
		if( mIdx == pIdx )
			return;

		// Relink
		mIsFinal = false;
		mIdx = pIdx;
	}

	// Relink
	bool mIsFinal;
	/*unsigned*/ int mIdx;

	// Index
	Point3D mPos;
	bool mCut;
};

//================================================================================
using LzServices::Sortable;
class CSEdge
{
public:
	// Sorted list: score = dist( cut, mN[0] )
    List< Sortable<size_t, double> > mCuts;
};
}

//================================================================================
static Vector<CSNode> sNods; //**** diviser en stack (taille statique) et heap (liste) dynamique pour les coupes
static Vector<CSEdge> sEdges;
//
static size_t sEdgNods_A;
static size_t sEdgNods_B;

//================================================================================
size_t FinalIdx( size_t pIdx )
{
	// Go down the rabbit hole...
	while( true )
	{
		CSNode & lNod = sNods[ pIdx ];
		if( lNod.mIsFinal )
			return pIdx;
		else
			pIdx = lNod.mIdx;
	}
}

//================================================================================
CSNode & FinalNode( size_t pIdx )
{
	// Go down the rabbit hole...
	while( true )
	{
		CSNode & lNod = sNods[ pIdx ];
		if( lNod.mIsFinal )
			return lNod;
		else
			pIdx = lNod.mIdx;

//***** bouger ce truc ds CSNode
//LzLogM("", "infinie final node")
	}
}

//================================================================================
static void AddCut( size_t pEdgIdx, size_t pCutIdx )
{
//const Point3D & lNod = sNods[ pEdgIdx ].mPos;
//const Point3D & lCut = sNods[ pCutIdx ].mPos;
	const Point3D & lNod = FinalNode( pEdgIdx ).mPos;
	const Point3D & lCut = FinalNode( pCutIdx ).mPos;

    Sortable<size_t, double> lInt(pCutIdx, lNod.DistanceTo(lCut));

	sEdges[ pEdgIdx ].mCuts.AddIntoIncrList( lInt, /*pUnique= */true );
}

//================================================================================
static bool InterEdges( size_t pEA, size_t pEB, double pTolerance )
{
	// Calc node indices
    const size_t iA[] = { pEA, pEA==sEdgNods_A-1            ? 0          : pEA+1 };
    const size_t iB[] = { pEB, pEB==sEdgNods_A+sEdgNods_B-1 ? sEdgNods_A : pEB+1 };

	// Read nodes
	const Point3D & lA0 = FinalNode( iA[0] ).mPos;
	const Point3D & lA1 = FinalNode( iA[1] ).mPos;
	//
	const Point3D & lB0 = FinalNode( iB[0] ).mPos;
	const Point3D & lB1 = FinalNode( iB[1] ).mPos;
//LzLogM("", "Final node ("<<iA[0]<<") = "<<FinalIdx(iA[0]))
//LzLogM("", "Final node ("<<iA[1]<<") = "<<FinalIdx(iA[1]))
//LzLogM("", "Final node ("<<iB[0]<<") = "<<FinalIdx(iB[0]))
//LzLogM("", "Final node ("<<iB[1]<<") = "<<FinalIdx(iB[1]))

	// Build lines
//const Line3D lLineA( lA0, lA1 );
//const Line3D lLineB( lB0, lB1 );
	Line3D lLineA, lLineB;
	try
	{
		lLineA = Line3D( lA0, lA1 );
		lLineB = Line3D( lB0, lB1 );
	}
	catch( ... )
	{
		// Ignore the intersection between those edges
		return false;
	}

	// Compute intersection
	Point3D lInter;
	if( lLineA.PseudoIntersection_NoExc( lLineB, lInter ) )
	{
		// Compute status for each edge
		// -2 = out of edge, -1 = edge cut, 0 = on node 0, 1 = on node 1

		int lSt_A = -2; // Invalid intersection
		{
			if( lInter.DistanceTo( lA0 ) <= pTolerance )
				lSt_A = 0; // On node 0
			else
			if( lInter.DistanceTo( lA1 ) <= pTolerance )
				lSt_A = 1; // On node 1
			else
			if( (lA0 - lInter)*(lA1 - lInter) < 0 )
				lSt_A = -1; // On edge
		}

		int lSt_B = -2; // Invalid intersection
		{
			if( lInter.DistanceTo( lB0 ) <= pTolerance )
				lSt_B = 0; // On node 0
			else
			if( lInter.DistanceTo( lB1 ) <= pTolerance )
				lSt_B = 1; // On node 1
			else
			if( (lB0 - lInter)*(lB1 - lInter) < 0 )
				lSt_B = -1; // On edge
		}

		// No intersection?
		if( lSt_A==-2 || lSt_B==-2 )
		{
			// One of the intersections is not valid
			return false;
		}
		else
		// Both edge cuts?
		if( lSt_A==-1 && lSt_B==-1 )
		{
			// Create new cut node
            size_t lCutIdx = sNods.Size();
			sNods.PushBack( CSNode(lInter, true) );

			// Insert intersection in both edges
			AddCut( pEA, lCutIdx );
			AddCut( pEB, lCutIdx );
		}
		else
		// Only one edge cut?
		if( lSt_A==-1 || lSt_B==-1 )
		{
			// Link to existing node

			if( lSt_A >= 0 )
			{
				// Intersection is node lSt_A on edge A

				// Mark node as cut
//sNods[ iFinA[lSt_A] ].mCut = true;
				FinalNode( iA[lSt_A] ).mCut = true;

				// Insert cut on edge B
				AddCut( pEB, iA[lSt_A] );
			}
			else
			{
				// Intersection is node lSt_B on edge B

				// Mark node as cut
//sNods[ iFinB[lSt_B] ].mCut = true;
				FinalNode( iB[lSt_B] ).mCut = true;

				// Insert cut on edge A
				AddCut( pEA, iB[lSt_B] );
			}
		}
		else
		// Node cuts only?
		{
			// Intersection is node lSt_A on edge A == node lSt_B on edge B

//// Mark both as cuts
//sNods[ iA[lSt_A] ].mCut = true;
//sNods[ iB[lSt_B] ].mCut = true;

			// No need to relink if both nodes are already linked to the same node
			if( FinalIdx(iA[lSt_A]) != FinalIdx(iB[lSt_B]) )
			{
				// Never relink a final node onto itself or you will create an infinite rabbit hole!

				// Create new cut node
                size_t lMidIdx = sNods.Size();
				const Point3D & lPtA = FinalNode( iA[lSt_A] ).mPos;
				const Point3D & lPtB = FinalNode( iB[lSt_B] ).mPos;

				// Push as final node
				sNods.PushBack( CSNode(lPtA.MidPointTo(lPtB), true) );

//if( pEA==30 && pEB==205 )
//	LzLogM("", "Hop ------------------------------------------")

//TxLogM(pEA<<" - "<<pEB<<": node cuts only: "<<iA[lSt_A]<<", "<<iB[lSt_B]<<" : relinking to "<<lMidIdx)
				// Relink
				FinalNode( iA[lSt_A] ).RelinkTo( lMidIdx );
//LzLogM("", "Now final index "<<iA[lSt_A]<<" = "<<FinalIdx( iA[lSt_A] ))
//LzLogM("", "--- Final index "<<iB[lSt_B]<<" = "<<FinalIdx( iB[lSt_B] ))

	CSNode & lFinB = FinalNode( iB[lSt_B] );

				lFinB.RelinkTo( lMidIdx );
//LzLogM("", "---- Now final index "<<iB[lSt_B]<<" = "<<FinalIdx( iB[lSt_B] ))
			}
		}

		return true;
	}
	else
	{
		// Overlap?

// TODO
//
//**** NO INTERSECTION ****
//
// TODO
		return false;
	}
}

//================================================================================
static void HarvestStrips( size_t pFrom, size_t pTo, List< List<size_t> > & pStrips )
{
	// Automaton for strip extraction
	enum class State { NOP, ADD };
	
	// Unwrap all indices in A - cuts included, without repeating first node
    List<size_t> lUnwrap;
    for( size_t e_idx=pFrom ; e_idx<pTo ; e_idx++ )
	{
		// Remove repeating indices, in case node relinking produces and edge with identical start and end nodes

		// Add first node
        size_t lFinIdx = FinalIdx( e_idx );
		if( lUnwrap.Count()==0 || lFinIdx != lUnwrap.GetTail() )
			lUnwrap.AddTail( FinalIdx( e_idx ) );

		// Add cuts
        const List< Sortable<size_t, double> > & lCuts = sEdges[ e_idx ].mCuts;
		BrowseList( iC, lCuts )
		{
//*** DEBUG
if( FinalIdx( lCuts.GetAt(iC).mT ) >= sNods.Size() )
    LzLogException("", "OVERFLOW 1")
//*** DEBUG

            size_t lFinIdx = FinalIdx( lCuts.GetAt(iC).mT );
			if( lFinIdx != lUnwrap.GetTail() )
				lUnwrap.AddTail( lFinIdx );
		}
	}

////**** TODO: Remove repeating indices, 
////**** TODO: Remove repeating indices, in case node relinking produces and edge with identical start and end nodes
////
//void * iU0 = lUnwrap
//BrowseList( iU0, lUnwrap )
//{
//	void * iU1 = lUnwrap.NextPos( iU0 );
//	if( iU1 == nullptr )
//		break;
//
//	if( lUnwrap.GetAt( iU0 ) == lUnwrap.GetAt( iU1 ) )
//		LzLogM("", "**** FOUND repeating index "<<lUnwrap.GetAt( iU0 ))
//}
////
////**** TODO: Remove repeating indices, in case node relinking produces and edge with identical start and end nodes
////**** TODO: Remove repeating indices, in case node relinking produces and edge with identical start and end nodes

//*** DEBUG
BrowseList( iU, lUnwrap )
	if( lUnwrap.GetAt(iU) >= sNods.Size() )
        LzLogException("", "OVERFLOW")
//*** DEBUG

	//
	// lUnwrap contains only final indices
	//
		
	// Set automaton
	State lState = State::NOP;
	void * lStartN;
	void * iN = lUnwrap.HeadPos();
	while( true )
	{
		// Current node
        size_t lIdx = lUnwrap.GetAt( iN );
			
		// Found a cut?
		if( sNods[ lIdx ].mCut ) // No need to use FinalNode as lUnwrap contains only final indices
		{
			// Start harvesting?
			if( lState == State::NOP )
			{
				// Init new strip
				lState = State::ADD;

				// Start pos?
				if( pStrips.Count() == 0 )
					lStartN = iN;

				// Init new strip
				pStrips.AddTail( { lIdx } );

				// Next node
				lUnwrap.Next( iN );
				if( iN == nullptr )
					iN = lUnwrap.HeadPos();
			}
			// Stop harvesting?
			else
			{
				// Finish current strip
				pStrips.GetTail().AddTail( lIdx );
				lState = State::NOP;

				// Back at start?
				if( iN == lStartN )
					break;
			}
		}
		else
		{
			// Harvesting a strip?
			if( lState == State::ADD )
				pStrips.GetTail().AddTail( lIdx );

			// Next node
			lUnwrap.Next( iN );
			if( iN == nullptr )
				iN = lUnwrap.HeadPos();
		}
	}
}

//================================================================================
static void AddFilterStrips( Outliner::Operator /*pOp*/,
                          const List< List<size_t> > & pFrom,
						  const List<Point3D> & pOut,
                          List< List<size_t> > & pTo )
{
	//// Del previous
	//pTo.DelAll(); --> nope, method accumulates strips

	// Filter strips
	BrowseList( iS, pFrom )
	{
		// Get the strip
        const List<size_t> & lS = pFrom.GetAt( iS );

		// Compute total length-weighted signed distance
		double lWeiSgnDst = 0.0;
		{
			void * iP0 = lS.HeadPos();
			while( true )
			{
				// Finished ?
				void * iP1 = lS.NextPos( iP0 );
				if( iP1 == nullptr )
					break;
			
				// Edge points
				const Point3D & lP0 = sNods[ lS.GetAt(iP0) ].mPos;
				const Point3D & lP1 = sNods[ lS.GetAt(iP1) ].mPos;

				// Edge mid
				const Point3D lMid = lP0.MidPointTo( lP1 );

				// Signed distance
				Point3D lTmp;
				double lSgnDst = Outliner::SgnDistToClosedOut( pOut, lMid, lTmp );

				// Accum
				lWeiSgnDst += lP0.DistanceTo(lP1) * lSgnDst;

				// Next
				iP0 = iP1;
			}
		}

//*** DEBUG
//LzLogM("", "lWeiSgnDst = "<<lWeiSgnDst)
//*** DEBUG

		// Keep only negative outlines
		if( lWeiSgnDst <= 0.0 )
			pTo.AddTail( lS );
	}
}

//================================================================================
static void AddCoordsOutline( const List<size_t> & pOut, List< List<Point3D> > & pTo )
{
	// Do not use Append to avoid copying lists
	pTo.AddTail( {} );
	BrowseList( iO, pOut )
		pTo.GetTail().AddTail( sNods[ pOut.GetAt(iO) ].mPos );
}

//================================================================================
void Outliner::CSG( Operator pOp,
					const List<Point3D> & pA, const List<Point3D> & pB,
					List< List<Point3D> > & pC,
					double pTolerance/*=1e-6*/,
					CSG_debug_info * ppDebugInfo/*=nullptr*/  )
{
	// Check
	if( !IsClosed(pA) || !IsClosed(pB) )
        LzLogException("", "Outlines must be closed to perform CSG!")

	// Check distance between consecutive nodes
	//
	// *** TODO
	//
	// Outline edges must be longer than tolerance
	// Outlines must be non-self-intersecting (within tolerance)
	//
	// *** TODO

	// Clean debug data
	if( ppDebugInfo )
		ppDebugInfo->Free();

	// Clean previous
	pC.DelAll();

	// Set count of nodes and edges in each outline
	sEdgNods_A = pA.Count() - 1; // Last is duplicate
	sEdgNods_B = pB.Count() - 1; // Last is duplicate
    const size_t lTot_EdgNods = sEdgNods_A + sEdgNods_B;

	//---------------------------------------
	// Set initial nodes
	//---------------------------------------

	sNods.Resize( lTot_EdgNods );
////////////////////////////////sNods.SetCapacity( lTot_EdgNods );
	{
        size_t iN = 0;

		// From A
		BrowseList( iP, pA )
		{
			// Skip last point: repeated first point
			if( pA.NextPos(iP) == nullptr )
				break;

			// Store node and set final index
			sNods[ iN ] = pA.GetAt( iP );
//sNods[ iN ].SetIdx( iN );
			iN++;
		}

		// From B
		BrowseList( iP, pB )
		{
			// Skip last point: repeated first point
			if( pB.NextPos(iP) == nullptr )
				break;

			// Store node and set final index
			sNods[ iN ] = pB.GetAt( iP );
//sNods[ iN ].SetIdx( iN );
			iN++;
		}
	}

	//---------------------------------------
	// Create edges
	//---------------------------------------

	sEdges.Resize( lTot_EdgNods );
    for( size_t e=0 ; e<lTot_EdgNods ; e++ )
		sEdges[e].mCuts.DelAll();

	//---------------------------------------
	// Intersect edges
	//---------------------------------------

	bool lHaveCuts = false;
    for( size_t ea=0          ; ea<sEdgNods_A            ; ea++ )
    for( size_t eb=sEdgNods_A ; eb<sEdgNods_A+sEdgNods_B ; eb++ )
		lHaveCuts |= InterEdges(ea, eb, pTolerance);

//***** DEBUG
if( ppDebugInfo )
{
    ppDebugInfo->mNods.resize( sNods.Size() );
    for( size_t n=0 ; n<sNods.Size() ; n++ )
	{
		ppDebugInfo->mNods[n] = sNods[n].mPos;
		if( sNods[n].mCut )
			ppDebugInfo->mCuts.AddTail( n );
	}
}
//***** DEBUG

	//---------------------------------------
	// Extract strips
	//---------------------------------------
	//
	// Strips start with a cut index and finish with a cut index
	// In between, there is no cut index
	//
	// There can be more strips than cuts e.g. in star shaped contours
	// if node relinking has brought all center nodes together
	//
    List< List<size_t> > lStrips_A, lStrips_B;
	if( lHaveCuts )
	{
		// Need to cut initial outlines into strips
		HarvestStrips( 0,          sEdgNods_A,            lStrips_A );
		HarvestStrips( sEdgNods_A, sEdgNods_A+sEdgNods_B, lStrips_B );
	}
	else
	{
		// No cuts: strips are initial outlines

		// A
		lStrips_A.AddTail( {} );
        for( size_t n=0 ; n<sEdgNods_A ; n++ )
			lStrips_A.GetTail().AddTail( n );
		//
		// Close strip
		lStrips_A.GetTail().AddTail( 0 );

		// B
		lStrips_B.AddTail( {} );
        for( size_t n=sEdgNods_A ; n<lTot_EdgNods ; n++ )
			lStrips_B.GetTail().AddTail( n );
		//
		// Close strip
		lStrips_B.GetTail().AddTail( sEdgNods_A );
	}

//***** DEBUG
if( ppDebugInfo )
{
	ppDebugInfo->mStrips_A = lStrips_A;
	ppDebugInfo->mStrips_B = lStrips_B;
}
//***** DEBUG

	//---------------------------------------
	// Compute signed distances
	//---------------------------------------

	//---------------------------------------
	// Filter strips (inverted if needed)
	//---------------------------------------

	//
	// *** TODO: CHANGE ORIENTATION of strips if substracting
	//
	//
	// *** TODO: REMOVE duplicate strips before filtering and stitching
	// ***       duplicate strips are the ones with same begin and end index
	//

    List< List<size_t> > lFilStr_A, lFilStr_B;
	AddFilterStrips( pOp, lStrips_A, pB, lFilStr_A );
	AddFilterStrips( pOp, lStrips_B, pA, lFilStr_B );

//***** DEBUG
if( ppDebugInfo )
{
	ppDebugInfo->mFilStr_A = lFilStr_A;
	ppDebugInfo->mFilStr_B = lFilStr_B;
}
//***** DEBUG

	// Merge lists
    List< List<size_t> > lFiltrStr;
	lFiltrStr.Append( lFilStr_A );
	lFiltrStr.Append( lFilStr_B );

	//---------------------------------------
	// Stitch strips into final contour
	//---------------------------------------


//LzLogM("", "*** Total filtered strips= "<<lFiltrStr.Count())


	// Remove all closed strips
	void * iS = lFiltrStr.HeadPos();
	while( iS )
	{
		// Get strip
        const List<size_t> & lS = lFiltrStr.GetAt( iS );

		// Is already closed and not degenerate? (needs at least 3 points + 1 repeating for closure)
		if( lS.GetHead()==lS.GetTail() && lS.Count()>=4 )
		{
			// Convert to coordinates and add to output
			AddCoordsOutline( lS, pC );

LzLogM("", "Stashing self-closed contour "<<lS.ListToString())

			// Remove from strips
			lFiltrStr.DelAtAndNext( iS );
		}
		else
			lFiltrStr.Next( iS );
	}


//LzLogM("", "*** Closed strips= "<<pC.Count())


	// If no more strips, my job here is done
	if( lFiltrStr.Count() == 0 )
		return;

	// Try stitching remaining strips together
    List<size_t> lCurrOut;
	enum class State { NOP, ADD };

	State lState = State::NOP;
	bool lChange;
	do
	{
		// Assume that nothing will change
		lChange = false;

		// Check all strips
		BrowseList( iS, lFiltrStr )
		{
			// Get strip
            List<size_t> & lS = lFiltrStr.GetAt( iS );

			// Not doing anything right now?
			if( lState == State::NOP )
			{
				// Try stitching this strip
				lCurrOut = lS;
				lFiltrStr.DelAtAndNext( iS );
				lState = State::ADD;

				// Things are changing
				lChange = true;
			}
			else
			{
				// Try to stitch this strip to current strip
				if( lCurrOut.GetTail() == lS.GetHead() )
				{
					// Stitch S after current
					lS.DelHead();
					lCurrOut.Append( lS );
					lFiltrStr.DelAtAndNext( iS );

					// Things are changing
					lChange = true;
				}
				else
				if( lCurrOut.GetHead() == lS.GetTail() )
				{
					// Stitch S before current
					lS.DelTail();
					lCurrOut.Prepend( lS );
					lFiltrStr.DelAtAndNext( iS );

					// Things are changing
					lChange = true;
				}

				// Completed a closed outline?
				if( lChange && lCurrOut.GetHead()==lCurrOut.GetTail() )
				{
					// Is not degenerate? (needs at least 3 points + 1 repeating for closure)
					if( lCurrOut.Count() >= 4 )
					{
						// Convert to coordinates and add to output
						AddCoordsOutline( lCurrOut, pC );
					}

LzLogM("", "Stashing contour "<<lCurrOut.ListToString())

					// Mark as finished
					lCurrOut.DelAll();
					lState = State::NOP;
				}
			}
		}

		// Finished loop but couldn't add anything to current strip? Drop it!
		if( lState==State::ADD && !lChange )
		{
//LzLogM("", "Dropping contour "<<lCurrOut.ListToString())

			lCurrOut.DelAll();
			lState = State::NOP;
		}
	}
	while( lChange );
}
#pragma endregion


//================================================================================
void Outliner::NormalOffsetOutline( List<Point3D> & pOutline, double pOffset )
{
    // Check
    if( !IsClosed(pOutline) )
        LzLogException("", "Cannot apply offset to an open outline!")

    // Check
    if( pOutline.Count() < 4 )
        LzLogException("", "Cannot apply offset to an outline with less than 3 points (+ 1 repeated head point)!")

    // Compulte LSQ plane
    Plane3D lPlane( pOutline );

    // Get normal
    const Vector3D & lNor = lPlane.Normal();

    // Find orientation
    Orientation lOr = GetOrientation( pOutline, lNor );

    // Offset points
    BrowseList( iP, pOutline )
    {
        // Break if last point
        if( iP == pOutline.TailPos() )
            break;

        // Compute local outgoing normal
        Vector3D lOutNor;
        {
            // Find point before
            void * lPos0;
            {
                if( iP == pOutline.HeadPos() )
                    lPos0 = pOutline.PrevPos( pOutline.TailPos() );
                else
                    lPos0 = pOutline.PrevPos( iP );
            }

            // Find point after
            void * lPos1 = pOutline.NextPos( iP ); // iP is not allowed to reach TailPos

            // Get points
            Vector3D AB = pOutline.GetAt( lPos1 ) - pOutline.GetAt( lPos0 );
            AB.Normalize();

            // Compute normal
            if( lOr == Orientation::CCW )
                lOutNor = AB ^ lNor;
            else
                lOutNor = lNor ^ AB;
        }

        // Move point
        pOutline.GetAt( iP ) += pOffset * lOutNor;
    }

    // Move tail point
    pOutline.GetTail() = pOutline.GetHead();
}


#pragma region "Crop"
/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------
 *Point3D Bisect( const Point3D & pA, const bool pKeepA, const Point3D & pB, const bool pKeepB, std::function<bool(const Point3D & pPt)> pKeep, double pEpsilon )
 *{
 *    // Check
 *    if( pKeepA == pKeepB )
 *        LzLogException("", "A and B cannot both be in or out!")
 *
 *    // Create moving copies
 *    Point3D lA = pA;
 *    Point3D lB = pB;
 *
 *    while( lA.DistanceTo(lB) > pEpsilon )
 *    {
 *        // Evaluate mid point
 *        Point3D lC = lA.MidPointTo( lB );
 *        bool lKeepC = pKeep(lC);
 *
 *        // Move A or B
 *        if( lKeepC == pKeepA )
 *            lA = lC;
 *        else
 *        if( lKeepC == pKeepB )
 *            lB = lC;
 *    }
 *
 *    return lA.MidPointTo( lB );
 *}*/

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------
 *void Outliner::Crop( const List<Point3D> & pOut, List< List<Point3D> > & pTo, std::function<bool(const Point3D & pPt)> pKeep )
 *{
 *    LzLogException("", "*** Outliner::Crop: TODO!")
 *
 *    // Local epsilon
 *    static const double lEps = 1e-6;
 *
 *    // Clean previous
 *    pTo.DelAll();
 *
 *    // Check
 *    if( pOut.Count() == 0 )
 *        return;
 *
 *    // Create local copy
 *    List<Point3D> lOut = pOut;
 *
 *    // Add bisect points
 *    void * iP = lOut.HeadPos();
 *    while( iP )
 *    {
 *        // Curr and next point
 *        void * iNextP = lOut.NextPos( iP );
 *        if( !iNextP )
 *            break;
 *
 *        // Curr and next point
 *        const Point3D lA = lOut.GetAt( iP );
 *        const Point3D lB = lOut.GetAt( iNextP );
 *        //
 *        bool lKeepA = pKeep(lA);
 *        bool lKeepB = pKeep(lB);
 *
 *        // Have intersection
 *        if( lKeepA != lKeepB )
 *        {
 *            const Point3D lC = Bisect( lA, lKeepA, lB, lKeepB, pKeep, lEps );
 *
 *            // Insert point only if far enough from both ends
 *            if( lC.DistanceTo(lA)>lEps && lC.DistanceTo(lB)>lEps )
 *                lOut.AddBefore( iNextP, lC );
 *
 *LzLogM("", "Added bisect point")
 *        }
 *
 *        // Next point
 *        iP = iNextP;
 *    }
 *
 *pTo.AddTail( pOut );
 *pTo.GetHead().DelTail();
 *
 *}*/

//================================================================================
/*
#if 1 // TEST CODE
{
    static const double DEG_2_RAD = LzMath::PI / 180.0;

    //////////////////////////////////////////////////////
    if( 0 )
    {
        LzLogN("", "Test: circle")

        List<Point3D> lOut;
        {
            for( size_t a=0 ; a<360 ; a+=45 )
                lOut.AddTail( Point3D( 10*cos(DEG_2_RAD*a), 10*sin(DEG_2_RAD*a), 0 ) );

            // Close outline
            lOut.AddTail( lOut.GetHead() );
        }

        // Clip
        List< List<Point3D> > lClips;
        Outliner::Crop( lOut, Plane3D(Point3D(0,0,0), Vector3D(0,1,0)), lClips );
//        Outliner::Crop( lOut, Plane3D(Point3D(0,0,0), Vector3D(1,0,0)), lClips );
//        Outliner::Crop( lOut, Plane3D(Point3D(5,0,0), Vector3D(1,0,0)), lClips );
//        Outliner::Crop( lOut, Plane3D(Point3D(0,5,0), Vector3D(0,1,0)), lClips );
//        Outliner::Crop( lOut, Plane3D(Point3D(-20,0,0), Vector3D(1,0,0)), lClips );
//        Outliner::Crop( lOut, Plane3D(Point3D(+20,0,0), Vector3D(1,0,0)), lClips );

        // Log
        LzLogM("", lClips.Count()<<" clip(s)")
        BrowseList( iC, lClips )
        {
            LzLogN("", "Clip")
            const List<Point3D> & lCl = lClips.GetAt( iC );
            BrowseList( iP, lCl )
            {
                LzLogM("", lCl.GetAt(iP).ToString())
            }
        }
    }
    //////////////////////////////////////////////////////
    if( 1 )
    {
        LzLogN("", "Test: star")

        List<Point3D> lOut;
        {
            lOut.AddTail( Point3D(+5,  0, 0) );
            lOut.AddTail( Point3D(+1, +1, 0) );
            lOut.AddTail( Point3D( 0, +5, 0) );
            lOut.AddTail( Point3D(-1, +1, 0) );
            lOut.AddTail( Point3D(-5,  0, 0) );
            lOut.AddTail( Point3D(-1, -1, 0) );
            lOut.AddTail( Point3D( 0, -5, 0) );
            lOut.AddTail( Point3D(+1, -1, 0) );
            // Close outline
            lOut.AddTail( lOut.GetHead() );
        }

        // TEST: open outline
        lOut.DelTail();

        // Clip
        List< List<Point3D> > lClips;
//        Outliner::Crop( lOut, Plane3D(Point3D(0,0,0), Vector3D(0,1,0)), lClips );
//        Outliner::Crop( lOut, Plane3D(Point3D(0,0.5,0), Vector3D(0,1,0)), lClips );
        Outliner::Crop( lOut, Plane3D(Point3D(2,2,0), Vector3D(1,1,0)), lClips );

        // Log
        LzLogM("", lClips.Count()<<" clip(s)")
        BrowseList( iC, lClips )
        {
            LzLogN("", "Clip")
            const List<Point3D> & lCl = lClips.GetAt( iC );
            BrowseList( iP, lCl )
            {
                LzLogM("", lCl.GetAt(iP).ToString())
            }
        }
    }

    return;
}
#endif
*/
void Outliner::Crop( const List<Point3D> & pOutline, const Plane3D & pCutPlane, List< List<Point3D> > & pRes, double pTol/*=1e-6*/ )
{
    // Clear previous
    pRes.DelAll();

    //-------------------------------
    // Slice list and assign status
    //-------------------------------

    List<Point3D> lPoints;
    List<double> lDists;
    enum class Status { In, Surf, Out };
    List<Status> lStats;
    //
    BrowseList( iP, pOutline )
    {
        // Get point
        const Point3D & lP1 = pOutline.GetAt( iP );

        // Calc signed dist and status
        double lD1;
        Status lS1;
        {
            const double lSgnD1 = pCutPlane.SignedDistanceTo( lP1 );
            lD1 = fabs( lSgnD1 );
            {
                if( lD1 <= pTol )
                    lS1 = Status::Surf;
                else
                if( lSgnD1 < 0 )
                    lS1 = Status::Out;
                else
                    lS1 = Status::In;
            }
        }

        // Surface or first point?
        if( lS1==Status::Surf || iP==pOutline.HeadPos() )
        {
            // Store current point
            lPoints.AddTail( lP1 );
            lDists.AddTail( lD1 );
            lStats.AddTail( lS1 );
        }
        else
        {
            //--> NOT surface AND NOT first point

            // Read previous status
            Status lS0 = lStats.GetTail();

            // Have a transition?
            if( lS0!=Status::Surf && lS0!=lS1 )
            {
                // Read previous point and absolute distance
                const Point3D & lP0 = lPoints.GetTail();
                const double lD0 = lDists.GetTail();

                // Store inter point
                lPoints.AddTail( lP0 + (lD0 / (lD0 + lD1))*(lP1 - lP0) );
                lDists.AddTail( 0.0 );
                lStats.AddTail( Status::Surf );
            }

            //--> Code above is skipped if: previous status == surface or same status as current
            //--> Assuming ONLY ONE surface point between two consecutive In/Out points

            // Store current point
            lPoints.AddTail( lP1 );
            lDists.AddTail( lD1 );
            lStats.AddTail( lS1 );
        }
    }


    //-------------------------------
    // Extract strips
    //-------------------------------

    List<Point3D> lCurrPoints;
    List<Status> lCurrStats; // Will contain only S and I statuses
    BrowseTwoLists( iP, lPoints, iS, lStats )
    {
        // Read current status
        Status lS = lStats.GetAt( iS );

        // Transitions
        if( lS == Status::Out )
        {
            // Skip
            continue;
        }
        else
        if( lS == Status::In )
        {
            // Add to current strip
            lCurrPoints.AddTail( lPoints.GetAt(iP) );
            lCurrStats.AddTail( Status::In );
        }
        else
        // lS == Status::Surf
        {
            // Have non empty current strip?
            if( lCurrPoints.Count() )
            {
                //--> Yes, curr points is not empty

                // Check if transitionning from S to S or I to S
                if( lCurrStats.GetTail() == Status::Surf )
                {
                    //--> S to S

                    // Replace previous S with more recent S
                    lCurrPoints = { lPoints.GetAt(iP) };
                    lCurrStats = { Status::Surf };
                }
                else
                {
                    //--> I to S

                    // Add finish point
                    lCurrPoints.AddTail( lPoints.GetAt(iP) );
                    lCurrStats.AddTail( Status::Surf );

                    // Stash
                    pRes.AddTail( lCurrPoints );

                    // Purge current
                    lCurrPoints.DelAll();
                    lCurrStats.DelAll();
                }
            }
            else
            {
                //--> Nope, curr points is empty

                // Init new strip with Surf point
                lCurrPoints = { lPoints.GetAt(iP) };
                lCurrStats = { Status::Surf };
            }
        }
    }

    // Stash last (running) outline, if not reduced to a Surf (starting) point
//    if( lCurrPoints.Count() && lCurrStats.GetTail()==Status::In )
    if( lCurrPoints.Count() > 1 )
        pRes.AddTail( lCurrPoints );


    //-------------------------------
    // Merge closed strip
    //-------------------------------

    if( IsClosed(pOutline) && pRes.Count()>=2 && lStats.GetHead()==Status::In )
    {
        //--> First strip passing through head point, and another strip passing through (the same) tail point

        // Remove last point from tail strip
        pRes.GetTail().DelTail();

        // Append head strip to tail strip
        pRes.GetTail().Append( pRes.GetHead() );

        // Remove head strip
        pRes.DelHead();
    }
}
#pragma endregion


//================================================================================
Point3D Outliner::FindCurviPoint( const List<Point3D> & pOutline, double pS )
{
    // Check
    if( pOutline.Count() == 0 )
        LzLogException("", "Cannot find curvilinear point on an empty outline!")

    // Check
    if( pS < 0 )
        LzLogException("", "Cannot find point for a negative curvilinear absciss!")

    double lSoFar = 0.0;
    BrowseList( iP0, pOutline )
    {
        // Next pos
        void * iP1 = pOutline.NextPos( iP0 );

        // Reached end of list?
        if( !iP1 )
            break;

        // Get points
        const Point3D & lP0 = pOutline.GetAt( iP0 );
        const Point3D & lP1 = pOutline.GetAt( iP1 );
        const double lLen01 = lP0.DistanceTo( lP1 );

        // Within reach?
        if( lSoFar + lLen01 >= pS )
        {
            // Yes!

            // Compute remainder
            const double lRemS = pS - lSoFar; // RemS in [0, Len01]

            // Unit dir
            Vector3D lP0P1 = lP1 - lP0;
            lP0P1.Normalize();

            // Compute point
            return lP0 + lRemS*lP0P1;
        }
        else
        {
            // No
            lSoFar += lLen01;
        }
    }

    // Reached end of list before finding curvilinear abscissa
    return pOutline.GetTail();
}

//================================================================================
void Outliner::UniformResample( const List<Point3D> & pOutline, size_t pNbSegs, List<Point3D> & pRes )
{
    // Check
    if( &pOutline == &pRes )
        LzLogException("", "Cannot resample outline to itself!")

    // Check
    if( pOutline.Count() < 2 )
        LzLogException("", "Cannot resample outline with less than 2 points!")

    // Check
    if( pNbSegs < 1 )
        LzLogException("", "Cannot resample outline into less than 1 segment!")

    // Compute step
    const double lStep = OutlineLength( pOutline ) / pNbSegs;

    // Clear previous
    pRes.DelAll();

    // Add start
    pRes.AddTail( pOutline.GetHead() );

    // Add end segments
    for( size_t s=1 ; s<pNbSegs ; s++ )
        pRes.AddTail( FindCurviPoint(pOutline, s*lStep) );

    // Add end
    pRes.AddTail( pOutline.GetTail() );
}
#pragma endregion
}
