#include "Tree3D.h"
#include "Plane3D.h"
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>
#if defined(USING_QT) && !defined(NO_OPENGL)
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/gl.h>
#endif


namespace LzGeom
{
#pragma region "Construction / destruction"
//================================================================================
Tree3D::Tree3D() : mpRootNode( nullptr ) /*: mNbLeaf(0), mConcurrentLeaves(0)*/
{
}

//================================================================================
Tree3D::~Tree3D()
{
    Free();
}

//================================================================================
void Tree3D::Free()
{
    delete mpRootNode;
    mpRootNode = nullptr;

    mPoints.Free();
    mNew2Old.Free();
    mOld2New.Free();
}
#pragma endregion


#pragma region "Set-up"
//================================================================================
void Tree3D::Create( const vector<Point3D> & pPoints, bool pLogStats/*=true*/, double pMinDiameter/*=0.001*/, unsigned int pMaxLeafCount/*=10*/, unsigned int pMaxTreeDepth/*=10*/ )
{
    // Vector from vector
    Vector<Point3D> lPts( (unsigned int)pPoints.size() );
    for( unsigned int p = 0 ; p < pPoints.size() ; p++ )
        lPts[p] = pPoints[p];

    // Go
    Create( lPts, pLogStats, pMinDiameter, pMaxLeafCount, pMaxTreeDepth );
}

//================================================================================
void Tree3D::Create( const Vector<Point3D> & pPoints, bool pLogStats/*=true*/, double pMinDiameter/*=0.001*/, unsigned int pMaxLeafCount/*=10*/, unsigned int pMaxTreeDepth/*=10*/ )
{
    // Check
    if( !pPoints.Size() )
        LzLogException("",  "Cannot build a tree on empty points cloud!" );

    // Free previous
    Free();

try
{

    // Log
    LzLogTimeNode("", "Creating 3D-tree")

    // Copy points localy for later spatial quick-sorting
    mPoints = pPoints;

    // Create vector of indices
    mNew2Old.Resize( mPoints.Size() );
    for( unsigned int p=0 ; p<mPoints.Size() ; p++ )
        mNew2Old[p] = p;

	// Statistics
    mNbLeaf = 0;
    mSizeMax = 0;
    mSizeMin = pMaxLeafCount;
    mDepthMax = 0;

    // Create tree
    mpRootNode = new Node;
    mpRootNode->UpdateTree( mPoints,								// Data points
							mNew2Old,								// Mapping between new (quick sorted) and old (initial vertices)
							mDepthMax, mNbLeaf, mSizeMax, mSizeMin,
							0, mPoints.Size()-1,					// Begin - End
							0,										// Curr depth
							pMaxLeafCount, pMinDiameter, pMaxTreeDepth );

    // Compute inverse index mapping
    mOld2New.Resize( mPoints.Size() );
    for( unsigned int p = 0 ; p < mPoints.Size() ; p++ )
        mOld2New[ mNew2Old[p] ] = p;

// Statistics
    mMoyLeaf = ( double )pPoints.Size() / ( double )mNbLeaf;

    //code couleur

//// Change indices of triangles
//for (unsigned int iTri = 0; iTri < ppMesh->mTriangles.size(); iTri++)
//{
//  for (unsigned int iVertex = 0; iVertex < 3; iVertex++)
//  {
//      for (unsigned int iV = 0; iV < mpPointsIndices->size(); iV++)
//          if (mpPointsIndices->at(iV) == ppMesh->mTriangles[iTri].mIdxV[iVertex])
//          {
//              ppMesh->mTriangles[iTri].mIdxV[iVertex] = iV;
//              break;
//          }
//  }
//}

// Compute topology
//mTopo.SetMesh( *ppMesh );

    // Log
    if( pLogStats )
    {
        LzLogN("", "Statistics")
        LzLogM("", "Depth max.= " << mDepthMax << ".")
        LzLogM("", "Average points per leaf= " << mMoyLeaf << ".")
        LzLogM("", "" << mNbLeaf << " leave(s).")
        LzLogM("", "Leaf count min.= " << mSizeMin << ".")
        LzLogM("", "Leaf count max.= " << mSizeMax << ".")
    }
    LzLogM("", "Finished.")
}
catch ( const std::exception & e )
{
    LzLogErr("", e.what())
    Free();
    LzLogException("", "Could not init 3D-tree!");
}
}
#pragma endregion


#pragma region "Computation"
//================================================================================
unsigned int Tree3D::GetNearestVertex( const Point3D & pPt, NearestMode pMode ) const
{
    // Returns final leaf containing actual data
    const Node * lpNode = mpRootNode->GetNearestLeaf( pPt );

    // Find nearest vertex within final leaf
    unsigned int lNearest = lpNode->GetNearestPoint( pPt, mPoints );

    // Explore other potential candidates
    double lMinDist = pPt.DistanceTo( mPoints[lNearest] );

    // Wanna go fast?
    if( pMode == NearestMode::Fast )
        return mNew2Old[ lNearest ];

    // Nope, let's be accurate!
//************************************* reset statistics ????????????????????????????????????????????????????????????????
//mConcurrentLeaves = 0;
unsigned int mConcurrentLeaves = 0;
    mpRootNode->GetNextLeaf( mConcurrentLeaves, lNearest, lMinDist, pPt, mPoints );
//************************************* reset statistics ????????????????????????????????????????????????????????????????

    // Convert nearest index from qsorted to initial index
    return mNew2Old[ lNearest ];
}

//================================================================================
void Tree3D::GetRangeVertices( const Point3D & pPt, double pRange, List<size_t> & pRangeIdx ) const
{
	// Clean previous
	pRangeIdx.DelAll();

	// Collect range vertices
	mpRootNode->CollectRangeVertices( pPt, pRange, pRangeIdx, mPoints, mNew2Old );
}

//================================================================================
const Point3D & Tree3D::OldPoint3D( size_t pOldIdx ) const
{
    // Check
    if( pOldIdx > mOld2New.Size() )
        LzLogException("",  "Index " << pOldIdx << " is out of range! " << mOld2New.Size() << " vertice(s) in table." );

    return mPoints[ mOld2New[pOldIdx] ];
}
#pragma endregion


#pragma region "Debug"
#if defined(USING_QT) && !defined(NO_OPENGL)
//================================================================================
void Tree3D::DrawLeaves() const
{
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    glDisable( GL_LIGHTING );
    glLineWidth( 2 );

    if( mpRootNode )
        mpRootNode->DrawLeaves( mMoyLeaf, mSizeMax, mSizeMin );

    glPopAttrib();
}

//================================================================================
void Tree3D::DrawPlanes( /* bool pIsMeters */ ) const
{
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    glEnable( GL_LIGHTING );
    glDisable( GL_CULL_FACE );

    if( mpRootNode )
        mpRootNode->DrawPlanes();

    glPopAttrib();
}

//================================================================================
void Tree3D::Draw2D( double pWidth, double pHeight ) const
{
    if( mpRootNode )
    {
        glPushAttrib( GL_ALL_ATTRIB_BITS );
        glDisable( GL_LIGHTING );

        glLineWidth( 1 );
        glPointSize( 5 );

        //const Point3D lRoot = ;

        //if( mDepthMax )
        mpRootNode->DrawTree2D( Point3D( pWidth / 2.0, pHeight, 0 ), pWidth / 2.0, -pHeight / mDepthMax );
        //else
        //  mpRootNode->DrawTree2D( lRoot, 0, 0 );

        glPopAttrib();
    }
}
#endif
#pragma endregion


//================================================================================
//double Tree3D::GetDistance(const Point3D & pPt, const Mesh * ppMesh, Point3D & pBestPoint)
//{
//  const unsigned int lNearest = GetNearestVertex(pPt);
//  double lMinDist = pPt.DistanceTo(mpPoints->at(lNearest));
//
//  // Find triangles
//  const List<int> & lNeighborTris = mTopo.TopVers()[ lNearest ].mT;
//  BrowseList( iT, lNeighborTris )
//  {
//      const LzTriModel::TopTri & lTri = mTopo.TopTris()[ lNeighborTris.GetAt(iT) ];
//
//      const Point3D & lA = mpPoints->at( lTri.mV[0] );
//      const Point3D & lB = mpPoints->at( lTri.mV[1] );
//      const Point3D & lC = mpPoints->at( lTri.mV[2] );
//
//      Point3D lProjP;
//      double lDist = LzTriModel::DistanceComputer::PositiveDistanceToTriangle( pPt, lA, lB, lC, lProjP);
//      if(lDist < lMinDist)
//      {
//          pBestPoint = lProjP;
//          lMinDist = lDist;
//      }
//  }
//  return lMinDist;
//}

//================================================================================
//void Tree3D::DrawCandidatesTri(const unsigned int pNearest)
//{
//  const List<int> & lNeighborTris = mTopo.TopVers()[pNearest].mT;
//  BrowseList(iT, lNeighborTris)
//  {
//      const LzTriModel::TopTri & lTri = mTopo.TopTris()[lNeighborTris.GetAt(iT)];
//
//      const Point3D & lA = mpPoints->at(lTri.mV[0]);
//      const Point3D & lB = mpPoints->at(lTri.mV[1]);
//      const Point3D & lC = mpPoints->at(lTri.mV[2]);
//
//      // Candidate Triangles
//      glColor3d(1, 0, 0);
//      //glLineWidth(5);
//      glBegin(GL_LINES);
//      glVertex3dv(lA.mV);
//      glVertex3dv(lB.mV);
//      glEnd();
//      glBegin(GL_LINES);
//      glVertex3dv(lA.mV);
//      glVertex3dv(lC.mV);
//      glEnd();
//      glBegin(GL_LINES);
//      glVertex3dv(lB.mV);
//      glVertex3dv(lC.mV);
//      glEnd();
//  }
//}


#pragma region "Node"
#pragma region "Construction / Destruction"
//================================================================================
Tree3D::Node::Node()
    : mpNodeLeft( nullptr ), mpNodeRight( nullptr ), mSeparator( 0 ), mIdxStart( 0 ), mIdxEnd( 0 ), mSize( 0 )
{
}

//================================================================================
Tree3D::Node::~Node()
{
    delete mpNodeLeft;
    delete mpNodeRight;
}
#pragma endregion


#pragma region "Computations"
//================================================================================
void Tree3D::Node::UpdateTree( Vector<Point3D> & pPoints,
							   Vector<unsigned int> & pIndices,
                               unsigned int & pDepthMax, unsigned int & pNbLeaf, unsigned int & pSizeMax, unsigned int & pSizeMin,
                               int pBegin, int pEnd,
							   unsigned int pCurrDepth,
                               unsigned int pMaxLeafCount, double pMinDiameter, unsigned int pMaxDepth )
{
    // Check
    if( pEnd < pBegin )
        LzLogException("", "Illegal indices! Begin= "<<pBegin<<", End= "<<pEnd<<".");

    // Begin <= End
    mIdxStart = pBegin;
    mIdxEnd = pEnd;

    // Compute bbox
    mBBox.Reset();
    for( int p=mIdxStart ; p<=mIdxEnd ; p++ )
        mBBox.Update( pPoints[p] );

    // Choose the best axis for the separation
    mAxis = BestAxis( pMinDiameter ); // returns -1 if bbox is smaller than pMinDiameter

    // If we cannot split the data, we get one leave
    if( mAxis<0
	 || (mIdxEnd - mIdxStart)<=(int)pMaxLeafCount
	 || pCurrDepth==pMaxDepth )
    {
        mSize = mIdxEnd - mIdxStart + 1;

        // Update stats
        if( pDepthMax < pCurrDepth ) pDepthMax = pCurrDepth;
        pNbLeaf++;
        if( pSizeMax < mSize ) pSizeMax = mSize;
        if( pSizeMin > mSize ) pSizeMin = mSize;
    }
    else
    {
        // Quick sort according to the best axis
        LzMath::ToolBox::QuickSort( pPoints.Buffer(), pBegin, pEnd, mAxis, pIndices.Buffer() );

        // Find the pivot
        int lMaxLeftIdx;
        mSeparator = Pivot( pPoints, lMaxLeftIdx );

// Log
//if( pCurrDepth < 3 )
//{
//  LzLogNode("", "Split");
//  LzLogM("", "Axis= "+mAxis+".");
//  LzLogM("", "Separator= "+mSeparator+".");
//  LzLogM("", "Left count = +"+(lMaxLeftIdx-pBegin+1)+".");
//  LzLogM("", "Right count= +"+(pEnd-lMaxLeftIdx)+".");
//}


        // Then go down the tree recursively
        mpNodeLeft = new Node;
        mpNodeLeft->UpdateTree( pPoints, pIndices,
                                pDepthMax, pNbLeaf, pSizeMax, pSizeMin,
                                pBegin, lMaxLeftIdx, pCurrDepth+1,
                                pMaxLeafCount, pMinDiameter, pMaxDepth );

        mpNodeRight = new Node;
        mpNodeRight->UpdateTree( pPoints, pIndices,
                                 pDepthMax, pNbLeaf, pSizeMax, pSizeMin,
                                 lMaxLeftIdx+1, pEnd, pCurrDepth+1,
                                 pMaxLeafCount, pMinDiameter, pMaxDepth );
    }
}

//================================================================================
const Tree3D::Node * Tree3D::Node::GetNearestLeaf( const Point3D & pPt ) const
{
    if( mSize )
        return this;
    else
    {
        if( pPt.mV[mAxis] <= mSeparator )
            return mpNodeLeft->GetNearestLeaf( pPt );
        else
            return mpNodeRight->GetNearestLeaf( pPt );
    }
}

//================================================================================
unsigned int Tree3D::Node::GetNearestPoint( const Point3D & pPt, const Vector<Point3D> & pPoints ) const
{
    double lMinDist = pPt.DistanceTo( pPoints[mIdxStart] );
    unsigned int lBestIndex = mIdxStart;

    for( int p = mIdxStart + 1 ; p <= mIdxEnd ; p++ )
    {
        double lDist = pPt.DistanceTo( pPoints[p] );

        if( lMinDist > lDist )
        {
            lMinDist = lDist;
            lBestIndex = p;
        }
    }

    return lBestIndex;
}

//================================================================================
void Tree3D::Node::GetNextLeaf( unsigned int & pConcurrentLeaves, unsigned int & pNearest, double & pMinDist, const Point3D & pPt, const Vector<Point3D> & pPoints ) const
{
    if( mSize )
    {
//**** pour ne pas rechercher dans la feuille trouvee en mode fast... pas tres elegant
		if( !mBBox.PointIsIn(pPt, false) ) 
//**** pour ne pas rechercher dans la feuille trouvee en mode fast... pas tres elegant
		{
            for( int iPoints = mIdxStart ; iPoints <= mIdxEnd ; iPoints++ )
            {
                double lDist = pPt.DistanceTo( pPoints[iPoints] );
                if( lDist < pMinDist )
                {
                    pMinDist = lDist;
                    pNearest = iPoints;
                }
            }

            // Update stats
            pConcurrentLeaves++;
        }
    }
    else
    {
        if( mpNodeLeft->IsChoice( pPt, pMinDist ) )
            mpNodeLeft->GetNextLeaf( pConcurrentLeaves, pNearest, pMinDist, pPt, pPoints );

        if( mpNodeRight->IsChoice( pPt, pMinDist )  )
            mpNodeRight->GetNextLeaf( pConcurrentLeaves, pNearest, pMinDist, pPt, pPoints );
    }
}

//================================================================================
void Tree3D::Node::CollectRangeVertices( const Point3D & pPt, double pRange, List<size_t> & pRangeIdx, const Vector<Point3D> & pPoints, const Vector<unsigned int> & pNew2Old ) const
{
    if( mSize )
    {
//**** optimisation by checking intersection with this node's bbox
		if( IsChoice( pPt, pRange ) )
//**** optimisation by checking intersection with this node's bbox
		{
			for( int p=mIdxStart ; p<=mIdxEnd ; p++ )
			{
				if( pPt.DistanceTo(pPoints[p]) < pRange )
					pRangeIdx.AddTail( pNew2Old[p] );
			}
		}
    }
    else
    {
        if( mpNodeLeft->IsChoice( pPt, pRange ) )
            mpNodeLeft->CollectRangeVertices( pPt, pRange, pRangeIdx, pPoints, pNew2Old );

        if( mpNodeRight->IsChoice( pPt, pRange )  )
            mpNodeRight->CollectRangeVertices( pPt, pRange, pRangeIdx, pPoints, pNew2Old );
    }
}

#pragma endregion


/*
//================================================================================
void Tree3D::Node::LogDisplay(const unsigned int pDepth, const vector<Point3D> * ppPoints) const
{
    if (mpNodeRight)
    {
        //DISPLAY RIGHT DATA
        if ((!mpNodeRight->mpNodeLeft && !mpNodeRight->mpNodeRight) || mpNodeRight->mAxis<0)
        {
            LzLogMsg("", "Depth : " + (pDepth + 1));
            for (int iData = mpNodeRight->mIdxStart; iData <= mpNodeRight->mIdxEnd; iData++)
                LzLogMsg("", " (" + ppPoints->at(iData).mV[0] + "," + ppPoints->at(iData).mV[1] + "," + ppPoints->at(iData).mV[2] + "),");
        }

        // OR GO TO THE RIGHT SUBTREE
        else
            mpNodeRight->LogDisplay(pDepth + 1, ppPoints);
    }

    //DISPLAY SEPARATOR
    LzLogMsg("", "Depth : " + pDepth);
    LzLogMsg("", mSeparator + ": axis " + mAxis);

    if (mpNodeLeft)
    {
        //DISPLAY LEFT DATA
        if ((!mpNodeLeft->mpNodeLeft && !mpNodeLeft->mpNodeRight) || mpNodeLeft->mAxis<0)
        {
            LzLogMsg("", "Depth : " + (pDepth + 1));
            for (int iData = mpNodeLeft->mIdxStart; iData <= mpNodeLeft->mIdxEnd; iData++)
                LzLogMsg("", " (" + ppPoints->at(iData).mV[0] + "," + ppPoints->at(iData).mV[1] + "," + ppPoints->at(iData).mV[2] + "),");
        }
        //OR GO TO THE LEFT SUBTREE
        else
            mpNodeLeft->LogDisplay(pDepth + 1, ppPoints);
    }
}

//================================================================================
void Tree3D::Node::Draw() const
{
    DrawHyperplane();
    if (!mpNodeLeft->mSize && mpNodeLeft->mAxis >= 0)
        mpNodeLeft->Draw();
    if (!mpNodeRight->mSize && mpNodeRight->mAxis >= 0)
        mpNodeRight->Draw();
}
*/

//================================================================================
int Tree3D::Node::BestAxis( double pMinDiameter ) const
{
    if( mBBox.SizeX() > pMinDiameter
	 || mBBox.SizeY() > pMinDiameter
	 || mBBox.SizeZ() > pMinDiameter )
        return mBBox.WidestDim();
	else
		return -1;
}

//================================================================================
double Tree3D::Node::Pivot( const Vector<Point3D> & pPoints, int & pMaxLeftIdx ) const
{
#if 0
    //
    // Dichotomic and possibly unbalanced tree
    //
    double lPivot = ( pPoints[mIdxStart].mV[mAxis] + pPoints[mIdxEnd].mV[mAxis] ) / 2.0;

    for( pMaxLeftIdx = mIdxStart ; pMaxLeftIdx <= mIdxEnd ; pMaxLeftIdx++ )
    {
        if( pPoints[pMaxLeftIdx].mV[mAxis] > lPivot )
            break;
    }

    pMaxLeftIdx--;

    return lPivot;
#else
    const int lSize = mIdxEnd - mIdxStart + 1;
    int lIdxLeft = mIdxStart + lSize / 2 - 1;
    int lIdxRight = mIdxStart + lSize / 2;

    double lPivot = ( pPoints[lIdxLeft].mV[mAxis] + pPoints[lIdxRight].mV[mAxis] ) / 2.0;
    pMaxLeftIdx = lIdxRight - 1;

    // If pivot is irrelevant (don't split the data), look around where should be the best pivot
    while( lIdxLeft  >= mIdxStart
            && lIdxRight <= mIdxEnd
            && pPoints[lIdxRight].mV[mAxis] - pPoints[lIdxLeft].mV[mAxis] < 1e-12 )
    {
        lIdxRight++;

        if( !LzMath::ToolBox::IsZero( pPoints[lIdxRight].mV[mAxis] - lPivot, 1e-12 ) )
        {
            lPivot = ( pPoints[lIdxRight].mV[mAxis] + pPoints[lIdxRight - 1].mV[mAxis] ) / 2.0;
            pMaxLeftIdx = lIdxRight - 1;
            break;
        }

        lIdxLeft--;

        if( !LzMath::ToolBox::IsZero( pPoints[lIdxLeft].mV[mAxis] - lPivot, 1e-12 ) )
        {
            lPivot = ( pPoints[lIdxLeft].mV[mAxis] + pPoints[lIdxLeft + 1].mV[mAxis] ) / 2.0;
            pMaxLeftIdx = lIdxLeft;
        }
    }

    // Check indices
    if( lIdxLeft < mIdxStart || lIdxRight > mIdxEnd )
        LzLogException("",  "Could not find pivot! Left/right indices are out of range (left= " << lIdxLeft << ", right= " << lIdxRight << ", begin=" << mIdxStart << ", end= " << mIdxEnd << ")" );

    return lPivot;
#endif
}

/*
//================================================================================
void Tree3D::Node::DrawHyperplane() const
{
    glEnable(GL_LIGHTING);
    switch (mAxis)
    {
    case 0:
        glColor3d(1, 0, 0);
        LzGeom::Plane3D(1, 0, 0, -mSeparator).Draw(Point3D(mSeparator, mBBox.Center().Y(), mBBox.Center().Z()), mBBox.SizeZ(), mBBox.SizeY());
        break;
    case 1:
        glColor3d(0, 1, 0);
        LzGeom::Plane3D(0, 1, 0, -mSeparator).Draw(Point3D(mBBox.Center().X(), mSeparator, mBBox.Center().Z()), mBBox.SizeZ(), mBBox.SizeX());
        break;
    case 2:
        glColor3d(0, 0, 1);
        LzGeom::Plane3D(0, 0, 1, -mSeparator).Draw(Point3D(mBBox.Center().X(), mBBox.Center().Y(), mSeparator), mBBox.SizeX(), mBBox.SizeY());
        break;
    default:
        LzLogMsg("", "ERROR : Hyperplane is crashing");
    }
}
*/

//================================================================================
//bool Tree3D::Node::IsChoice1( const Point3D & pPt, double pMinDist ) const
//{
//    if ( pPt.X() > mBBox.MaxX() )
//    {
//        const double lDeltaX = pPt.X() - mBBox.MaxX();
//
//        if ( pPt.Y() > mBBox.MaxY() )
//        {
//            const double lDeltaX2 = lDeltaX * lDeltaX;
//            const double lDeltaY = pPt.Y() - mBBox.MaxY();
//            const double lDeltaY2 = lDeltaY * lDeltaY;
//            const double lMinD2 = pMinDist * pMinDist;
//
//            // X>XM , Y>YM, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZM * lDeltaZM ) );
//            }
//            // X>XM , Y>YM, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z() ;
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZm * lDeltaZm ) );
//            }
//            // X>XM , Y>YM, Zm<=Z<=ZM
//            else
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 ) );
//        }
//        else if ( pPt.Y() < mBBox.MinY() )
//        {
//            const double lDeltaX2 = lDeltaX * lDeltaX;
//            const double lDeltaY = mBBox.MinY() - pPt.Y();
//            const double lDeltaY2 = lDeltaY * lDeltaY;
//            const double lMinD2 = pMinDist * pMinDist;
//
//            // X>XM , Y<Ym, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZM * lDeltaZM ) );
//            }
//            // X>XM , Y<Ym, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZm * lDeltaZm ) );
//            }
//            // X>XM , Y<Ym, Zm<=Z<=ZM
//            else
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 ) );
//        }
//        else
//        {
//            // X>XM , Ym<=Y<=YM, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaX * lDeltaX + lDeltaZM * lDeltaZM ) );
//            }
//            // X>XM , Ym<=Y<=YM, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaX * lDeltaX + lDeltaZm * lDeltaZm ) );
//            }
//            // X>XM , Ym<=Y<=YM, Zm<=Z<=ZM
//            else
//                return ( pMinDist >= lDeltaX );
//        }
//    }
//    else if ( pPt.X() < mBBox.MinX() )
//    {
//        const double lDeltaX = mBBox.MinX() - pPt.X();
//
//        if ( pPt.Y() > mBBox.MaxY() )
//        {
//            const double lMinD2 = pMinDist * pMinDist;
//            const double lDeltaX2 = lDeltaX * lDeltaX;
//            const double lDeltaY = pPt.Y() - mBBox.MaxY();
//            const double lDeltaY2 = lDeltaY * lDeltaY;
//
//            // X<Xm , Y>YM, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZM * lDeltaZM ) );
//            }
//            // X<Xm , Y>YM, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZm * lDeltaZm ) );
//            }
//            // X<Xm , Y>YM, Zm<=Z<=ZM
//            else
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 ) );
//        }
//        else if ( pPt.Y() < mBBox.MinY() )
//        {
//            const double lMinD2 = pMinDist * pMinDist;
//            const double lDeltaX2 = lDeltaX * lDeltaX;;
//            const double lDeltaY = mBBox.MinY() - pPt.Y();
//            const double lDeltaY2 = lDeltaY * lDeltaY;
//
//            // X<Xm , Y<Ym, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZM * lDeltaZM ) );
//            }
//            // X<Xm , Y<Ym, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 + lDeltaZm * lDeltaZm ) );
//            }
//            // X<Xm , Y<Ym, Zm<=Z<=ZM
//            else
//                return ( lMinD2 >= ( lDeltaX2 + lDeltaY2 ) );
//        }
//        else
//        {
//            // X<Xm , Ym<=Y<=YM, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaX * lDeltaX + lDeltaZM * lDeltaZM ) );
//            }
//            // X<Xm , Ym<=Y<=YM, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaX * lDeltaX + lDeltaZm * lDeltaZm ) );
//            }
//            // X<Xm , Ym<=Y<=YM, Zm<=Z<=ZM
//            else
//                return ( pMinDist >= lDeltaX );
//        }
//    }
//    else
//    {
//        if ( pPt.Y() > mBBox.MaxY() )
//        {
//            const double lDeltaY = pPt.Y() - mBBox.MaxY();
//
//            // Xm<=X<=XM , Y>YM, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaY * lDeltaY + lDeltaZM * lDeltaZM ) );
//            }
//            // Xm<=X<=XM , Y>YM, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaY * lDeltaY + lDeltaZm * lDeltaZm ) );
//            }
//            // Xm<=X<=XM , Y>YM, Zm<=Z<=ZM
//            else
//                return ( pMinDist >= lDeltaY );
//        }
//        else if ( pPt.Y() < mBBox.MinY() )
//        {
//            const double lDeltaY = mBBox.MinY() - pPt.Y();
//
//            // Xm<=X<=XM , Y<Ym, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//            {
//                const double lDeltaZM = pPt.Z() - mBBox.MaxZ();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaY * lDeltaY + lDeltaZM * lDeltaZM ) );
//            }
//            // Xm<=X<=XM , Y<Ym, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//            {
//                const double lDeltaZm = mBBox.MinZ() - pPt.Z();
//                return ( ( pMinDist * pMinDist ) >= ( lDeltaY * lDeltaY + lDeltaZm * lDeltaZm ) );
//            }
//            // Xm<=X<=XM , Y<Ym, Zm<=Z<=ZM
//            else
//                return pMinDist >= lDeltaY;
//        }
//        else
//        {
//            // Xm<=X<=XM, Ym<=Y<=YM, Z>ZM
//            if ( pPt.Z() > mBBox.MaxZ() )
//                return ( pMinDist >= ( pPt.Z() - mBBox.MaxZ() ) );
//            // Xm<=X<=XM , Ym<=Y<=YM, Z<Zm
//            else if ( pPt.Z() < mBBox.MinZ() )
//                return ( pMinDist >= ( mBBox.MinZ() - pPt.Z() ) );
//            // Xm<=X<=XM , Ym<=Y<=YM, Zm<Z<ZM
//            else
//                return true;
//        }
//    }
//}

//================================================================================
bool Tree3D::Node::IsChoice( const Point3D & pPt, double pMinDist ) const
{
    const double lMinDist2 = pMinDist * pMinDist;
	double lBBoxDist2 = 0.0;

	for( unsigned int d=0 ; d<3 ; d++ )
	{
		if( pPt.mV[d] > mBBox.Max(d) )
		{
			double lDelta = pPt.mV[d] - mBBox.Max(d);
			lBBoxDist2 += lDelta * lDelta;
		}
		else
		if( pPt.mV[d] < mBBox.Min(d) )
		{
			double lDelta = mBBox.Min(d) - pPt.mV[d];
			lBBoxDist2 += lDelta * lDelta;
		}
	}

	return lMinDist2 >= lBBoxDist2;
}


#pragma region "Display"
#if defined(USING_QT) && !defined(NO_OPENGL)
//================================================================================
void Tree3D::Node::DrawLeaves( double pMoyLeaf, unsigned int pSizeMax, unsigned int pSizeMin ) const
{
    if( !mSize )
    {
        mpNodeLeft->DrawLeaves( pMoyLeaf, pSizeMax, pSizeMin );
        mpNodeRight->DrawLeaves( pMoyLeaf, pSizeMax, pSizeMin );
    }
    else
    {
        const int lGap = ( int )abs( mSize - pMoyLeaf );

//******************************* a calculer en externe
        int lGapMax = ( ( ( double )pSizeMax - pMoyLeaf ) > ( pMoyLeaf - ( double )pSizeMin ) ) ? pSizeMax - ( int )pMoyLeaf : ( int )pMoyLeaf - pSizeMin;
        if( !lGapMax )
            lGapMax = 1;
//******************************* a calculer en externe

        LzMath::ToolBox::SetScaleColor( -lGapMax, lGapMax, lGap );
        mBBox.Draw();
    }
}

//================================================================================
void Tree3D::Node::DrawPlanes() const
{
    if( !mSize )
    {
        switch ( mAxis )
        {
        case 0:
            glColor3d( 1, 0, 0 );
            LzGeom::Plane3D( 1, 0, 0, -mSeparator ).Draw( Point3D( mSeparator, mBBox.Center().Y(), mBBox.Center().Z() ), mBBox.SizeZ(), mBBox.SizeY() );
            break;
        case 1:
            glColor3d( 0, 1, 0 );
            LzGeom::Plane3D( 0, 1, 0, -mSeparator ).Draw( Point3D( mBBox.Center().X(), mSeparator, mBBox.Center().Z() ), mBBox.SizeZ(), mBBox.SizeX() );
            break;
        case 2:
            glColor3d( 0, 0, 1 );
            LzGeom::Plane3D( 0, 0, 1, -mSeparator ).Draw( Point3D( mBBox.Center().X(), mBBox.Center().Y(), mSeparator ), mBBox.SizeX(), mBBox.SizeY() );
            break;
        default:
            LzLogMsg("", "ERROR: Hyperplane is crashing")
        }


        mpNodeLeft->DrawPlanes();
        mpNodeRight->DrawPlanes();
    }
}

//================================================================================
void Tree3D::Node::DrawTree2D( const Point3D & pPos, double pCurrW, double pOffY ) const
{
    if( !mSize )
    {
        const Point3D lPosL = pPos + Vector3D( -pCurrW / 2.0, pOffY, 0 );
        const Point3D lPosR = pPos + Vector3D( +pCurrW / 2.0, pOffY, 0 );

        // Root
        if( mAxis == 0 ) glColor3d( 1, 0, 0 );
        else if( mAxis == 1 ) glColor3d( 0, 1, 0 );
        else if( mAxis == 2 ) glColor3d( 0, 0, 1 );
        else
            LzLogException("",  "Unexpected axis! Axis= " << (int)mAxis << "." );

        glBegin( GL_POINTS );
        glVertex3dv( pPos.mV );
        glEnd();

        // Branches
        glColor3d( 0, 0.5, 0.5 );
        glBegin( GL_LINE_STRIP );
        glVertex2dv( lPosL.mV );
        glVertex2dv( pPos.mV );
        glVertex2dv( lPosR.mV );
        glEnd();

        // Recurse
        mpNodeLeft->DrawTree2D(  lPosL, pCurrW / 2.0, pOffY );
        mpNodeRight->DrawTree2D( lPosR, pCurrW / 2.0, pOffY );
    }
    else
    {
//glEnable( GL_LIGHTING );
//  DrawSphere( pPos, 5 );
//glDisable( GL_LIGHTING );

        glColor3d( 0, 1, 1 );
        glBegin( GL_POINTS );
        glVertex2dv( pPos.mV );
        glEnd();
    }
}
#endif
#pragma endregion
   

#pragma endregion
}


