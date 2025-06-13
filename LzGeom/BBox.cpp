#include "BBox.h"
#include <LzMath/ToolBox.h>
#include <LzServices/LzLog.h>
#include <limits>
#if defined(USING_QT) && !defined(NO_OPENGL)
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/gl.h>
#endif


//*************** conflicting with std::numeric_limits<double>::max()
#undef max
//*************** conflicting with std::numeric_limits<double>::max()


namespace LzGeom
{
//================================================================================
BBox::BBox()
{
    Reset();
}

//================================================================================
BBox::BBox( const Point3D & pPt )
{
    Reset();
    Update( pPt );
}

//================================================================================
BBox::BBox( const vector<Point3D> & pPts )
{
    Reset();
    Update( pPts );
}

//================================================================================
BBox::BBox( const Vector<Point3D> & pPts )
{
    Reset();
    Update( pPts );
}

//================================================================================
BBox::BBox( const List<Point3D> & pPts )
{
    Reset();
    Update( pPts );
}

//================================================================================
BBox::BBox( double pMinX, double pMaxX, double pMinY, double pMaxY, double pMinZ, double pMaxZ )
{
    mMin[0] = pMinX;
    mMax[0] = pMaxX;

    mMin[1] = pMinY;
    mMax[1] = pMaxY;

    mMin[2] = pMinZ;
    mMax[2] = pMaxZ;

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Reset()
{
    for( unsigned int i=0 ; i<3 ; i++ )
    {
        mMin[i] = +std::numeric_limits<double>::max();
        mMax[i] = -std::numeric_limits<double>::max();
    }

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Update( const Point3D & pPt )
{
    for( int i=0 ; i<3 ; i++ )
    {
        double lV = pPt.mV[i];
        if( mMin[i] > lV ) mMin[i] = lV;
        if( mMax[i] < lV ) mMax[i] = lV;
    }

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Update( const vector<Point3D> & pPts )
{
    for( vector<Point3D>::const_iterator iPt=pPts.begin() ; iPt!=pPts.end() ; iPt++ )
    {
        for( int i=0 ; i<3 ; i++ )
        {
            double lV = iPt->mV[i];
            if( mMin[i] > lV ) mMin[i] = lV;
            if( mMax[i] < lV ) mMax[i] = lV;
        }
    }

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Update( const Vector<Point3D> & pPts )
{
    for( size_t p=0 ; p<pPts.Size() ; p++ )
    {
		const Point3D & lPt = pPts[p];
        for( int i=0 ; i<3 ; i++ )
        {
            double lV = lPt.mV[i];
            if( mMin[i] > lV ) mMin[i] = lV;
            if( mMax[i] < lV ) mMax[i] = lV;
        }
    }

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Update( const List<Point3D> & pPts )
{
	BrowseList( iP, pPts )
	{
		const Point3D & lPt = pPts.GetAt( iP );
        for( int i=0 ; i<3 ; i++ )
        {
            double lV = lPt.mV[i];
            if( mMin[i] > lV ) mMin[i] = lV;
            if( mMax[i] < lV ) mMax[i] = lV;
        }
	}

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Update( const Point3D * ppPts, unsigned int pCount )
{
    for( unsigned int p=0 ; p<pCount ; p++ )
    {
        for( unsigned int i=0 ; i<3 ; i++ )
        {
            double lV = ppPts[p].mV[i];
            if( mMin[i] > lV ) mMin[i] = lV;
            if( mMax[i] < lV ) mMax[i] = lV;
        }
    }

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
double BBox::Size( unsigned int pDim ) const
{
    // Check
    if( pDim > 2 )
        LzLogException("", "Accessing invalid bbox dimension: " << pDim << "!")

    return mMax[pDim] - mMin[pDim];
}

//================================================================================
unsigned int BBox::WidestDim() const
{
    unsigned int lMaxDim = 0;
    double lMaxSize = Size( lMaxDim );

    for( unsigned int d = 1 ; d < 3 ; d++ )
    {
        double lSize = Size( d );

        // Wider?
        if( lMaxSize < lSize )
        {
            lMaxSize = lSize;
            lMaxDim = d;
        }
    }

    return lMaxDim;
}

//================================================================================
unsigned int BBox::ShortestDim() const
{
    unsigned int lMinDim = 0;
    double lMinSize = Size( lMinDim );

    for( unsigned int d = 1 ; d < 3 ; d++ )
    {
        double lSize = Size( d );

        // Shorter?
        if( lMinSize > lSize )
        {
            lMinSize = lSize;
            lMinDim = d;
        }
    }

    return lMinDim;
}

////================================================================================
//double & BBox::Min( unsigned int pDim )
//{
//	// Check
//    if( pDim > 2 )
//        LzLogException("", "Accessing invalid bbox dimension: " << pDim << "!")

//	return mMin[pDim];
//}

//================================================================================
double BBox::Min( unsigned int pDim ) const
{
	// Check
    if( pDim > 2 )
        LzLogException("", "Accessing invalid bbox dimension: " << pDim << "!")

	return mMin[pDim];
}

////================================================================================
//double & BBox::Max( unsigned int pDim )
//{
//	// Check
//	if( pDim > 2 )
//        LzLogException("", "Accessing invalid bbox dimension: " << pDim << "!")

//	return mMax[pDim];
//}

//================================================================================
double BBox::Max( unsigned int pDim ) const
{
	// Check
	if( pDim > 2 )
        LzLogException("", "Accessing invalid bbox dimension: " << pDim << "!")

	return mMax[pDim];
}

//================================================================================
void BBox::Scale( double pS )
{
    // Scale bounds
    for( int i = 0 ; i < 3 ; i++ )
    {
        mMin[i] *= pS;
        mMax[i] *= pS;
    }

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::Scale( double pS0, double pS1, double pS2 )
{
    // Scale bounds
    mMin[0] *= pS0;
    mMax[0] *= pS0;
    mMin[1] *= pS1;
    mMax[1] *= pS1;
    mMin[2] *= pS2;
    mMax[2] *= pS2;

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::ExpandBy( double pFactMarginX, double pFactMarginY, double pFactMarginZ )
{
#if 1
    double lMarginX = SizeX() * pFactMarginX;
    double lMarginY = SizeY() * pFactMarginY;
    double lMarginZ = SizeZ() * pFactMarginZ;

	AddMargin( lMarginX, lMarginY, lMarginZ );
#else
    // Add margin
    double lMarginX = SizeX() * pFactSizeX;
    mMin[0] -= lMarginX / 2;
    mMax[0] += lMarginX / 2;

    double lMarginY = SizeY() * pFactSizeY;
    mMin[1] -= lMarginY / 2;
    mMax[1] += lMarginY / 2;

    double lMarginZ = SizeZ() * pFactSizeZ;
    mMin[2] -= lMarginZ / 2;
    mMax[2] += lMarginZ / 2;

    // Update hexa points
    ComputeHexaPoints();
#endif
}

//================================================================================
void BBox::AddMargin( double pMarginX, double pMarginY, double pMarginZ )
{
    mMin[0] -= pMarginX / 2;
    mMax[0] += pMarginX / 2;

    mMin[1] -= pMarginY / 2;
    mMax[1] += pMarginY / 2;

    mMin[2] -= pMarginZ / 2;
    mMax[2] += pMarginZ / 2;

    // Update hexa points
    ComputeHexaPoints();
}

//================================================================================
void BBox::RigidTransform( const RigidTr3D & pTr )
{
	// Compute transformed hexa points
	Point3D lTrHexaPoints[8];
	for( int h=0 ; h<8 ; h++ )
		lTrHexaPoints[h] = pTr * mHexaPoints[h];

	// Recomute bbox
	Reset();
	Update( lTrHexaPoints, 8 );
}

//================================================================================
const BBox & BBox::operator+=( const Vector3D & pV )
{
	// Compute transformed hexa points
	Point3D lTrHexaPoints[8];
	for( int h=0 ; h<8 ; h++ )
		lTrHexaPoints[h] = mHexaPoints[h] + pV;

	// Recomute bbox
	Reset();
	Update( lTrHexaPoints, 8 );

    return *this;
}

//================================================================================
const BBox & BBox::operator-=( const Vector3D & pV )
{
	// Compute transformed hexa points
	Point3D lTrHexaPoints[8];
	for( int h=0 ; h<8 ; h++ )
		lTrHexaPoints[h] = mHexaPoints[h] - pV;

	// Recomute bbox
	Reset();
	Update( lTrHexaPoints, 8 );

    return *this;
}

//================================================================================
void BBox::Mirror( const Plane3D & pMirror )
{
	// Compute mirrored hexa points
	Point3D lMirrHexaPoints[8];
	for( int h=0 ; h<8 ; h++ )
		lMirrHexaPoints[h] = pMirror.Symmetrical( mHexaPoints[h] );

	// Recomute bbox
	Reset();
	Update( lMirrHexaPoints, 8 );
}

//================================================================================
Point3D BBox::RandomPoint() const
{
	// Create uniform distributions
	std::uniform_real_distribution<double> lDist_X( MinX(), MaxX() ); // Inclusive interval
	std::uniform_real_distribution<double> lDist_Y( MinY(), MaxY() ); // Inclusive interval
	std::uniform_real_distribution<double> lDist_Z( MinZ(), MaxZ() ); // Inclusive interval

    return Point3D( lDist_X(LzServices::RandomEngine()), lDist_Y(LzServices::RandomEngine()), lDist_Z(LzServices::RandomEngine()) );
}

//================================================================================
const Point3D & BBox::HexaPoint( unsigned int pI ) const
{
    // Check
    if( pI > 7 )
        LzLogException("",  "Accessing invalid hexa point: " << pI << "!" );

    return mHexaPoints[pI];
}

//================================================================================
bool BBox::PointIsIn( const Point3D & pPt, bool pIncludingBoundary, double pBoundaryMargin/*=1e-10*/ ) const
{
    if( pIncludingBoundary )
        return ( pPt.X()>mMin[0]-pBoundaryMargin && pPt.X()<mMax[0]+pBoundaryMargin )
            && ( pPt.Y()>mMin[1]-pBoundaryMargin && pPt.Y()<mMax[1]+pBoundaryMargin )
            && ( pPt.Z()>mMin[2]-pBoundaryMargin && pPt.Z()<mMax[2]+pBoundaryMargin );
    else
        return ( pPt.X()>mMin[0] && pPt.X()<mMax[0] )
            && ( pPt.Y()>mMin[1] && pPt.Y()<mMax[1] )
            && ( pPt.Z()>mMin[2] && pPt.Z()<mMax[2] );
}

//================================================================================
void BBox::Log( const std::string & pTag/*="BBox"*/ ) const
{
    LzLogNode("", pTag)
    LzLogM("", "X: " << mMin[0] << " to " << mMax[0] << ", size X= " << Size( 0 ))
    LzLogM("", "Y: " << mMin[1] << " to " << mMax[1] << ", size Y= " << Size( 1 ))
    LzLogM("", "Z: " << mMin[2] << " to " << mMax[2] << ", size Z= " << Size( 2 ))
    LzLogM("", "")
    LzLogM("", "Center: " << 0.5 * ( mMin[0] + mMax[0] ) << ", " << 0.5 * ( mMin[1] + mMax[1] ) << ", " << 0.5 * ( mMin[2] + mMax[2] ) << ".")
}

//================================================================================
string BBox::toString() const
{
    return "["+std::to_string(mMin[0])+" to "+std::to_string(mMax[0])+"; "
              +std::to_string(mMin[1])+" to "+std::to_string(mMax[1])+"; "
              +std::to_string(mMin[2])+" to "+std::to_string(mMax[2])+"]";
}


#if defined(USING_QT) && !defined(NO_OPENGL)
//================================================================================
void BBox::Draw() const
{
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    glDisable( GL_LIGHTING );
    glDisable( GL_CULL_FACE );

    glBegin( GL_QUAD_STRIP );
    {
        glVertex3d( mMax[0], mMax[1], mMax[2] );
        glVertex3d( mMin[0], mMax[1], mMax[2] );

        glVertex3d( mMax[0], mMin[1], mMax[2] );
        glVertex3d( mMin[0], mMin[1], mMax[2] );

        glVertex3d( mMax[0], mMin[1], mMin[2] );
        glVertex3d( mMin[0], mMin[1], mMin[2] );

        glVertex3d( mMax[0], mMax[1], mMin[2] );
        glVertex3d( mMin[0], mMax[1], mMin[2] );

        glVertex3d( mMax[0], mMax[1], mMax[2] );
        glVertex3d( mMin[0], mMax[1], mMax[2] );
    }
    glEnd();

    glPopAttrib();
}

//================================================================================
void BBox::DrawFill() const
{
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    glBegin( GL_QUADS );
    {
        // Max X
        glNormal3d( +1, 0, 0 );
        glVertex3d( mMax[0], mMax[1], mMax[2] );
        glVertex3d( mMax[0], mMin[1], mMax[2] );
        glVertex3d( mMax[0], mMin[1], mMin[2] );
        glVertex3d( mMax[0], mMax[1], mMin[2] );

        // Min X
        glNormal3d( -1, 0, 0 );
        glVertex3d( mMin[0], mMax[1], mMax[2] );
        glVertex3d( mMin[0], mMax[1], mMin[2] );
        glVertex3d( mMin[0], mMin[1], mMin[2] );
        glVertex3d( mMin[0], mMin[1], mMax[2] );

        // Max Y
        glNormal3d( 0, +1, 0 );
        glVertex3d( mMax[0], mMax[1], mMax[2] );
        glVertex3d( mMax[0], mMax[1], mMin[2] );
        glVertex3d( mMin[0], mMax[1], mMin[2] );
        glVertex3d( mMin[0], mMax[1], mMax[2] );

        // Min Y
        glNormal3d( 0, -1, 0 );
        glVertex3d( mMax[0], mMin[1], mMax[2] );
        glVertex3d( mMin[0], mMin[1], mMax[2] );
        glVertex3d( mMin[0], mMin[1], mMin[2] );
        glVertex3d( mMax[0], mMin[1], mMin[2] );

        // Max Z
        glNormal3d( 0, 0, +1 );
        glVertex3d( mMax[0], mMax[1], mMax[2] );
        glVertex3d( mMin[0], mMax[1], mMax[2] );
        glVertex3d( mMin[0], mMin[1], mMax[2] );
        glVertex3d( mMax[0], mMin[1], mMax[2] );

        // Min Z
        glNormal3d( 0, 0, -1 );
        glVertex3d( mMax[0], mMax[1], mMin[2] );
        glVertex3d( mMax[0], mMin[1], mMin[2] );
        glVertex3d( mMin[0], mMin[1], mMin[2] );
        glVertex3d( mMin[0], mMax[1], mMin[2] );
    }
    glEnd();

    glPopAttrib();
}
#endif

//================================================================================
void BBox::ComputeHexaPoints()
{
    mHexaPoints[0] = Point3D( mMin[0], mMin[1], mMin[2] );
    mHexaPoints[1] = Point3D( mMax[0], mMin[1], mMin[2] );
    mHexaPoints[2] = Point3D( mMax[0], mMax[1], mMin[2] );
    mHexaPoints[3] = Point3D( mMin[0], mMax[1], mMin[2] );

    mHexaPoints[4] = Point3D( mMin[0], mMin[1], mMax[2] );
    mHexaPoints[5] = Point3D( mMax[0], mMin[1], mMax[2] );
    mHexaPoints[6] = Point3D( mMax[0], mMax[1], mMax[2] );
    mHexaPoints[7] = Point3D( mMin[0], mMax[1], mMax[2] );
}
}
