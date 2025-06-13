#include "Line3D.h"
#include "Plane3D.h"
#include <LzMath/ToolBox.h>
#if defined(USING_QT) && !defined(NO_OPENGL)
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/gl.h>
#endif

namespace LzGeom
{
#pragma region "Construction / destruction"
//================================================================================
Line3D::Line3D()
: mRoot( 0, 0, 0 ), mDir( 1, 0, 0 )
{
    // Already normalized
}

//================================================================================
Line3D::Line3D( const Point3D & pRoot, const Vector3D & pDir )
: mRoot( pRoot ), mDir( pDir )
{
    Normalize();
}

//================================================================================
Line3D::Line3D( const Point3D & pA, const Point3D & pB )
: mRoot( pA ), mDir( pB - pA )
{
    Normalize();
}

//================================================================================
Line3D::Line3D( const Vector<Point3D> & pPointCloud )
{
	FromPointCloud( pPointCloud );
}

//================================================================================
Line3D::Line3D( const List<Point3D> & pPointCloud )
{
    Vector<Point3D> lCloud;
    pPointCloud.ToVector( lCloud );

    FromPointCloud( lCloud );
}

//================================================================================
void Line3D::FromPointCloud( const Vector<Point3D> & pPointCloud )
{
    int m = pPointCloud.Size();
    if( m < 2 )
        LzLogException("", "Unable to compute least-squares line! Only "<<m<<" point(s) provided.")

    Point3D lMean;
    Vector3D lV[3];
    LzMath::ToolBox::PCA( pPointCloud, lMean, lV );

	mRoot = lMean;
	mDir = lV[0];

    Normalize();
}

//================================================================================
void Line3D::Normalize()
{
    mDir.Normalize();
}
#pragma endregion


#pragma region "Computing"
//================================================================================
Line3D Line3D::operator+( const Vector3D & pV ) const
{
    return Line3D( mRoot + pV, mDir );
}

//================================================================================
Line3D Line3D::operator-( const Vector3D & pV ) const
{
    return Line3D( mRoot - pV, mDir );
}

//================================================================================
const Line3D & Line3D::operator+=( const Vector3D & pV )
{
    mRoot += pV;
    return *this;
}

//================================================================================
const Line3D & Line3D::operator-=( const Vector3D & pV )
{
    mRoot -= pV;
    return *this;
}

//================================================================================
bool Line3D::IsParallelTo( const Line3D & pLine ) const
{
    Vector3D lCross = mDir ^ pLine.mDir;
    return LzMath::ToolBox::IsZero( lCross.Norm() );
}

//================================================================================
double Line3D::DistanceTo( const Point3D & pPoint ) const
{
    return pPoint.DistanceTo( Projection( pPoint ) );
}

//================================================================================
double Line3D::DistanceTo( const Line3D & pLine, bool pMinMode/*=false*/ ) const
{
if( pMinMode )
{
    // If
    //
    // a = this->mRoot
    // u = this->mDir
    // A(t) = a + t * u
    //
    // b = pLine.mRoot
    // v = pLine.mDir
    // B(s) = b - s * v (note the '-' instead of '+')
    //
    // (s, t) --> min of: || A(t) - B(s) ||^2
    //
    // D = (a - b), Du = D * u, Dv = D * v,
    // u2 = u * u, v2 = v * v, uv = u * v
    //
    // then
    //
    // | u2  uv | | t |   | -Dv |
    // | uv  v2 | | s | = | -Du |
    //
    const double u2 = mDir.SquaredNorm();
    const double v2 = pLine.mDir.SquaredNorm();
    const double uv = mDir * pLine.mDir;

    // Compute determinant
    const double lDet = u2*v2 - uv*uv;

    // Check parallelism using the determinant
    if( LzMath::ToolBox::IsZero(lDet) )
    {
        // Lines are parallel - but distance can be computed nonetheless

        // Project my root on other line
        Point3D lProjRoot = pLine.Projection( mRoot );
        return lProjRoot.DistanceTo( mRoot );
    }
    else
    {
        // Compute missing quantities
        const Vector3D D = mRoot - pLine.mRoot;
        const double Du = D * mDir;
        const double Dv = D * pLine.mDir;

        // Inverse matrix is
        //           |  v2  -uv |
        // (1 / Det) | -uv   u2 |
        //
        const double t = ( -v2*Du +uv*Dv ) / lDet;
        const double s = ( +uv*Du -u2*Dv ) / lDet;

        // Find nearest points
        const Point3D A =       mRoot + t *       mDir;
        const Point3D B = pLine.mRoot - s * pLine.mDir;

        // Compute distance
        return A.DistanceTo(B);
    }
}
else
{
    // Compute distance using intersections
    if( IsParallelTo( pLine ) )
    {
        // Lines are parallel - but distance can be computed nonetheless

        // Project my root on other line
        Point3D lProjRoot = pLine.Projection( mRoot );
        return lProjRoot.DistanceTo( mRoot );
    }
    else
    {
        Point3D lA, lB;
        PseudoIntersection( pLine, lA, lB );
        return lA.DistanceTo( lB );
    }
}
}

//================================================================================
Point3D Line3D::Projection( const Point3D & pPoint ) const
{
    double d = ( pPoint - mRoot ) * mDir;
    return mRoot + d * mDir;
}

//================================================================================
Point3D Line3D::Intersection( const Plane3D & pPlane ) const
{
    return pPlane.Intersection( *this );
}

//================================================================================
bool Line3D::Intersection_NoExc( const Plane3D & pPlane, Point3D & pTo ) const
{
	return pPlane.Intersection_NoExc( *this, pTo );
}

//================================================================================
Point3D Line3D::PseudoIntersection( const Line3D & pLine ) const
{
    Point3D lA, lB;
    PseudoIntersection( pLine, lA, lB );

    return lA + 0.5 * ( lB - lA );
}

//================================================================================
bool Line3D::PseudoIntersection_NoExc( const Line3D & pLine, Point3D & pTo ) const
{
    Point3D lA, lB;
    if( PseudoIntersection_NoExc(pLine, lA, lB) )
	{
        pTo = lA.MidPointTo( lB );
		return true;
	}
	else
		return false;
}

//================================================================================
void Line3D::PseudoIntersection( const Line3D & pOtherLine, Point3D & pMyPoint, Point3D & pOtherPoint ) const
{
    // Compute transverse vector
    Vector3D lTrans = mDir ^ pOtherLine.mDir;

    // My plane
    Plane3D lMyPlane( mRoot, mDir, lTrans );
//**** BUG interversion des points! 2019-06-18
//pMyPoint = lMyPlane.Intersection( pOtherLine );
pOtherPoint = lMyPlane.Intersection( pOtherLine );
//**** BUG interversion des points! 2019-06-18

    // Other plane
    Plane3D lOtherPlane( pOtherLine.mRoot, pOtherLine.mDir, lTrans );
//**** BUG interversion des points! 2019-06-18
//pOtherPoint = lOtherPlane.Intersection( *this );
pMyPoint = lOtherPlane.Intersection( *this );
//**** BUG interversion des points! 2019-06-18
}

//================================================================================
bool Line3D::PseudoIntersection_NoExc( const Line3D & pOtherLine, Point3D & pMyPoint, Point3D & pOtherPoint ) const
{
    // Compute transverse vector
    Vector3D lTrans = mDir ^ pOtherLine.mDir;

    // My plane
	{
		Vector3D lNormal = mDir ^ lTrans;

		// Check
        if( LzMath::ToolBox::IsZero( lNormal.Norm() ) )
			return false;

		Plane3D lMyPlane( mRoot, lNormal );

		// Check
//**** BUG interversion des points! 2019-06-18
//if( !lMyPlane.Intersection_NoExc(pOtherLine, pMyPoint) )
//	return false;
if( !lMyPlane.Intersection_NoExc(pOtherLine, pOtherPoint) )
	return false;
//**** BUG interversion des points! 2019-06-18
	}

    // Other plane
	{
		Vector3D lNormal = pOtherLine.mDir ^ lTrans;

		// Check
        if( LzMath::ToolBox::IsZero( lNormal.Norm() ) )
			return false;

		Plane3D lOtherPlane( pOtherLine.mRoot, lNormal );

		// Check
//**** BUG interversion des points! 2019-06-18
//if( !lOtherPlane.Intersection_NoExc(*this, pOtherPoint) )
//	return false;
if( !lOtherPlane.Intersection_NoExc(*this, pMyPoint) )
	return false;
//**** BUG interversion des points! 2019-06-18
	}

	return true;
}
#pragma endregion


#pragma region "Drawing"
#if !defined(NO_OPENGL)
//================================================================================
void Line3D::Draw( const Point3D & pCenter, double pLen ) const
{
    // // Build referential
    // Point3D C = Projection( pCenter );

    // double lHalfLen = pLen / 2;

    // glBegin( GL_LINES );
    // {
    //     glVertex3dv( ( C - lHalfLen * mDir ).mV );
    //     glVertex3dv( ( C + lHalfLen * mDir ).mV );
    // }
    // glEnd();
}
#endif
#pragma endregion
}
