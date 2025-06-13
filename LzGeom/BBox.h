#pragma once

#include "LzLib_Geom.h"
#include "Point3D.h"
#include "Plane3D.h"
#include <LzServices/Vector.h>
#include <vector>


// DLL exports
namespace LzGeom
{
    class DLL_EXPORT_LzGeom BBox;
}


namespace LzGeom
{
using LzServices::Vector;
using std::vector;
using std::string;


class BBox
{
public:
    BBox();
    BBox( const Point3D & pPt );
    BBox( const vector<Point3D> & pPts );
    BBox( const Vector<Point3D> & pPts );
    BBox( const List<Point3D> & pPts );
    BBox( double pMinX, double pMaxX, double pMinY, double pMaxY, double pMinZ, double pMaxZ );
    virtual ~BBox() {}

    void Reset();
    bool IsValid() const
    {
        return mMin[0]<=mMax[0] && mMin[1]<=mMax[1] && mMin[2]<=mMax[2];
    }

    void Update( const Point3D & pPt );
    void Update( const vector<Point3D> & pPts );
    void Update( const Vector<Point3D> & pPts );
    void Update( const List<Point3D> & pPts );
    void Update( const Point3D * ppPts, unsigned int pCount );

    Point3D Root() const { return Point3D( mMin[0], mMin[1], mMin[2] ); }
    Point3D Center() const { return Point3D(0.5*(mMin[0] + mMax[0]), 0.5*(mMin[1] + mMax[1]), 0.5*(mMin[2] + mMax[2])); }

    double Size( unsigned int pDim ) const;
    double SizeX() const { return Size( 0 ); }
    double SizeY() const { return Size( 1 ); }
    double SizeZ() const { return Size( 2 ); }

    double Volume() const { return SizeX() * SizeY() * SizeZ(); }

    unsigned int WidestDim() const;
    unsigned int ShortestDim() const;

//    double & Min( unsigned int pDim );
//    double & MinX() { return Min( 0 ); }
//    double & MinY() { return Min( 1 ); }
//    double & MinZ() { return Min( 2 ); }

    double Min( unsigned int pDim ) const;
    double MinX() const { return Min( 0 ); }
    double MinY() const { return Min( 1 ); }
    double MinZ() const { return Min( 2 ); }

//    double & Max( unsigned int pDim );
//    double & MaxX() { return Max( 0 ); }
//    double & MaxY() { return Max( 1 ); }
//    double & MaxZ() { return Max( 2 ); }

    double Max( unsigned int pDim ) const;
    double MaxX() const { return Max( 0 ); }
    double MaxY() const { return Max( 1 ); }
    double MaxZ() const { return Max( 2 ); }

    void Scale( double pS );
    void Scale( double pS0, double pS1, double pS2 );

	// Expands BBox by a certain factor. E.g. pFactMarginX = 10% = 0.1 ==> total added margin in X = 0.1 * old SizeX
    void ExpandBy( double pFactMarginX, double pFactMarginY, double pFactMarginZ );
    void ExpandBy( double pFactAllMargins ) { ExpandBy( pFactAllMargins, pFactAllMargins, pFactAllMargins ); }
	//
	// Expands BBox by a certain absolute margin. E.g. pMarginX = 20.0 ==> new SizeX = 20.0 + old SizeX
    void AddMargin( double pMarginX, double pMarginY, double pMarginZ );
    void AddMargin( double pMargin ) { AddMargin( pMargin, pMargin, pMargin ); }

	// Rigid transformations and mirroring
	// /!\ BBox remains cartesian ==> careful with rotations/mirror planes!
	void RigidTransform( const RigidTr3D & pTr );
	const BBox & operator+=( const Vector3D & pV );
	const BBox & operator-=( const Vector3D & pV );
	void Mirror( const Plane3D & pMirror );

    Point3D RandomPoint() const;
    const Point3D & HexaPoint( unsigned int pI ) const;
    bool PointIsIn( const Point3D & pPt, bool pIncludingBoundary, double pBoundaryMargin=1e-10 ) const;

    void Log( const std::string & pTag = "BBox" ) const;
    string toString() const;

#if defined(USING_QT) && !defined(NO_OPENGL)
    // Rendering
    void Draw() const;
    void DrawFill() const;
#endif

    void ComputeHexaPoints(); // Must be called explicitly if bounds are changed through non-const accessors Min() or Max()

protected:
    double mMin[3];
    double mMax[3];

    Point3D mHexaPoints[8];
};
}
