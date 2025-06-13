#pragma once

#include "Point3D.h"
#include "Vector3D.h"
#include <LzServices/Vector.h>
#include <LzServices/List.h>


namespace LzGeom
{
class Plane3D;


using LzServices::Vector;
using LzServices::List;


class Line3D
{
#pragma region "Construction / destruction"
public:
    Line3D();
    Line3D( const Point3D & pRoot, const Vector3D & pDir );
    Line3D( const Point3D & pA, const Point3D & pB );
	Line3D( const Vector<Point3D> & pPointCloud );
	Line3D( const List<Point3D> & pPointCloud );
protected:
    void FromPointCloud( const Vector<Point3D> & pPointCloud );
    void Normalize();
#pragma endregion


#pragma region "Computing"
public:
    const Point3D & Root() const
    {
        return mRoot;
    }
	// Dir is always normalized
    const Vector3D & Dir() const
    {
        return mDir;
    }
    Line3D operator+( const Vector3D & pV ) const;
    Line3D operator-( const Vector3D & pV ) const;
    const Line3D & operator+=( const Vector3D & pV );
    const Line3D & operator-=( const Vector3D & pV );
    bool IsParallelTo( const Line3D & pLine ) const;
    double DistanceTo( const Point3D & pPoint ) const;
    double DistanceTo( const Line3D & pLine, bool pMinMode=false ) const;
    Point3D Projection( const Point3D & pPoint ) const;
	//
    Point3D Intersection( const Plane3D & pPlane ) const;
    bool Intersection_NoExc( const Plane3D & pPlane, Point3D & pTo ) const;
	//
    Point3D PseudoIntersection( const Line3D & pLine ) const;
	bool PseudoIntersection_NoExc( const Line3D & pLine, Point3D & pTo ) const;
	//
	void PseudoIntersection( const Line3D & pOtherLine, Point3D & pMyPoint, Point3D & pOtherPoint ) const;
	bool PseudoIntersection_NoExc( const Line3D & pOtherLine, Point3D & pMyPoint, Point3D & pOtherPoint ) const;
#pragma endregion


#pragma region "Drawing"
#if !defined(NO_OPENGL)
public:
    void Draw( const Point3D & pCenter, double pLen ) const;
#endif
#pragma endregion


#pragma region "Data"
protected:
    Point3D mRoot;
    Vector3D mDir;

    friend class Plane3D;
#pragma endregion
};
}
