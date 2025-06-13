#pragma once

#include "Point3D.h"
#include "Vector3D.h"
#include <LzServices/List.h>
#include <LzServices/Vector.h>
#include <LzMath/Matrix.h>
#include <vector>


namespace LzGeom
{
using LzServices::List;
using LzServices::Vector;
using LzMath::Matrix;
using std::vector;

class Line3D;


class Plane3D
{
#pragma region "Construction / destruction"
public:
    Plane3D();
    Plane3D( std::function<void(Plane3D * pThis)> pInit );
    Plane3D( double pA, double pB, double pC, double pD );
    Plane3D( const double * pP );
    Plane3D( const Point3D & pPoint, const Vector3D & pNormal );
    Plane3D( const Point3D & pPoint, const Vector3D & pX, const Vector3D & pY );
    Plane3D( const Point3D & pPoint1, const Point3D & pPoint2, const Point3D & pPoint3 ); // CCW side is +
    Plane3D( const vector<Point3D> & pPointCloud );
    Plane3D( const Vector<Point3D> & pPointCloud );
    Plane3D( const List<Point3D> & pPointCloud );
protected:
    void FromPointCloud( const vector<Point3D> & pPointCloud );
    void Normalize();
#pragma endregion


#pragma region "Access"
public:
    double A() const { return mP[0]; }
    double B() const { return mP[1]; }
    double C() const { return mP[2]; }
    double D() const { return mP[3]; }
//double & P( unsigned int pP ) { return mP[pP]; } ******** NON car on ne maitrise plus la contrainte de Normalize()
//double P( unsigned int pP ) const { return mP[pP]; }
    const double * P() const { return mP; }
    const Vector3D & X() const { return mX; }
    const Vector3D & Y() const { return mY; }
    const Vector3D & Normal() const { return mNormal; }
#pragma endregion


#pragma region "Modifier"
public:
    void SetPlaneCoeffs( double pA, double pB, double pC, double pD );
#pragma endregion


#pragma region "Computing"
public:
    //-----------------------------------------------------------------------------------------
    // Always call Normalize() in methods that change the normal (i.e. mP[0], mP[1] or mP[2])
    //-----------------------------------------------------------------------------------------
    void FlipNormal();
    double SignedDistanceTo( const Point3D & pPoint ) const;
    Point3D Symmetrical( const Point3D & pPoint ) const;
    Vector3D Symmetrical( const Vector3D & pVector ) const;
    Point3D Projection( const Point3D & pPoint ) const;
    Vector3D Projection( const Vector3D & pV ) const;
	//
    Point3D Projection( const Point3D & pPoint, const Vector3D & pProjDir ) const;
    bool Projection_NoExc( const Point3D & pPoint, const Vector3D & pProjDir, Point3D & pTo ) const;
	//
    Line3D Intersection( const Plane3D & pOther ) const;
    bool Intersection_NoExc( const Plane3D & pOther, Line3D & pTo ) const;
	//
    Point3D Intersection( const Line3D & pLine ) const;
    bool Intersection_NoExc( const Line3D & pLine, Point3D & pTo ) const;
	//
    void SetPositiveSide( const Vector3D & pToPos );
    void SetPositiveSide( const Point3D & pInPos );
    void SetNegativeSide( const Point3D & pInNeg );
    const Plane3D & operator*=( const RigidTr3D & pT );
    const Plane3D & operator*=( double pScale );
    const Plane3D & operator/=( double pScale );
    Plane3D operator+( const Vector3D & pV ) const;
    Plane3D operator-( const Vector3D & pV ) const;
    const Plane3D & operator+=( const Vector3D & pV );
    const Plane3D & operator-=( const Vector3D & pV );
    bool operator!=( const Plane3D & pOther ) const;
    Matrix GetSymmetryMat4x4() const;
    //
    bool IsSameAs( const Plane3D & pOther, double pTol ) const;
#pragma endregion


#pragma region "Drawing"
public:
#if !defined(NO_OPENGL)
    void Draw( const Point3D & pCenter, double pDimX, double pDimY = -1 ) const;
    void Draw( const Point3D & pCenter, double pDimX, double pDimY, unsigned int pStripsX, unsigned int pStripsY ) const;
#endif
    void DrawToArray( const Point3D & pCenter, double pDimX, double pDimY, vector<Point3D> & pQuad ) const;
#pragma endregion


#pragma region "Debug"
public:
    std::string ToString() const;
    std::string ToString_ABCD() const;
#pragma endregion


#pragma region "Data"
protected:
    double mP[4];
    Vector3D mX;        // Orthonormal base X
    Vector3D mY;        // Orthonormal base Y
    Vector3D mNormal;   // Orthonormal base Z: must always be normalized
#pragma endregion
};
}
