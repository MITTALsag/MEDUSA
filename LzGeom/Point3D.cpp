#include "Point3D.h"
#include "Vector3D.h"
#include "RigidTr3D.h"
#include <LzServices/LzLog.h>
#include <cmath>


namespace LzGeom
{
//================================================================================
Point3D::Point3D() : Coord3D( 0, 0, 0, 1 )
{
}

//================================================================================
Point3D::Point3D( double pX, double pY, double pZ ) : Coord3D( pX, pY, pZ, 1 )
{
}

//================================================================================
Point3D::Point3D( const double pV[3] ) : Coord3D( pV[0], pV[1], pV[2], 1 )
{
}

//================================================================================
double Point3D::SquaredDistanceTo( const Point3D & pB ) const
{
    double DX = X() - pB.X();
    double DY = Y() - pB.Y();
    double DZ = Z() - pB.Z();
    return DX * DX + DY * DY + DZ * DZ;
}

//================================================================================
Point3D Point3D::MidPointTo( const Point3D & pB ) const
{
    return Point3D( 0.5*( mV[0] + pB.mV[0] ),
                    0.5*( mV[1] + pB.mV[1] ),
                    0.5*( mV[2] + pB.mV[2] ) );
}

//================================================================================
double Point3D::DistanceTo( const Point3D & pB ) const
{
    return sqrt( SquaredDistanceTo( pB ) );
}

//================================================================================
Vector3D Point3D::operator-( const Point3D & pB ) const
{
    return Vector3D( X() - pB.X(), Y() - pB.Y(), Z() - pB.Z() );
}

//================================================================================
Point3D Point3D::operator+( const Vector3D & pV ) const
{
    return Point3D( X() + pV.X(), Y() + pV.Y(), Z() + pV.Z() );
}

//================================================================================
Point3D Point3D::operator+( const Point3D & pP ) const
{
    return Point3D( X() + pP.X(), Y() + pP.Y(), Z() + pP.Z() );
}

//================================================================================
Point3D Point3D::operator-( const Vector3D & pV ) const
{
    return Point3D( X() - pV.X(), Y() - pV.Y(), Z() - pV.Z() );
}

//================================================================================
const Point3D & Point3D::operator+=( const Vector3D & pV )
{
    X() += pV.X();
    Y() += pV.Y();
    Z() += pV.Z();
    return *this;
}

//================================================================================
const Point3D & Point3D::operator-=( const Vector3D & pV )
{
    X() -= pV.X();
    Y() -= pV.Y();
    Z() -= pV.Z();
    return *this;
}

//================================================================================
const Point3D & Point3D::operator*=( const RigidTr3D & pT )
{
    *this = pT * *this;
    return *this;
}

//================================================================================
const Point3D & Point3D::operator*=( const Matrix & pMat )
{
    *this = pMat * (*this);
    return *this;
}

//================================================================================
const Point3D & Point3D::operator*=( double pScale )
{
    mV[0] *= pScale;
    mV[1] *= pScale;
    mV[2] *= pScale;

    return *this;
}

//================================================================================
const Point3D & Point3D::operator/=( double pScale )
{
    mV[0] /= pScale;
    mV[1] /= pScale;
    mV[2] /= pScale;

    return *this;
}

//================================================================================
Point3D operator*( const Matrix & pMat, const Point3D & pV )
{
    Point3D lProd;
    pV.MultiplyTo( pMat, lProd );
    return lProd;
}

//================================================================================
Point3D Point3D::operator*( double pScalar ) const
{
    Point3D lPt( *this );
    lPt *= pScalar;

    return lPt;
}

}
