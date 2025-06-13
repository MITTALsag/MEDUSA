#pragma once

#include "LzLib_Geom.h"
#include "Coord3D.h"


// DLL exports
namespace LzGeom
{
    class DLL_EXPORT_LzGeom Point3D;
}


namespace LzGeom
{
class Vector3D;
class RigidTr3D;


class Point3D : public Coord3D
{
public:
    // Construction / destruction
    Point3D();
    Point3D( double pX, double pY, double pZ );
    Point3D( const double pV[3] );


#pragma region "Xml parsing"
#ifdef USE_CLI
    virtual void ReadFromXml( XmlTextReader ^pReader ) override { Coord3D::ReadFromXml( "Point3D", pReader ); }
    virtual void WriteToXml( XmlWriter ^pWriter ) const override { Coord3D::WriteToXml( "Point3D", pWriter ); }
#endif
#pragma endregion


    // Operators
    void Reset()
    {
        mV[0] = mV[1] = mV[2] = 0.0;
        mV[3] = 1.0;
    }
    double SquaredDistanceTo( const Point3D & pB ) const;
    double DistanceTo( const Point3D & pB ) const;
    Point3D MidPointTo( const Point3D & pB ) const;
    Vector3D operator-( const Point3D & pB ) const;
    Point3D operator+( const Vector3D & pV ) const;
    Point3D operator+( const Point3D & pV ) const;
    Point3D operator-( const Vector3D & pV ) const;
    Point3D operator*( double pScalar ) const;
    const Point3D & operator+=( const Vector3D & pV );
    const Point3D & operator-=( const Vector3D & pV );
    const Point3D & operator*=( const RigidTr3D & pT );
    const Point3D & operator*=( const Matrix & pMat );
    const Point3D & operator*=( double pScale );
    const Point3D & operator/=( double pScale );
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
    friend Point3D operator*( const Matrix & pMat, const Point3D & pA );
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
};
}
