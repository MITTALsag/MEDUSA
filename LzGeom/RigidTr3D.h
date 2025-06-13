#pragma once

#include <LzMath/Matrix.h>
//#include <LzMath/Quaternion.h>
#include <string>


namespace LzGeom
{
using LzMath::Matrix;
using std::string;


class Coord3D;
class Point3D;
class Vector3D;
class Plane3D;
class Line3D;


class RigidTr3D
{
//*** Factory plutot que constructeur? ==> si oui, alors copy const avec move / copy elision
//*** Marche avec types alias (et pas vraies classes)
//public:
//    static RigidTr3D Build( std::function<RigidTr3D(void)> pBuild );
//*** Factory plutot que constructeur? ==> si oui, alors copy const avec move / copy elision

public:
    // Construction / destruction
    RigidTr3D();
    RigidTr3D( std::function<void(RigidTr3D * pThis)> pInit );
    //
    enum class EltsOrder { RowByRow/*RowMajor*/, ColByCol/*ColMajor*/ };
    RigidTr3D( const float pMat[16], EltsOrder pOrder );
    RigidTr3D( const double pMat[16], EltsOrder pOrder );
    RigidTr3D( double pDegX, double pDegY, double pDegZ, const Vector3D & pTr );
    RigidTr3D( double pDeg, const Vector3D & pRotDir, const Vector3D & pTr );

    // X, Y and Z must form a direct orthonormal base otherwise R will not be an orthonormal matrix and 'Inverse' will produce uncertain results
    RigidTr3D( const Vector3D & pX, const Vector3D & pY, const Vector3D & pZ, const Point3D & pO );
    RigidTr3D( const Vector3D & pX, const Vector3D & pY, const Point3D & pO );

    // Operators
    // SetOrthonormal builds an orthonormal ref from unnecessarily orthonormal pX and pY
	void SetOrthonormal( Vector3D pX, Vector3D pY, const Point3D & pO );
    void ResetToZero();
    void ResetToId();
    void SetFromArray( const double pMat[16], EltsOrder pOrder );
    void SetFromArray( const float pMat[16], EltsOrder pOrder );
    //
    void SetRotation( double pDegX, double pDegY, double pDegZ );
    void SetRotation( double pDeg, const Vector3D & pRotDir );
/*
 **** TEST BEFORE USE
 *    template<class T> void SetRotationFromArray( T pMat[9], EltsOrder pOrder )
 *    {
 *        if( pOrder == EltsOrder::RowByRow )
 *        {
 *            for( int i = 0 ; i < 3 ; i++ )
 *            for( int j = 0 ; j < 3 ; j++ )
 *              mR[i * 3 + j] = (double)pMat[i * 3 + j];
 *        }
 *        else
 *        {
 *            for( int i = 0 ; i < 3 ; i++ )
 *            for( int j = 0 ; j < 3 ; j++ )
 *              mR[i * 3 + j] = (double)pMat[j * 3 + i];
 *        }
 *    }
 **** TEST BEFORE USE
 */
    //
    void SetTranslation( const Vector3D & pT );
    Vector3D GetTranslation() const;
    //
    double & T( unsigned int pIJ );
    double T( unsigned int pI ) const;
    //
    double & R( unsigned int pI, unsigned int pJ );
    double R( unsigned int pI, unsigned int pJ ) const;
    double R_Determinant() const;
    //
    void GetRotDirAndAngle( Vector3D & pDir, double & pDeg ) const;
    //
    void GetEulerAngles( double & pDegX, double & pDegY, double & pDegZ ) const;
    //
    // /!\ Linearly interpolating Euler angles is dangerous!
    // pProp: proportion between 0 and 1. 0: returns *this, 1: returns other.
    // Any intermediate value returns linear interpolation between the two.
    RigidTr3D EulerInterpolate( double pProp, const RigidTr3D & pToOther ) const;
    //
    // Interpolates using rotation direction and angle.
    // This might also produce artefacts.
    RigidTr3D RotDirInterpolate( double pProp, const RigidTr3D & pToOther ) const;
    //
    bool IsIdentity( double pTol ) const;
    bool IsSameAs( const RigidTr3D & pOther, double pTol ) const;
    RigidTr3D Inverse() const;
    void MultiplyTo( const Coord3D & pIn, Coord3D & pOut ) const;
    Point3D operator*( const Point3D & pA ) const;
    Vector3D operator*( const Vector3D & pV ) const;
    Plane3D operator*( const Plane3D & pP ) const;
    Line3D operator*( const Line3D & pL ) const;
    RigidTr3D operator*( const RigidTr3D & pT ) const;
    void ToMatrix4x4( Matrix & pMat ) const;
    void ToMatrix3x3( Matrix & pMat ) const;
//    void ToQuaternion( Quaternion& pQuat ) const;
    void FromMatrix4x4( const Matrix &pMat, bool pCheckRotDet = true );
    void FromMatrix3x3( const Matrix &pMat, bool pCheckRotDet = true );
#if !defined(NO_OPENGL)
    void glMultMatrixd() const;
#endif

    // Debug
    void Log( const string & pTxt = "" ) const;
    string ToOldString() const;
    string ToNewString( const string & pSep=" ", double pEpsilon=-1.0 ) const;
    void FromNewString( const string & pStr, bool pCheckRotDet = true, double pTol = -1 );
    void FromOldString( const string & pStr, bool pCheckRotDet = true, double pTol = -1 );

//protected:
	// Data
    double mR[9]; // Row major storage order: R00, R01, R02, R10, R11, R12, ...
    //
    //           R00  R01  R02
    // where R = R10  R11  R12
    //           R20  R21  R22
    //
    double mT[3];
};
}
