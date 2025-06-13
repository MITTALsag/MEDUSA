#pragma once

#include "Coord3D.h"
#include <LzServices/Vector.h>


namespace LzGeom
{
using LzServices::Vector;


class RigidTr3D;


#ifdef USE_CLI
using System::String;
#endif


class Vector3D : public Coord3D
{
#pragma region "Construction & destruction"
public:
    Vector3D();
    Vector3D( double pX, double pY, double pZ );
    Vector3D( const double pV[3] );
    #pragma endregion


#pragma region "Xml parsing"
#ifdef USE_CLI
public:
    virtual void ReadFromXml( XmlTextReader ^pReader ) override { Coord3D::ReadFromXml( "Vector3D", pReader ); }
    virtual void WriteToXml( XmlWriter ^pWriter ) const override { Coord3D::WriteToXml( "Vector3D", pWriter ); }
#endif
#pragma endregion


#pragma region "Operations"
    void Reset()
    {
        mV[0] = mV[1] = mV[2] = mV[3] = 0.0;
    }
    double operator*( const Vector3D & pV ) const;
    double SquaredNorm() const;
    double Norm() const;
    void Normalize();
    bool Normalize_NoExc( bool pTraceErr=true );

    void TensorProd( const Vector3D & pOther, Matrix & pTo ) const;
    Vector3D MeanVectorWith( const Vector3D & pV ) const;
    Vector3D operator^( const Vector3D & pV ) const;
    Vector3D operator/( double pK ) const;
    Vector3D operator/=( double pK );
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
    Vector3D operator*( double pK ) const;
    friend Vector3D operator*( double pK, const Vector3D & pV );
    friend Vector3D operator*( const Matrix & pMat, const Vector3D & pV );

//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
    Vector3D operator+() const;
    Vector3D operator-() const;
    Vector3D operator+( const Vector3D & pV ) const;
    Vector3D operator-( const Vector3D & pV ) const;
    const Vector3D & operator+=( const Vector3D & pV );
    const Vector3D & operator-=( const Vector3D & pV );
    const Vector3D & operator*=( double pK );

    const Vector3D & operator*=( const RigidTr3D & pT );
    const Vector3D & operator*=( const Matrix & pMat );
    Vector3D AffineMult_3x3( const Matrix & pMat ) const; // Preserves orthogonality in affine transform. If pMat is a rotation then same effect as operator*.
    Vector3D AffineMult_4x4( const Matrix & pMat ) const; // Preserves orthogonality in affine transform. If pMat is a rotation then same effect as operator*.

    // pZ = (*this) vector normalized
    void BuildOrthonormalBase( Vector3D & pX, Vector3D & pY, Vector3D & pZ ) const;

#pragma region "3D angles"
    // The positive orientation of plane (this,pB) is defined by pPositiveOrient.

    // Returns a signed angle in radians in [-Pi,+Pi].
    double SignedAngleTo( const Vector3D & pB, const Vector3D & pPositiveOrient ) const;
    bool SignedAngleTo_NoExc( const Vector3D & pB, const Vector3D & pPositiveOrient, double & pAngle ) const;

    // Returns a positive angle in radians in [0,2Pi].
    double PositiveAngleTo( const Vector3D & pB, const Vector3D & pPositiveOrient ) const;
    bool PositiveAngleTo_NoExc( const Vector3D & pB, const Vector3D & pPositiveOrient, double & pAngle ) const;

    // Returns a positive angle in radians in [0,Pi].
    double ShortPositiveAngleTo( const Vector3D & pB ) const;
    bool ShortPositiveAngleTo_NoExc( const Vector3D & pB, double & pAngle ) const;

	// Correlation
	// Finds max absolute scalar product of this with a bunch of vectors
	unsigned int StrongestAbsScalar( const Vector<Vector3D> & pDirs ) const;
#pragma endregion
#pragma endregion
};
}
