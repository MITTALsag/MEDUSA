#pragma once

#include "LzLib_Geom.h"
#include <LzMath/Matrix.h>
#include <string>


// DLL exports
namespace LzGeom
{
    class DLL_EXPORT_LzGeom Coord3D;
}


namespace LzGeom
{
using LzMath::Matrix;
using std::string;


class Coord3D
{
#pragma region "Construction & destruction"
public:
    Coord3D( double pX, double pY, double pZ, double pW );
#pragma endregion


#pragma region "Operators"
public:
    bool operator==( const Coord3D & pOther ) const
    {
        return mV[0] == pOther.mV[0]
            && mV[1] == pOther.mV[1]
            && mV[2] == pOther.mV[2];
    }
    bool operator!=( const Coord3D & pOther ) const
    {
        return mV[0] != pOther.mV[0]
            || mV[1] != pOther.mV[1]
            || mV[2] != pOther.mV[2];
    }
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
//      Coord3D operator*( double pK ) const;
//      friend Coord3D operator*( double pK, const Coord3D & pV );
//
//// To Coord3D ??
//      friend Coord3D operator*( const Matrix & pMat, const Coord3D & pV );
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??

    // Index of dimension with maximal absolute value
    unsigned int MaxAbsValIdx() const;
    
    // Remplacer par un operateur * ???
    void MultiplyTo( const Matrix & pMat, Coord3D & pV ) const; // pV = pMat * this

    // Trace
    string ToString() const;

    // Specific access
    double X() const { return mV[0]; }
    double Y() const { return mV[1]; }
    double Z() const { return mV[2]; }
    double W() const { return mV[3]; }
    double & X() { return mV[0]; }
    double & Y() { return mV[1]; }
    double & Z() { return mV[2]; }
    double & W() { return mV[3]; }
    //
    double operator[]( unsigned int pIndex ) const;
    double & operator[]( unsigned int pIndex );
#pragma endregion


#pragma region "Data"
public:
    double mV[4];
#pragma endregion
};

// Fast operators
#define _LzGeomCoord3D_A_is_B_add_C_mult_D(A,B,C,D) \
{                           \
    (A).mV[0] = (B).mV[0] + (C)*(D).mV[0];      \
    (A).mV[1] = (B).mV[1] + (C)*(D).mV[1];      \
    (A).mV[2] = (B).mV[2] + (C)*(D).mV[2];      \
}
}
