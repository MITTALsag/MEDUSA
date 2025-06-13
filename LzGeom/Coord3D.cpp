#include "Coord3D.h"
#include <LzServices/LzLog.h>
#include <iomanip> // std::precision
//#include <cstdlib>
#include <cmath>


namespace LzGeom
{
#pragma region "Construction & destruction"
//================================================================================
Coord3D::Coord3D( double pX, double pY, double pZ, double pW )
{
    mV[0] = pX;
    mV[1] = pY;
    mV[2] = pZ;
    mV[3] = pW;
}
#pragma endregion


#pragma region "Operators"
//================================================================================
//  Coord3D Coord3D::operator*( double pK ) const
//  {
//      return Coord3D( X()*pK, Y()*pK, Z()*pK, W() ); //************ W()*pK
//  }
//
//================================================================================
//  Coord3D operator*( double pK, const Coord3D & pV )
//  {
//      return pV*pK;
//  }
//
//================================================================================
//  Coord3D operator*( const Matrix & pMat, const Coord3D & pV )
//  {
//#if 1
//      Coord3D lProd;
//      if( pMat.Rows()!=pMat.Cols() || (pMat.Rows()!=3 && pMat.Rows()!=4) )
//      {
//          LzLogErr("", "Invalid matrix dimensions! ("+pMat.Rows()+"x"+pMat.Cols()+")");
//          return lProd;
//      }
//
//      // Mult
//      for( unsigned int i=0 ; i<pMat.Rows() ; i++ )
//      for( unsigned int j=0 ; j<pMat.Cols() ; j++ )
//          lProd.mV[i] += pMat.Elt(i,j) * pV.mV[j];
//#pragma message("***** a tester : Vector3D operator*( const Matrix & pMat, const Vector3D & pV )")
//#else
//      Vector3D lProd;
//      pV.MultiplyTo( pMat, lProd );
//#endif
//      return lProd;
//  }

//================================================================================
unsigned int Coord3D::MaxAbsValIdx() const
{
    // Assuming 0 is the max
    unsigned int lMaxIdx = 0;
    double lMaxAbs = std::abs( mV[0] );

    // 1 is better?
    const double lAbs1 = std::abs( mV[1] );
    if( lMaxAbs < lAbs1 )
    {
        lMaxIdx = 1;
        lMaxAbs = lAbs1;
    }

    // 2 is better?
    const double lAbs2 = std::abs( mV[2] );
    if( lMaxAbs < lAbs2 )
		return 2;

    // Return whatever was found before
    return lMaxIdx;
}

//================================================================================
void Coord3D::MultiplyTo( const Matrix & pMat, Coord3D & pV ) const
{
    // Check
    if( pMat.Rows()!=pMat.Cols() || (pMat.Rows()!=3 && pMat.Rows()!=4) )
        LzLogException("", "Invalid matrix dimensions! (" << pMat.Rows() << "x" << pMat.Cols() << ")")

    // Reset all coordinates
    pV.mV[0] = pV.mV[1] = pV.mV[2] = 0;

    // If 3x3 matrix then preserve the last coordinate, else reset
    if( pMat.Rows() == 3 )
        pV.mV[3] =  mV[3];
    else
        pV.mV[3] =  0;

    // Mult
    for( unsigned int i = 0 ; i < pMat.Rows() ; i++ )
    for( unsigned int j = 0 ; j < pMat.Cols() ; j++ )
        pV.mV[i] += pMat.Elt( i, j ) * mV[j];
}

//================================================================================
string Coord3D::ToString() const
{
	std::stringstream lStr;
    lStr << std::setprecision(15) << "[ " << mV[0] << ", " << mV[1] << ", " << mV[2] << ", " << mV[3] << " ]";

	return lStr.str();
}

//================================================================================
double Coord3D::operator[]( unsigned int pIndex ) const
{
    if ( pIndex > 3 )
        LzLogException("", "Index "<<pIndex<<"out of bounds")

    return mV[pIndex];
}

//================================================================================
double & Coord3D::operator[]( unsigned int pIndex )
{
    if ( pIndex > 3 )
        LzLogException("", "Index "<<pIndex<<"out of bounds")

    return mV[pIndex];
}
}
