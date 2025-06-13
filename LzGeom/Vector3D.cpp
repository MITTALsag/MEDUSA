#include "Vector3D.h"
#include "Point3D.h"
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>
#include <cstdlib>
#include <cmath>


namespace LzGeom
{
#pragma region "Construction & destruction"
//================================================================================
Vector3D::Vector3D() : Coord3D( 0, 0, 0, 0 )
{
}

//================================================================================
Vector3D::Vector3D( double pX, double pY, double pZ ) : Coord3D( pX, pY, pZ, 0 )
{
}

//================================================================================
Vector3D::Vector3D( const double pV[3] ) : Coord3D( pV[0], pV[1], pV[2], 0 )
{
}
#pragma endregion


#pragma region "Operations"
//================================================================================
double Vector3D::operator*( const Vector3D & pV ) const
{
    return X() * pV.X() + Y() * pV.Y() + Z() * pV.Z();
}

//================================================================================
double Vector3D::SquaredNorm() const
{
    return operator*( *this );
}

//================================================================================
double Vector3D::Norm() const
{
    return std::sqrt( SquaredNorm() );
}

//================================================================================
void Vector3D::Normalize()
{
    double lNorm = Norm();

	// Check
    if( LzMath::ToolBox::IsZero(lNorm) )
        LzLogException("", "Null vector!");

    // Normalize
    mV[0] /= lNorm;
    mV[1] /= lNorm;
    mV[2] /= lNorm;
}

//================================================================================
bool Vector3D::Normalize_NoExc( bool pTraceErr/*=true*/ )
{
    double lNorm = Norm();

	// Check
    if( LzMath::ToolBox::IsZero(lNorm) )
	{
        if( pTraceErr )
            LzLogErr("", "Null vector! Leaving unchanged, returning 'false'...")
		return false;
	}

    // Normalize
    mV[0] /= lNorm;
    mV[1] /= lNorm;
    mV[2] /= lNorm;

	return true;
}

//================================================================================
void Vector3D::TensorProd( const Vector3D & pOther, Matrix & pTo ) const
{
    pTo.SetDims( 3, 3 );
    pTo.LoadValue( 0 );

    for( int i = 0 ; i < 3 ; i++ )
        for( int j = 0 ; j < 3 ; j++ )
            pTo.Elt( i, j ) += mV[i] * pOther.mV[j];
}

//================================================================================
Vector3D Vector3D::MeanVectorWith( const Vector3D & pV ) const
{
    return Vector3D( 0.5 * ( mV[0] + pV.mV[0] ), 0.5 * ( mV[1] + pV.mV[1] ), 0.5 * ( mV[2] + pV.mV[2] ) );
}

//================================================================================
Vector3D Vector3D::operator^( const Vector3D & pV ) const
{
    return Vector3D( Y() * pV.Z() - Z() * pV.Y(), Z() * pV.X() - X() * pV.Z(), X() * pV.Y() - Y() * pV.X() );
}

//================================================================================
Vector3D Vector3D::operator/( double pK ) const
{
	// Check
    if( LzMath::ToolBox::IsZero( pK ) )
        LzLogException("",  "Dividing vector by 0!")

	return Vector3D( X() / pK, Y() / pK, Z() / pK );
}

//================================================================================
Vector3D Vector3D::operator/=( double pK )
{
	// Check
    if( LzMath::ToolBox::IsZero( pK ) )
        LzLogException("",  "Dividing vector by 0!")

	this->X() /= pK;
    this->Y() /= pK;
    this->Z() /= pK;

    return *this;
}

//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
//================================================================================
Vector3D Vector3D::operator*( double pK ) const
{
    return Vector3D( X() * pK, Y() * pK, Z() * pK );
}

//================================================================================
Vector3D operator*( double pK, const Vector3D & pV )
{
    return pV * pK;
}

//================================================================================
Vector3D operator*( const Matrix & pMat, const Vector3D & pV )
{
#if 1
    Vector3D lProd;
    pV.MultiplyTo( pMat, lProd );
    return lProd;
#else
    // Check
    if( pMat.Rows() != pMat.Cols() || ( pMat.Rows() != 3 && pMat.Rows() != 4 ) )
        LzLogException("",  "Invalid matrix dimensions! (" + pMat.Rows() + "x" + pMat.Cols() + ")" );

    Vector3D lProd;
    for( unsigned int i = 0 ; i < pMat.Rows() ; i++ )
        for( unsigned int j = 0 ; j < pMat.Cols() ; j++ )
            lProd.mV[i] += pMat.Elt( i, j ) * pV.mV[j];

    return lProd;
#endif
}
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??
//**************************************** Move code to Coord3D ??

//================================================================================
Vector3D Vector3D::operator+() const
{
    return Vector3D( +X(), +Y(), +Z() );
}

//================================================================================
Vector3D Vector3D::operator-() const
{
    return Vector3D( -X(), -Y(), -Z() );
}

//================================================================================
Vector3D Vector3D::operator+( const Vector3D & pV ) const
{
    return Vector3D( X() + pV.X(), Y() + pV.Y(), Z() + pV.Z() );
}

//================================================================================
Vector3D Vector3D::operator-( const Vector3D & pV ) const
{
    return Vector3D( X() - pV.X(), Y() - pV.Y(), Z() - pV.Z() );
}

//================================================================================
const Vector3D & Vector3D::operator+=( const Vector3D & pV )
{
    X() += pV.X();
    Y() += pV.Y();
    Z() += pV.Z();
    return *this;
}

//================================================================================
const Vector3D & Vector3D::operator-=( const Vector3D & pV )
{
    X() -= pV.X();
    Y() -= pV.Y();
    Z() -= pV.Z();
    return *this;
}

//================================================================================
const Vector3D & Vector3D::operator*=( double pK )
{
    X() *= pK;
    Y() *= pK;
    Z() *= pK;
    return *this;
}

//================================================================================
const Vector3D & Vector3D::operator*=( const RigidTr3D & pT )
{
    *this = pT * *this;
    return *this;
}

//================================================================================
const Vector3D & Vector3D::operator*=( const Matrix & pMat )
{
    *this = pMat * *this;
    return *this;
}

//================================================================================
Vector3D Vector3D::AffineMult_3x3( const Matrix & pMat ) const
{
    // Check
    if( pMat.Rows() != pMat.Cols() || pMat.Rows() != 3 )
        LzLogException("",  "Invalid matrix dimensions! (" << pMat.Rows() << "x" << pMat.Cols() << ")" );

    Matrix B, C;

    // Compute transform for normals: C = ( (pMat)_3x3 )^(-t)
    pMat.M3x3_InverseTo( B );
    B.TransposeTo( C );

    return C * *this;
}

//================================================================================
Vector3D Vector3D::AffineMult_4x4( const Matrix & pMat ) const
{
    // Check
    if( pMat.Rows() != pMat.Cols() || pMat.Rows() != 4 )
        LzLogException("",  "Invalid matrix dimensions! (" << pMat.Rows() << "x" << pMat.Cols() << ")" );

    Matrix A( 3, 3 ), B, C;

    // Need to extract 3x3 submatric first
    pMat.SubMatrixTo( 0, 0, A );

    // Compute transform for normals: C = ( (pMat)_3x3 )^(-t)
    A.M3x3_InverseTo( B );
    B.TransposeTo( C );

    return C * *this;
}

//================================================================================
void Vector3D::BuildOrthonormalBase( Vector3D & pX, Vector3D & pY, Vector3D & pZ ) const
{
    // Normalize pZ
    pZ = *this;
    pZ.Normalize();

    // Get rid of smallest coord among X and Y
    int iNull = std::abs( pZ.mV[0] ) < std::abs( pZ.mV[1] ) ? 0 : 1 ;

    pX.mV[iNull]       = 0;
    pX.mV[( iNull + 1 ) % 3] = +pZ.mV[( iNull + 2 ) % 3];
    pX.mV[( iNull + 2 ) % 3] = -pZ.mV[( iNull + 1 ) % 3];

    // Normalize pX
    pX.Normalize();

    // pY is normalized
    pY = pZ ^ pX;
}


#pragma region "3D angles"
//================================================================================
double Vector3D::SignedAngleTo( const Vector3D & pB, const Vector3D & pPositiveOrient ) const
{
    double lNormAB = Norm() * pB.Norm();

    if( LzMath::ToolBox::IsZero( lNormAB ) )
        LzLogException("", "One of the vectors is null!" )

    // Cos
    double Cos = (*this * pB) / lNormAB;

    // Sin
    Vector3D lN = *this ^ pB;
    double Sin = lN.Norm() / lNormAB;
    if( lN * pPositiveOrient < 0 )
        Sin *= -1;

    // Angle
    if( Cos >= +LzMath::LzEpsilon )
    {
        // Cos is > 0 (and non epsilonesque)
        // Angle is between -90 and +90 deg
        return std::atan( Sin / Cos );
    }
    else
    if( Cos <= -LzMath::LzEpsilon )
    {
        // Cos is < 0 (and non epsilonesque)
        // Angle is between 90 and 180, or between -180 and -90
        double beta = std::atan( Sin / Cos );
        return Sin>=0 ? beta+LzMath::PI : beta-LzMath::PI ;
    }
    else
    {
        // Cos is 0 (epsilonesque), angle is either +180 or -180
        return Sin>=0 ? +LzMath::PI/2 : -LzMath::PI/2 ;
    }
}
#if 0 //*** From Twinlib
//================================================================================
double SignedAngleTo( const Vector3 & pA, const Vector3 & pB, const Vector3 & pPositiveOrient )
{
  // Norm product
  const double lNormAB = pA.norm() * pB.norm();

  // Check whether one of the vectors is 0
  if( IsZero(lNormAB, EPS) )
    THROW_EXCEPTION("Could not compute signed angle! One of the vectors is null.");

  // Check
  if( pPositiveOrient.norm() < EPS )
    THROW_EXCEPTION("Could not compute signed angle! Orientation vector is null.");

  // Cos
  const double Cos = pA.dot(pB) / lNormAB;

  // A and B are collinear?
  if( Cos > +1.0 - EPS )
  {
    // A and B are positively aligned
    return 0;
  }
  else
  if( Cos < -1.0 + EPS )
  {
    // A and B are negatively aligned
    return +PI;
  }
  else
  {
    // A and B are NOT collinear

    // Check coplanarity
    const Vector3 lN = pA.cross(pB); // Non null since A and B are not collinear
    const double lNorm_dot_Orient = lN.dot(pPositiveOrient);
    //
    if( std::abs(lNorm_dot_Orient) < EPS )
      THROW_EXCEPTION("Could not compute signed angle! Orientation vector is coplanar with the 2 vectors.")

    // Sin
    double Sin = lN.norm() / lNormAB; // Positive sine; the sign is adjusted below
    if( lNorm_dot_Orient < 0 )
      Sin *= -1;

    // Right angle?
    if( std::abs(Cos) < EPS )
    {
      return Sin >= 0 ? +PI/2 : -PI/2;
    }
    else
    {
      return std::atan2( Sin, Cos );
    }
  }
}
#endif

//================================================================================
bool Vector3D::SignedAngleTo_NoExc( const Vector3D & pB, const Vector3D & pPositiveOrient, double & pAngle ) const
{
    double lNormAB = Norm() * pB.Norm();

    if( LzMath::ToolBox::IsZero( lNormAB ) )
        return false;

    // Cos
    double Cos = ( *this * pB ) / lNormAB;

    // Sin
    Vector3D lN = *this ^ pB;
    double Sin = lN.Norm() / lNormAB;
    if( lN * pPositiveOrient < 0 )
        Sin *= -1;

    // Angle
    if( Cos >= +LzMath::LzEpsilon )
        pAngle = std::atan( Sin / Cos );
    else
    if( Cos <= -LzMath::LzEpsilon )
    {
        double beta = std::atan( Sin / Cos );
        pAngle = Sin>0 ? beta+LzMath::PI : beta-LzMath::PI ;
    }
    else
        pAngle = Sin>0 ? +LzMath::PI/2 : -LzMath::PI/2 ;

	return true;
}

//================================================================================
double Vector3D::PositiveAngleTo( const Vector3D & pB, const Vector3D & pPositiveOrient ) const
{
    double lSigned = SignedAngleTo( pB, pPositiveOrient );
    return lSigned >= 0 ? lSigned : 2 * LzMath::PI + lSigned ;
}

//================================================================================
bool Vector3D::PositiveAngleTo_NoExc( const Vector3D & pB, const Vector3D & pPositiveOrient, double & pAngle ) const
{
    double lSigned;
	bool lRes = SignedAngleTo_NoExc( pB, pPositiveOrient, lSigned );

if( !lRes )
	return false;

    pAngle = lSigned >= 0 ? lSigned : LzMath::DPI + lSigned ;
	return true;
}

//================================================================================
double Vector3D::ShortPositiveAngleTo( const Vector3D & pB ) const
{
    //
    // NB: should be implemented using arccos( this \dot pB / (normThis * normB) )
    //
    double lPosAngle = PositiveAngleTo( pB, (*this)^pB );

    if( lPosAngle <= LzMath::PI )
        return lPosAngle;
    else
        return LzMath::DPI - lPosAngle;
}

//================================================================================
bool Vector3D::ShortPositiveAngleTo_NoExc( const Vector3D & pB, double & pAngle ) const
{
    double lPosAngle;
	bool lRes = PositiveAngleTo_NoExc( pB, ( *this )^pB, lPosAngle );

if( !lRes )
	return false;

    if( lPosAngle <= LzMath::PI )
        pAngle = lPosAngle;
    else
        pAngle = 2 * LzMath::PI - lPosAngle;

	return true;
}

//================================================================================
unsigned int Vector3D::StrongestAbsScalar( const Vector<Vector3D> & pDirs ) const
{
	// Check
	if( pDirs.Size() == 0 )
        LzLogException("", "Cannot find strongest scalar product with 0 vectors!");

	// Find max absolute scalar product of this with a bunch of vectors
	double lMaxAbsScal = 0.0;
	unsigned int lBestIdx = 0;
	for( unsigned int d=0 ; d<pDirs.Size() ; d++ )
	{
		// Compute scalar product and update
		double lAbsScal = fabs( *this * pDirs[d] );
		if( lMaxAbsScal < lAbsScal )
		{
			lMaxAbsScal = lAbsScal;
			lBestIdx = d;
		}
	}

	return lBestIdx;
}
#pragma endregion
#pragma endregion
}
