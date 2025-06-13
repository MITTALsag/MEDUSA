#include "RigidTr3D.h"
#include "Coord3D.h"
#include "Point3D.h"
#include "Vector3D.h"
#include "Plane3D.h"
#include "Line3D.h"
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>
#if defined(USING_QT) && !defined(NO_OPENGL)
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/gl.h>
#endif
#if defined(USING_QT)
    #include <QString>
    #include <QRegExp>
    #include <QStringList>
#endif
#include <cmath>
#include <iomanip> // std::precision


namespace LzGeom
{
using LzMath::Matrix;


//================================================================================
static const double DEG2RAD = LzMath::PI / 180.0;
static const double RAD2DEG = 180.0 / LzMath::PI;

//================================================================================
RigidTr3D::RigidTr3D()
{
    // Set id
    ResetToId();
}

//================================================================================
RigidTr3D::RigidTr3D( std::function<void(RigidTr3D * pThis)> pInit )
{
    // Set id, so that pInit has a well defined Id transform to begin with
    ResetToId();

    // Call initializer
    pInit(this);
}

//================================================================================
RigidTr3D::RigidTr3D( const float pMat[16], EltsOrder pOrder )
{
    SetFromArray( pMat, pOrder );
}

//================================================================================
RigidTr3D::RigidTr3D( const double pMat[16], EltsOrder pOrder )
{
    SetFromArray( pMat, pOrder );
}

//================================================================================
RigidTr3D::RigidTr3D( double pDegX, double pDegY, double pDegZ, const Vector3D & pTr )
{
    SetRotation( pDegX, pDegY, pDegZ );
    mT[0] = pTr.X();
    mT[1] = pTr.Y();
    mT[2] = pTr.Z();
}

//================================================================================
RigidTr3D::RigidTr3D( double pDeg, const Vector3D & pRotDir, const Vector3D & pTr )
{
    SetRotation( pDeg, pRotDir );
    mT[0] = pTr.X();
    mT[1] = pTr.Y();
    mT[2] = pTr.Z();
}

//================================================================================
RigidTr3D::RigidTr3D( const Vector3D & pX, const Vector3D & pY, const Vector3D & pZ, const Point3D & pO )
{
    // XYZT
    for( int j=0 ; j<3 ; j++ )
    {
        mR[j*3 + 0] = pX.mV[j];
        mR[j*3 + 1] = pY.mV[j];
        mR[j*3 + 2] = pZ.mV[j];
        mT[j] = pO.mV[j];
    }
}

//================================================================================
RigidTr3D::RigidTr3D( const Vector3D & pX, const Vector3D & pY, const Point3D & pO )
{
    // Compute missing Z
    const Vector3D lZ = pX ^ pY;

    // XYZT
    for( int j=0 ; j<3 ; j++ )
    {
        mR[j*3 + 0] = pX.mV[j];
        mR[j*3 + 1] = pY.mV[j];
        mR[j*3 + 2] = lZ.mV[j];
        mT[j] = pO.mV[j];
    }
}

//================================================================================
void RigidTr3D::SetOrthonormal( Vector3D pX, Vector3D pY, const Point3D & pO )
{
	// Recover orthonormality
	pX.Normalize();
	//
	Vector3D lZ = pX^pY;
	lZ.Normalize();
	//
	pY = lZ^pX;

	// Set
	*this = RigidTr3D( pX, pY, lZ, pO );
}

//================================================================================
void RigidTr3D::ResetToZero()
{
    // /!\ Invalid rotation
    mR[0] = mR[1] = mR[2] = mR[3] = mR[4] = mR[5] = mR[6] = mR[7] = mR[8] = 0;

    // 0 translation
    mT[0] = mT[1] = mT[2] = 0;
}

//================================================================================
void RigidTr3D::ResetToId()
{
    ResetToZero();

    // Id rotation
    mR[0] = 1;
    mR[4] = 1;
    mR[8] = 1;
}

//================================================================================
void RigidTr3D::SetFromArray( const double pMat[16], EltsOrder pOrder )
{
    if( pOrder == EltsOrder::RowByRow )
    {
        for( int i=0 ; i<3 ; i++ )
        {
            for( int j=0 ; j<3 ; j++ )
                mR[i * 3 + j] = pMat[i * 4 + j];

            mT[i] = pMat[i * 4 + 3];
        }
    }
    else
    {
        for( int i=0 ; i<3 ; i++ )
        {
            for( int j=0 ; j<3 ; j++ )
                mR[i * 3 + j] = pMat[j * 4 + i];

            mT[i] = pMat[12 + i];
        }
    }
}

//================================================================================
void RigidTr3D::SetFromArray( const float pMat[16], EltsOrder pOrder )
{
    if( pOrder == EltsOrder::RowByRow )
    {
        for( int i=0 ; i<3 ; i++ )
        {
            for( int j=0 ; j<3 ; j++ )
                mR[i * 3 + j] = pMat[i * 4 + j];

            mT[i] = pMat[i * 4 + 3];
        }
    }
    else
    {
        for( int i=0 ; i<3 ; i++ )
        {
            for( int j=0 ; j<3 ; j++ )
                mR[i * 3 + j] = pMat[j * 4 + i];

            mT[i] = pMat[12 + i];
        }
    }
}

//================================================================================
void RigidTr3D::SetRotation( double pDegX, double pDegY, double pDegZ )
{
    // Convert to radians
    double lRX = DEG2RAD * pDegX;
    double lRY = DEG2RAD * pDegY;
    double lRZ = DEG2RAD * pDegZ;

    double cx = cos( lRX );
    double sx = sin( lRX );
    double cy = cos( lRY );
    double sy = sin( lRY );
    double cz = cos( lRZ );
    double sz = sin( lRZ );

    double lRotX[9] =
    {
        1,  0,   0,
        0, cx, -sx,
        0, sx,  cx
    };
    Matrix lMX( 3, 3, lRotX );

    double lRotY[] =
    {
        cy, 0, sy,
        0, 1,  0,
        -sy, 0, cy
    };
    Matrix lMY( 3, 3, lRotY );

    double lRotZ[] =
    {
        cz, -sz, 0,
        sz,  cz, 0,
        0,   0, 1
    };
    Matrix lMZ( 3, 3, lRotZ );

    // Combine
    Matrix lRotZY, lRotZYX;
    lMZ.Mult( lMY, lRotZY );
    lRotZY.Mult( lMX, lRotZYX );

    // Set rotation
    FromMatrix3x3( lRotZYX, false );
}

//================================================================================
void RigidTr3D::SetRotation( double pDeg, const Vector3D & pRotDir )
{
    // Compute norm
    double lNorm = pRotDir.Norm();

    // Check
    if( lNorm < 1e-6 )
        LzLogException("", "Rotation direction cannot be a null (or even small) vector! Norm= " << lNorm << ".")

    // Compute stuff
    double X = pRotDir.X() / lNorm;
    double Y = pRotDir.Y() / lNorm;
    double Z = pRotDir.Z() / lNorm;
    //
    double lRad = LzMath::PI * pDeg / 180.0;
    double C = cos( lRad );
    double S = sin( lRad );
    double OmC = 1 - C;
    //
    double XS = X * S;
    double YS = Y * S;
    double ZS = Z * S;
    //
    double XxY = X * Y;
    double XxZ = X * Z;
    double YxZ = Y * Z;

    // Fill-in matrix
    mR[0] = X * X * OmC + C;
    mR[1] = XxY * OmC - ZS;
    mR[2] = XxZ * OmC + YS;
    mR[3] = XxY * OmC + ZS;
    mR[4] = Y * Y * OmC + C;
    mR[5] = YxZ * OmC - XS;
    mR[6] = XxZ * OmC - YS;
    mR[7] = YxZ * OmC + XS;
    mR[8] = Z * Z * OmC + C;
}

//================================================================================
void RigidTr3D::SetTranslation( const Vector3D & pT )
{
    mT[0] = pT.X();
    mT[1] = pT.Y();
    mT[2] = pT.Z();
}

//================================================================================
Vector3D RigidTr3D::GetTranslation() const
{
    return Vector3D( mT[0], mT[1], mT[2] );
}

//================================================================================
double & RigidTr3D::T( unsigned int pI )
{
    // Check
    if( pI >= 3 )
        LzLogException("", "Cannot access translation at "<<pI<<",!")

    return mT[pI];
}

//================================================================================
double RigidTr3D::T( unsigned int pI ) const
{
    return const_cast<RigidTr3D *>(this)->T( pI );
}

//================================================================================
double & RigidTr3D::R( unsigned int pI, unsigned int pJ )
{
    // Check
    if( pI>=3 || pJ>=3 )
        LzLogException("", "Cannot access rotation matrix at "<<pI<<", "<<pJ<<"!")

    return mR[pI*3 + pJ];
}

//================================================================================
double RigidTr3D::R( unsigned int pI, unsigned int pJ ) const
{
    return const_cast<RigidTr3D *>(this)->R( pI, pJ );
}

//================================================================================
double RigidTr3D::R_Determinant() const
{
    Matrix lTr;
    ToMatrix4x4( lTr );

    return lTr.M4x4_Det();
}

//================================================================================
void RigidTr3D::GetRotDirAndAngle( Vector3D & pDir, double & pDeg ) const
{
    // http://www.euclideanspace.com/index.html
    // https://stackoverflow.com/questions/12463487/obtain-rotation-axis-from-rotation-matrix-and-translation-vector-in-opencv

    // Storage order: R00, R01, R02,   R10, R11, R12,  R20, R21, R22
    //                 0    1    2      3    4    5     6    7    8
    const double R00 = mR[0];
    const double R01 = mR[1];
    const double R02 = mR[2];

    const double R10 = mR[3];
    const double R11 = mR[4];
    const double R12 = mR[5];

    const double R20 = mR[6];
    const double R21 = mR[7];
    const double R22 = mR[8];

   // angle = acos(( R00 + R11 + R22 - 1)/2);
    const double lRad = acos( (R00 + R11 + R22 - 1) / 2.0 );
    pDeg = RAD2DEG * lRad;

    // Axis x,y,x:
    const double lDen = sqrt( (R21 - R12)*(R21 - R12) + (R02 - R20)*(R02 - R20) + (R10 - R01)*(R10 - R01) );

    // x = (R21 - R12)/sqrt((R21 - R12)^2+(R02 - R20)^2+(R10 - R01)^2);
    pDir.X() = (R21 - R12) / lDen;

    // y = (R02 - R20)/sqrt((R21 - R12)^2+(R02 - R20)^2+(R10 - R01)^2);
    pDir.Y() = (R02 - R20) / lDen;

    // z = (R10 - R01)/sqrt((R21 - R12)^2+(R02 - R20)^2+(R10 - R01)^2);
    pDir.Z() = (R10 - R01) / lDen;
}

//================================================================================
void RigidTr3D::GetEulerAngles( double & pDegX, double & pDegY, double & pDegZ ) const
{
    if( !LzMath::ToolBox::IsZero( mR[6] - 1 ) && !LzMath::ToolBox::IsZero( mR[6] + 1 ) ) //(mR[6] != -1) && (mR[6] != 1) )
    {
        // two solutions
        //double lTheta1 = -std::asin(mR[6]);
        //double lTheta2 = PI - lTheta1;
        //double lPsi1 = std::atan2((mR[7] / std::cos(lTheta1)), (mR[8] / std::cos(lTheta1)));
        //double lPsi2 = std::atan2((mR[7] / std::cos(lTheta2)), (mR[8] / std::cos(lTheta2)));
        //double lPhi1 = std::atan2((mR[3] / std::cos(lTheta1)), (mR[0] / std::cos(lTheta1)));
        //double lPhi2 = std::atan2((mR[3] / std::cos(lTheta2)), (mR[0] / std::cos(lTheta2)));
        // first solution chosen
        double lRadY = asin( mR[6] );
        pDegY = RAD2DEG * lRadY;
        const double lCosY = cos( lRadY );

        pDegX = RAD2DEG * atan2( ( mR[7] / lCosY ), ( mR[8] / lCosY ) );
        pDegZ = RAD2DEG * atan2( ( mR[3] / lCosY ), ( mR[0] / lCosY ) );
    }
    else
    {
        pDegZ = 0; // can be set to anything
        if( LzMath::ToolBox::IsZero( mR[6] + 1 ) ) // mR[6] == -1
        {
            pDegY = 90;
            pDegX = /*pDegZ +*/ RAD2DEG * atan2( mR[1], mR[2] );
        }
        else
        {
            pDegY = -90;
            pDegX = /*-pDegZ +*/ RAD2DEG * atan2( -mR[1], -mR[2] );
        }
    }
}

//================================================================================
RigidTr3D RigidTr3D::EulerInterpolate( double pProp, const RigidTr3D & pToOther ) const
{
    // Get Euler angles for this and other
    double lDeg_0[3];
    double lDeg_1[3];
    GetEulerAngles( lDeg_0[0], lDeg_0[1], lDeg_0[2] );
    pToOther.GetEulerAngles( lDeg_1[0], lDeg_1[1], lDeg_1[2] );

    // Interpolate angles and translation
    double lInterDeg[3];
    for( int i=0 ; i<3 ; i++ )
        lInterDeg[i] = (1-pProp)*lDeg_0[i] + pProp*lDeg_1[i];
    //
    const Vector3D lTr = (1-pProp)*GetTranslation() + pProp*pToOther.GetTranslation();

    // Build transform
    return RigidTr3D( lInterDeg[0], lInterDeg[1], lInterDeg[2], lTr );
}

//================================================================================
RigidTr3D RigidTr3D::RotDirInterpolate( double pProp, const RigidTr3D & pToOther ) const
{
    // Get rotation direction and angle
    Vector3D lDir_0, lDir_1;
    double lDeg_0, lDeg_1;
    GetRotDirAndAngle( lDir_0, lDeg_0 );
    pToOther.GetRotDirAndAngle( lDir_1, lDeg_1 );

    // Interpolate
    const Vector3D lDir = (1-pProp)*lDir_0 + pProp*lDir_1;
    const double lDeg = (1-pProp)*lDeg_0 + pProp*lDeg_1;
    Vector3D lTr = (1-pProp)*GetTranslation() + pProp*pToOther.GetTranslation();

    return RigidTr3D(lDeg, lDir, lTr);
}

//================================================================================
bool RigidTr3D::IsIdentity( double pTol ) const
{
#if 0
    NO, yields more iterations
    // Adopt Eigen approach, use Frobenius (L2) norm
    // http://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html#ae8443357b808cd393be1b51974213f9c

    // Squared norm
    double lSumSqr = 0.0;

    // Contributions from Rotation
    for( int r=0 ; r<9 ; r++ )
    {
        // Select test
        if( r % 4 )
        {
            lSumSqr += mR[r]*mR[r];
        }
        else
        {
            double diff = mR[r] - 1;
            lSumSqr = diff*diff;
        }
    }

    // Contributions from Translation
    for( int t=0 ; t<3 ; t++ )
        lSumSqr += mT[t]*mT[t];

    // Test
    return sqrt(lSumSqr) <= pTol;
#else
	// Rotation
	for( int r=0 ; r<9 ; r++ )
	{
		// Select test
		if( r % 4 )
		{
			if( fabs(mR[r]) > pTol )
				return false;
		}
		else
		{
			if( fabs(mR[r] - 1) > pTol )
				return false;
		}
	}

	// Translation
	for( int t=0 ; t<3 ; t++ )
	{
		if( fabs(mT[t]) > pTol )
			return false;
	}

	return true;
#endif
}

//================================================================================
bool RigidTr3D::IsSameAs( const RigidTr3D & pOther, double pTol ) const
{
    // Rotation
    for( int r=0 ; r<9 ; r++ )
    {
        if( fabs(mR[r] - pOther.mR[r]) > pTol )
            return false;
    }

    // Translation
    for( int t=0 ; t<3 ; t++ )
    {
        if( fabs(mT[t] - pOther.mT[t]) > pTol )
            return false;
    }

    return true;
}

//================================================================================
RigidTr3D RigidTr3D::Inverse() const
{
    RigidTr3D lInv;

    // Transpose rotation
    lInv.mR[0] = mR[0];
    lInv.mR[1] = mR[3];
    lInv.mR[2] = mR[6];
    lInv.mR[3] = mR[1];
    lInv.mR[4] = mR[4];
    lInv.mR[5] = mR[7];
    lInv.mR[6] = mR[2];
    lInv.mR[7] = mR[5];
    lInv.mR[8] = mR[8];

    // Inverse translation
    lInv.mT[0] = - ( lInv.mR[0] * mT[0] + lInv.mR[1] * mT[1] + lInv.mR[2] * mT[2] );
    lInv.mT[1] = - ( lInv.mR[3] * mT[0] + lInv.mR[4] * mT[1] + lInv.mR[5] * mT[2] );
    lInv.mT[2] = - ( lInv.mR[6] * mT[0] + lInv.mR[7] * mT[1] + lInv.mR[8] * mT[2] );

    return lInv;
}

//================================================================================
void RigidTr3D::MultiplyTo( const Coord3D & pIn, Coord3D & pOut ) const
{
    for( int i = 0 ; i < 3 ; i++ )
    {
        // Translation
        pOut.mV[i] = mT[i] * pIn.mV[3];

        // Rotation
        for( int j = 0 ; j < 3 ; j++ )
            pOut.mV[i] += mR[i * 3 + j] * pIn.mV[j];
    }

    // Copy W
    pOut.mV[3] = pIn.mV[3];
}

//================================================================================
Point3D RigidTr3D::operator*( const Point3D & pA ) const
{
    Point3D lProd;
    MultiplyTo( pA, lProd );
    return lProd;
}

//================================================================================
Vector3D RigidTr3D::operator*( const Vector3D & pV ) const
{
    Vector3D lProd;
    MultiplyTo( pV, lProd );
    return lProd;
}

//================================================================================
Plane3D RigidTr3D::operator*( const Plane3D & pP ) const
{
    Point3D lNewOrigin = (*this) * pP.Projection( Point3D(0, 0, 0) );
    Vector3D lNewNormal = (*this) * pP.Normal();

    return Plane3D( lNewOrigin, lNewNormal );
}

//================================================================================
Line3D RigidTr3D::operator*( const Line3D & pL ) const
{
    return Line3D( ( *this ) * pL.Root(), ( *this ) * pL.Dir() );
}

//================================================================================
RigidTr3D RigidTr3D::operator*( const RigidTr3D & pT ) const
{
    RigidTr3D lProd;

    // Rotation
    for( int i=0 ; i<3 ; i++ )
    for( int j=0 ; j<3 ; j++ )
    {
        lProd.mR[i * 3 + j] = 0;
        for( int k=0 ; k<3 ; k++ )
            lProd.mR[i * 3 + j] += mR[i * 3 + k] * pT.mR[k * 3 + j];
    }

    // Translation
    for( int i=0 ; i<3 ; i++ )
    {
        lProd.mT[i] = mT[i];
        for( int j=0 ; j<3 ; j++ )
            lProd.mT[i] += mR[i * 3 + j] * pT.mT[j];
    }

    return lProd;
}

//================================================================================
void RigidTr3D::ToMatrix4x4( Matrix &pMat ) const
{
    pMat.SetDims( 4, 4 );
    for( int i = 0 ; i < 3 ; i++ )
    {
        for( int j = 0 ; j < 3 ; j++ )
            pMat.Elt( i, j ) = mR[i * 3 + j];

        pMat.Elt( i, 3 ) = mT[i];
    }

    pMat.Elt( 3, 0 ) = pMat.Elt( 3, 1 ) = pMat.Elt( 3, 2 ) = 0;
    pMat.Elt( 3, 3 ) = 1;
}

//================================================================================
void RigidTr3D::ToMatrix3x3( Matrix & pMat ) const
{
    pMat.SetDims( 3, 3 );
    for( int i = 0 ; i < 3 ; i++ )
        for( int j = 0 ; j < 3 ; j++ )
            pMat.Elt( i, j ) = mR[i * 3 + j];
}

//================================================================================
//void RigidTr3D::ToQuaternion( Quaternion & pQuat ) const
//{
//    Matrix lMat( 3, 3 );
//    ToMatrix3x3( lMat );
//    pQuat.FromMatrix ( lMat );
//}

//================================================================================
void RigidTr3D::FromMatrix4x4( const Matrix &pMat, bool pCheckRotDet/*=true*/ )
{
    if( pMat.Rows() != 4 || pMat.Cols() != 4        // A 4x4 matrix, which last row should be:
            || !LzMath::ToolBox::IsZero( pMat.Elt( 3, 0 ) ) // 0
            || !LzMath::ToolBox::IsZero( pMat.Elt( 3, 1 ) ) // 0
            || !LzMath::ToolBox::IsZero( pMat.Elt( 3, 2 ) ) // 0
            || !LzMath::ToolBox::IsZero( pMat.Elt( 3, 3 ) - 1 ) ) // 1
        LzLogException("", "Matrix is not a 3D rigid transform!" )

    // Check rotation ?
    if( pCheckRotDet && !LzMath::ToolBox::IsZero( pMat.M4x4_Det() - 1 ) )
        LzLogException("", "Matrix does not define a direct rigid rotation!" )

    // Copy rotation and translation
    for( int i = 0 ; i < 3 ; i++ )
    {
        for( int j = 0 ; j < 3 ; j++ )
            mR[i * 3 + j] = pMat.Elt( i, j );

        mT[i] = pMat.Elt( i, 3 );
    }
}

//================================================================================
void RigidTr3D::FromMatrix3x3( const Matrix &pMat, bool pCheckRotDet/*=true*/ )
{
    if( pMat.Rows() != 3 || pMat.Cols() != 3 )  // A 3x3 matrix
        LzLogException("", "Matrix is not 3x3!" )

    // Check rotation ?
    if( pCheckRotDet && !LzMath::ToolBox::IsZero( pMat.M3x3_Det() - 1 ) )
        LzLogException("", "Matrix does not define a direct rigid rotation!" )

    // Copy rotation only
    for( int i = 0 ; i < 3 ; i++ )
    {
        for( int j = 0 ; j < 3 ; j++ )
            mR[i * 3 + j] = pMat.Elt( i, j );

        mT[i] = 0;
    }
}

//================================================================================
#if !defined(NO_OPENGL)
void RigidTr3D::glMultMatrixd() const
{
    // double M[16] =
    // {
    //     mR[0], mR[3], mR[6], 0,
    //     mR[1], mR[4], mR[7], 0,
    //     mR[2], mR[5], mR[8], 0,
    //     mT[0], mT[1], mT[2], 1
    // };
    // ::glMultMatrixd( M );
}
#endif

//================================================================================
void RigidTr3D::Log( const std::string & pTxt/*=""*/ ) const
{
    LzLogN("", pTxt )
    LzLogM("", "  | " << mR[0] << " " << mR[1] << " " << mR[2] << " |" )
    LzLogM("", "R=| " << mR[3] << " " << mR[4] << " " << mR[5] << " |" )
    LzLogM("", "  | " << mR[6] << " " << mR[7] << " " << mR[8] << " |" )
    LzLogM("", "" );
    LzLogM("", "  | " << mT[0] << " |" )
    LzLogM("", "T=| " << mT[1] << " |" )
    LzLogM("", "  | " << mT[2] << " |" )
}

//================================================================================
std::string RigidTr3D::ToOldString() const
{
	std::stringstream lStrStr;

    // Rotation
    for( unsigned int i = 0 ; i < 9 ; i++ )
        lStrStr << mR[i] << " ";

    // Translation
    lStrStr << mT[0] << " " << mT[1] << " " << mT[2];

    return lStrStr.str();

//String ^lStr = "";

//// Rotation
//for( unsigned int i = 0 ; i < 9 ; i++ )
//    lStr += mR[i] + " ";

//// Translation
//lStr += mT[0] + " " + mT[1] + " " + mT[2];

//return lStr;
}

//================================================================================
std::string RigidTr3D::ToNewString( const string & pSep/*=" "*/, double pEpsilon/*=-1.0*/ ) const
{
    // Rounding off values to 0
    auto Round = [&]( double pX ) -> double
        {
            if( pEpsilon>0 && std::abs(pX)<=pEpsilon )
                return 0.0;
            else
                return pX;
        };

    // Write transform
	std::stringstream lStrStr;
    lStrStr << std::setprecision(15);
    lStrStr << Round(mR[0]) << pSep << Round(mR[1]) << pSep << Round(mR[2]) << pSep << Round(mT[0]) << pSep <<
               Round(mR[3]) << pSep << Round(mR[4]) << pSep << Round(mR[5]) << pSep << Round(mT[1]) << pSep <<
               Round(mR[6]) << pSep << Round(mR[7]) << pSep << Round(mR[8]) << pSep << Round(mT[2]);

    return lStrStr.str();
}

//================================================================================
void RigidTr3D::FromOldString( const std::string & /*pStdStr*/, bool /*pCheckRotDet*//*=true*/, double /*pTol*//*=-1*/ )
{
#ifdef USE_CLI
	System::String ^pStr = gcnew System::String( pStdStr.c_str() );

	// Parse and check
    cli::array<wchar_t> ^lSpace = { ' ' };
    cli::array<String^> ^lStrs = pStr->Split( lSpace, System::StringSplitOptions::RemoveEmptyEntries );
    if( lStrs->Length != 12 )
        LzLogException("",  "Cannot parse rigid transform 3D! Invalid string: "<<StdStrFromCLI(pStr)<<"." );

    // Try reading
    _LzLogTry
    {
        // Rotation
        for( unsigned int i = 0 ; i < 9 ; i++ )
            mR[i] =LzServices::StrToDouble( lStrs[i] );

        // Translation
        mT[0] = LzServices::StrToDouble( lStrs[ 9] );
        mT[1] = LzServices::StrToDouble( lStrs[10] );
        mT[2] = LzServices::StrToDouble( lStrs[11] );
    }
    _LzLogCatchAndThrow( "Cannot parse rigid transform3D! Syntax error in string: "<<StdStrFromCLI(pStr)<<"." )

    // Check determinant?
    if( pCheckRotDet )
    {
        Matrix lT;
        ToMatrix4x4( lT );

        double lDet = lT.M4x4_Det();

        if( !LzMath::ToolBox::IsZero( lDet - 1, pTol ) )
            LzLogException("",  "Matrix does not define a direct rigid rotation! (det= " << lDet << ")" );
    }
#else
    LzLogException("", "*** TODO RigidTr3D::FromOldString( const std::string & pStdStr, bool pCheckRotDet/*=true*/, double pTol/*=-1*/ )");
#endif
}

//================================================================================
void RigidTr3D::FromNewString( const std::string & pStdStr, bool pCheckRotDet/*=true*/, double pTol/*=-1*/ )
{
	// Values
    double lR[9];
    double lT[3];

    // STDLIB
    vector<string> lTokens;
    LzServices::SplitString( pStdStr, {" ", "\t", ","}, lTokens );

    // Check
    if( lTokens.size() != 12 )
        LzLogException("", "Found "<<lTokens.size()<<" token(s), expected 12!")

    // Try reading
    _LzLogTry
    {
        // Rotation
        lR[0] = LzServices::StrToDouble( lTokens[ 0] );
        lR[1] = LzServices::StrToDouble( lTokens[ 1] );
        lR[2] = LzServices::StrToDouble( lTokens[ 2] );

        lT[0] = LzServices::StrToDouble( lTokens[ 3] );

        lR[3] = LzServices::StrToDouble( lTokens[ 4] );
        lR[4] = LzServices::StrToDouble( lTokens[ 5] );
        lR[5] = LzServices::StrToDouble( lTokens[ 6] );

        lT[1] = LzServices::StrToDouble( lTokens[ 7] );

        lR[6] = LzServices::StrToDouble( lTokens[ 8] );
        lR[7] = LzServices::StrToDouble( lTokens[ 9] );
        lR[8] = LzServices::StrToDouble( lTokens[10] );

        lT[2] = LzServices::StrToDouble( lTokens[11] );
    }
    _LzLogCatchAndThrow("Cannot parse rigid transform3D! Syntax error in string: "<<pStdStr<<".")

 	// Commit
	for( int i=0 ; i<9 ; i++ )
		mR[i] = lR[i];

	for( int i=0 ; i<3 ; i++ )
		mT[i] = lT[i];

	// Check determinant?
    if( pCheckRotDet )
    {
        const double lDet = R_Determinant();

        if( !LzMath::ToolBox::IsZero( lDet - 1, pTol ) )
            LzLogException("", "Matrix does not define a direct rigid rotation! (det= " << std::setprecision(15) << lDet << ")" )
    }
}
}
