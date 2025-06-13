#include "Plane3D.h"
#include "Vector3D.h"
#include "Line3D.h"
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>
#if defined(USING_QT) && !defined(NO_OPENGL)
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/gl.h>
#endif
//
#include <cmath>


namespace LzGeom
{
#pragma region "Construction / destruction"
//================================================================================
Plane3D::Plane3D()
{
    mP[0] = 1;
    mP[1] = 0;
    mP[2] = 0;
    mP[3] = 0;

    Normalize();
}

//================================================================================
Plane3D::Plane3D( std::function<void(Plane3D * pThis)> pInit )
{
    // Set default plane for pInit to work on
    mP[0] = 1;
    mP[1] = 0;
    mP[2] = 0;
    mP[3] = 0;

    Normalize();

    // Init
    pInit( this );
}

//================================================================================
Plane3D::Plane3D( double pA, double pB, double pC, double pD )
{
    mP[0] = pA;
    mP[1] = pB;
    mP[2] = pC;
    mP[3] = pD;

    Normalize();
}

//================================================================================
Plane3D::Plane3D( const double * pP )
{
    mP[0] = pP[0];
    mP[1] = pP[1];
    mP[2] = pP[2];
    mP[3] = pP[3];

    Normalize();
}

//================================================================================
Plane3D::Plane3D( const Point3D & pPoint, const Vector3D & pNormal )
{
    mP[0] = pNormal.X();
    mP[1] = pNormal.Y();
    mP[2] = pNormal.Z();
    mP[3] = -( mP[0] * pPoint.X() + mP[1] * pPoint.Y() + mP[2] * pPoint.Z() );

    Normalize();
}

//================================================================================
Plane3D::Plane3D( const Point3D & pPoint, const Vector3D & pX, const Vector3D & pY )
{
    Vector3D lNormal = pX ^ pY;

    mP[0] = lNormal.X();
    mP[1] = lNormal.Y();
    mP[2] = lNormal.Z();
    mP[3] = -( mP[0] * pPoint.X() + mP[1] * pPoint.Y() + mP[2] * pPoint.Z() );

    Normalize();
}

//================================================================================
Plane3D::Plane3D( const Point3D & pPoint1, const Point3D & pPoint2, const Point3D & pPoint3 )
{
    Vector3D lNormal = ( pPoint2 - pPoint1 ) ^ ( pPoint3 - pPoint1 );

    mP[0] = lNormal.X();
    mP[1] = lNormal.Y();
    mP[2] = lNormal.Z();
    mP[3] = -(mP[0]*pPoint1.X() + mP[1]*pPoint1.Y() + mP[2]*pPoint1.Z());

    Normalize();
}

//================================================================================
Plane3D::Plane3D( const vector<Point3D> & pPointCloud )
{
    FromPointCloud( pPointCloud );
}

//================================================================================
Plane3D::Plane3D( const Vector<Point3D> & pPointCloud )
{
    vector<Point3D> lPC( pPointCloud.Size() );
    for( unsigned int p = 0 ; p < pPointCloud.Size() ; p++ )
        lPC[p] = pPointCloud[p];

    FromPointCloud( lPC );
}

//================================================================================
Plane3D::Plane3D( const List<Point3D> & pPointCloud )
{
    vector<Point3D> lCloud;
    pPointCloud.ToVector( lCloud );

    FromPointCloud( lCloud );
}

//================================================================================
void Plane3D::SetPlaneCoeffs( double pA, double pB, double pC, double pD )
{
    mP[ 0 ] = pA;
    mP[ 1 ] = pB;
    mP[ 2 ] = pC;
    mP[ 3 ] = pD;

    Normalize();
}

//================================================================================
void Plane3D::FromPointCloud( const vector<Point3D> & pPointCloud )
{
    // Check: at least 3 points needed, however not guarantee that these points
    //        will indeed define a plane. The PCA will always return something
    //        anyway.
    if( pPointCloud.size() < 3 )
        LzLogException("",  "Unable to compute least-squares plane using only "<<pPointCloud.size()<<" point(s)!")

    Point3D lMean;
    Vector3D lV[3];
    LzMath::ToolBox::PCA( pPointCloud, lMean, lV );

    const Vector3D & lN = lV[2];
    mP[0] = lN.X();
    mP[1] = lN.Y();
    mP[2] = lN.Z();
    mP[3] = -( mP[0] * lMean.X() + mP[1] * lMean.Y() + mP[2] * lMean.Z() );

    Normalize();

//************** essayer de minimiser sous la contrainte non lineaire norm(n) = 1

#if 0
    // Fill coordinates matrix
    Matrix X( m, 4 );
    {
        unsigned int i = 0;
        vector<Point3D>::const_iterator iPt;
        for( iPt = pPointCloud.begin() ; iPt != pPointCloud.end() ; iPt++ )
        {
            X.Elt( i, 0 ) = iPt->X();
            X.Elt( i, 1 ) = iPt->Y();
            X.Elt( i, 2 ) = iPt->Z();
            X.Elt( i, 3 ) = 1;

            i++;
        }
    }

    Matrix XtX;
    X.SymmetricalMult_MtM( XtX );


//*****************************************************************
    {
        Matrix M( 5, 5 );
        M.SubMatrixFrom( 0, 0, XtX );

        // Contrainte de Lagrange
        M.Elt( 0, 4 ) = 1;
        M.Elt( 1, 4 ) = 0;
        M.Elt( 2, 4 ) = 0;
        M.Elt( 3, 4 ) = 0;

        M.Elt( 4, 0 ) = 1;
        M.Elt( 4, 1 ) = 0;
        M.Elt( 4, 2 ) = 0;
        M.Elt( 4, 3 ) = 0;

        M.Elt( 4, 4 ) = 0;

/////// faire 3 optimisations sous L avec chacune des coords non nulle et retenir solution avec moindre residu

        M.Log( "M" );

        Matrix TOTO;
        M.CF_InverseTo( TOTO );

        Matrix Y( 5 );
        Y.LoadValue( 0 );
        Y.Elt( 4 ) = 1;

        Matrix V;
        TOTO.Mult( Y, V );

        V.Log( "VX" );

        // Commit
        for( unsigned int i = 0 ; i < 4 ; i++ )
            mV[i] = V.Elt( i );

        Normalize();
        zboub !

    }
//*****************************************************************
//{
//  Matrix M( 5, 5 );
//  M.SubMatrixFrom( 0, 0, XtX );
//
//  // Contrainte de Lagrange
//  M.Elt(0,4) = 0;
//  M.Elt(1,4) = 1;
//  M.Elt(2,4) = 0;
//  M.Elt(3,4) = 0;
//
//  M.Elt(4,0) = 0;
//  M.Elt(4,1) = 1;
//  M.Elt(4,2) = 0;
//  M.Elt(4,3) = 0;
//
//  M.Elt(4,4) = 0;
//
//      /////// faire 3 optimisations sous L avec chacune des coords non nulle et retenir solution avec moindre residu
//
//      M.Log("M");
//
//  Matrix TOTO;
//  M.CF_InverseTo( TOTO );
//
//  Matrix Y( 5 );
//  Y.LoadValue( 0 );
//  Y.Elt(4) = 1;
//
//  Matrix V;
//  TOTO.Mult( Y, V );
//
//  V.Log("VY");
//
//  // Commit
//  for( unsigned int i=0 ; i<4 ; i++ )
//      mV[i] = V.Elt(i);
//
//  Normalize();
//
//  LzLogM("", "PLANE Y "+ToString());
//}
//*****************************************************************
//{
//  Matrix M( 5, 5 );
//  M.SubMatrixFrom( 0, 0, XtX );
//
//  // Contrainte de Lagrange
//  M.Elt(0,4) = 0;
//  M.Elt(1,4) = 0;
//  M.Elt(2,4) = 1;
//  M.Elt(3,4) = 0;
//
//  M.Elt(4,0) = 0;
//  M.Elt(4,1) = 0;
//  M.Elt(4,2) = 1;
//  M.Elt(4,3) = 0;
//
//  M.Elt(4,4) = 0;
//
//      /////// faire 3 optimisations sous L avec chacune des coords non nulle et retenir solution avec moindre residu
//
//      M.Log("M");
//
//  Matrix TOTO;
//  M.CF_InverseTo( TOTO );
//
//  Matrix Y( 5 );
//  Y.LoadValue( 0 );
//  Y.Elt(4) = 1;
//
//  Matrix V;
//  TOTO.Mult( Y, V );
//
//  V.Log("VZ");
//
//  // Commit
//  for( unsigned int i=0 ; i<4 ; i++ )
//      mV[i] = V.Elt(i);
//
//  Normalize();
//
//  LzLogM("", "PLANE Z "+ToString());
//}
#endif
}

//================================================================================
void Plane3D::Normalize()
{
    double lN = sqrt( mP[0]*mP[0] + mP[1]*mP[1] + mP[2]*mP[2] );

    if( LzMath::ToolBox::IsZero( lN ) )
        LzLogException("",  "Invalid plane equation!" );

    mP[0] /= lN;
    mP[1] /= lN;
    mP[2] /= lN;
    mP[3] /= lN;

    mNormal = Vector3D( mP[0], mP[1], mP[2] );

    mNormal.BuildOrthonormalBase( mX, mY, mNormal );
}
#pragma endregion


#pragma region "Computing"
//================================================================================
void Plane3D::FlipNormal()
{
    for( int i = 0 ; i < 4 ; i++ )
        mP[i] *= -1;

    Normalize();
}

//================================================================================
double Plane3D::SignedDistanceTo( const Point3D & pPoint ) const
{
    return mP[0] * pPoint.X() + mP[1] * pPoint.Y() + mP[2] * pPoint.Z() + mP[3];
}

//================================================================================
Point3D Plane3D::Symmetrical( const Point3D & pPoint ) const
{
    return pPoint + 2 * ( Projection( pPoint ) - pPoint );
}

//================================================================================
Vector3D Plane3D::Symmetrical( const Vector3D & pVector ) const
{
    return pVector - 2 * ( pVector * mNormal ) * mNormal;
}

//================================================================================
Point3D Plane3D::Projection( const Point3D & pPoint ) const
{
    // The plane is normalized
    double d = SignedDistanceTo( pPoint );
    return pPoint - d * mNormal;
}

//================================================================================
Vector3D Plane3D::Projection( const Vector3D & pV ) const
{
    const Point3D A;
    const Point3D B = A + pV;

    return Projection( B ) - Projection( A );
}

//================================================================================
Point3D Plane3D::Projection( const Point3D & pPoint, const Vector3D & pProjDir ) const
{
    double lScal = mNormal * pProjDir;

    if( LzMath::ToolBox::IsZero( lScal ) )
        LzLogException("",  "Cannot project on plane along a parallel direction!" );

    double k = SignedDistanceTo( pPoint ) / lScal;
    return pPoint - k * pProjDir;
}

//================================================================================
bool Plane3D::Projection_NoExc( const Point3D & pPoint, const Vector3D & pProjDir, Point3D & pTo ) const
{
    double lScal = mNormal * pProjDir;

	// Check
    if( LzMath::ToolBox::IsZero( lScal ) )
        return false;

	// Compute projection
    double k = SignedDistanceTo( pPoint ) / lScal;
    pTo = pPoint - k * pProjDir;

	return true;
}

//================================================================================
Line3D Plane3D::Intersection( const Plane3D & pOther ) const
{
    Line3D lLine;

    // Compute line direction
    lLine.mDir = mNormal ^ pOther.Normal();

    // Find a root point starting from origin
    lLine.mRoot = pOther.Projection( Projection( Point3D( 0, 0, 0 ) ), lLine.mDir ^ mNormal );
    lLine.Normalize();
    return lLine;
}

//================================================================================
bool Plane3D::Intersection_NoExc( const Plane3D & pOther, Line3D & pTo ) const
{
    // Compute line direction
    pTo.mDir = mNormal ^ pOther.Normal();

	// Check
	double lDirNorm = pTo.mDir.Norm();
    if( LzMath::ToolBox::IsZero(lDirNorm) )
		return false;

    // Find a root point starting from origin
	if( !pOther.Projection_NoExc(Projection(Point3D(0, 0, 0)), pTo.mDir^mNormal, pTo.mRoot) )
		return false;

	pTo.mDir /= lDirNorm;

    return true;
}

//================================================================================
Point3D Plane3D::Intersection( const Line3D & pLine ) const
{
    return Projection( pLine.mRoot, pLine.mDir );
}

//================================================================================
bool Plane3D::Intersection_NoExc( const Line3D & pLine, Point3D & pTo ) const
{
	if( Projection_NoExc(pLine.mRoot, pLine.mDir, pTo) )
		return true;
	else
		return false;
}

//================================================================================
void Plane3D::SetPositiveSide( const Vector3D & pToPos )
{
    if( mNormal * pToPos < 0 )
        *this = Plane3D( -mP[0], -mP[1], -mP[2], -mP[3] );
}

//================================================================================
void Plane3D::SetPositiveSide( const Point3D & pInPos )
{
    if( SignedDistanceTo( pInPos ) < 0 )
        *this = Plane3D( -mP[0], -mP[1], -mP[2], -mP[3] );
}

//================================================================================
void Plane3D::SetNegativeSide( const Point3D & pInNeg )
{
    if( SignedDistanceTo( pInNeg ) > 0 )
        *this = Plane3D( -mP[0], -mP[1], -mP[2], -mP[3] );
}

//================================================================================
const Plane3D & Plane3D::operator*=( const RigidTr3D & pT )
{
    *this = pT * ( *this );
    return *this;
}

//================================================================================
const Plane3D & Plane3D::operator*=( double pScale )
{
    mP[3] *= pScale;
    return *this;
}

//================================================================================
const Plane3D & Plane3D::operator/=( double pScale )
{
    mP[3] /= pScale;
    return *this;
}

//================================================================================
Plane3D Plane3D::operator+( const Vector3D & pV ) const
{
    return Plane3D( mP[0], mP[1], mP[2], mP[3] - pV * Normal() );
}

//================================================================================
Plane3D Plane3D::operator-( const Vector3D & pV ) const
{
    return Plane3D( mP[0], mP[1], mP[2], mP[3] + pV * Normal() );
}

//================================================================================
const Plane3D & Plane3D::operator+=( const Vector3D & pV )
{
    *this = *this + pV;
    return *this;
}

//================================================================================
const Plane3D & Plane3D::operator-=( const Vector3D & pV )
{
    *this = *this - pV;
    return *this;
}

//================================================================================
bool Plane3D::operator!=( const Plane3D & pOther ) const
{
#if 1
    // Assumption: both planes are normalized i.e. A^2 + B^2 + C^2 = 1

    Vector3D u(        mP[0],        mP[1],        mP[2] );
    Vector3D v( pOther.mP[0], pOther.mP[1], pOther.mP[2] );

    double s = u * v;

    if( s <= 0 )
        return true;

    double k = 1.0 / s;
//double k = u.SquaredNorm() / s;

    return !LzMath::ToolBox::IsZero( mP[0] - k * pOther.mP[0] )
           || !LzMath::ToolBox::IsZero( mP[1] - k * pOther.mP[1] )
           || !LzMath::ToolBox::IsZero( mP[2] - k * pOther.mP[2] )
           || !LzMath::ToolBox::IsZero( mP[3] - k * pOther.mP[3] );
#else
    // Different normals?
    if( !LzMath::ToolBox::IsZero( ( mNormal - pOther.mNormal ).Norm() ) )
        return true;

    // Same normals: compare 2 projected points
    const Point3D O1 = Projection( Point3D( 0, 0, 0 ) );
    const Point3D O2 = pOther.Projection( Point3D( 0, 0, 0 ) );

    // Same projections?
    return !LzMath::ToolBox::IsZero( O1.DistanceTo( O2 ) );
#endif
}

//================================================================================
Matrix Plane3D::GetSymmetryMat4x4() const
{
    // Assumption: this plane is normalized i.e. A^2 + B^2 + C^2 = 1

    // Fill matrix row by row
    Matrix lSym( 4, 4 );

    lSym.Elt( 0, 0 ) = 1 - 2 * mP[0] * mP[0];
    lSym.Elt( 0, 1 ) =   - 2 * mP[0] * mP[1];
    lSym.Elt( 0, 2 ) =   - 2 * mP[0] * mP[2];
    lSym.Elt( 0, 3 ) =   - 2 * mP[0] * mP[3];

    lSym.Elt( 1, 0 ) =   - 2 * mP[1] * mP[0];
    lSym.Elt( 1, 1 ) = 1 - 2 * mP[1] * mP[1];
    lSym.Elt( 1, 2 ) =   - 2 * mP[1] * mP[2];
    lSym.Elt( 1, 3 ) =   - 2 * mP[1] * mP[3];

    lSym.Elt( 2, 0 ) =   - 2 * mP[2] * mP[0];
    lSym.Elt( 2, 1 ) =   - 2 * mP[2] * mP[1];
    lSym.Elt( 2, 2 ) = 1 - 2 * mP[2] * mP[2];
    lSym.Elt( 2, 3 ) =   - 2 * mP[2] * mP[3];

    lSym.Elt( 3, 0 ) = 0;
    lSym.Elt( 3, 1 ) = 0;
    lSym.Elt( 3, 2 ) = 0;
    lSym.Elt( 3, 3 ) = 1;

    return lSym;
}

//================================================================================
bool Plane3D::IsSameAs( const Plane3D & pOther, double pTol ) const
{
    // Check
    for( int i=0 ; i<4 ; i++ )
    {
        if( fabs(mP[i] - pOther.mP[i]) > pTol )
            return false;
    }

    return true;
}
#pragma endregion


#pragma region "Drawing"
#if !defined(NO_OPENGL)
//================================================================================
void Plane3D::Draw( const Point3D & pCenter, double pDimX, double pDimY/*=-1*/ ) const
{
    // // Set defaults
    // if( pDimY == -1 )
    //     pDimY = pDimX;

    // vector<Point3D> lQuad;
    // DrawToArray( pCenter, pDimX, pDimY, lQuad );

    // glBegin( GL_QUADS );
    // {
    //     glNormal3dv( mNormal.mV );
    //     glVertex3dv( lQuad[0].mV );
    //     glVertex3dv( lQuad[1].mV );
    //     glVertex3dv( lQuad[2].mV );
    //     glVertex3dv( lQuad[3].mV );
    // }
    // glEnd();
}

//================================================================================
void Plane3D::Draw( const Point3D & pCenter, double pDimX, double pDimY, unsigned int pStripsX, unsigned int pStripsY ) const
{
    // if( !pStripsX || !pStripsY )
    //     return;

    // const double lStepX = pDimX / pStripsX;
    // const double lStepY = pDimY / pStripsY;

    // Point3D O = Projection( pCenter ) - ( pDimX / 2 ) * mX - ( pDimY / 2 ) * mY;

    // glNormal3dv( mNormal.mV );

    // for( unsigned int x = 0 ; x < pStripsX ; x++ )
    // {
    //     glBegin( GL_QUAD_STRIP );
    //     for( unsigned int y = 0 ; y < pStripsY + 1 ; y++ )
    //     {
    //         glVertex3dv( ( O + ( x + 0 )*lStepX * mX + y * lStepY * mY ).mV );
    //         glVertex3dv( ( O + ( x + 1 )*lStepX * mX + y * lStepY * mY ).mV );
    //     }
    //     glEnd();
    // }
}
#endif

//================================================================================
void Plane3D::DrawToArray( const Point3D & pCenter, double pDimX, double pDimY, vector<Point3D> & pQuad ) const
{
    // Build referential
    Point3D C = Projection( pCenter );
#if 1
    double lHalfDimX = pDimX / 2;
    double lHalfDimY = pDimY / 2;

    Vector3D XpY = lHalfDimX * mX + lHalfDimY * mY;
    Vector3D XmY = lHalfDimX * mX - lHalfDimY * mY;

    pQuad.resize( 4 );

    pQuad[0] = C - XpY;
    pQuad[1] = C + XmY;
    pQuad[2] = C + XpY;
    pQuad[3] = C - XmY;
#else
    double lHalfDimX = pDimX / 2;
    double lHalfDimY = pDimY / 2;

    pQuad.resize( 4 );

    pQuad[0] = C - lHalfDimX * mX - lHalfDimY * mY;
    pQuad[1] = C + lHalfDimX * mX - lHalfDimY * mY;
    pQuad[2] = C + lHalfDimX * mX + lHalfDimY * mY;
    pQuad[3] = C - lHalfDimX * mX + lHalfDimY * mY;
#endif
}
#pragma endregion


#pragma region "Debug"
//================================================================================
std::string Plane3D::ToString() const
{
	std::stringstream lStrStr;
	lStrStr << "(" << A() << ")X + (" << B() << ")Y + (" << C() << ")Z + (" << D() << ") = 0";
    return lStrStr.str();
}

//================================================================================
std::string Plane3D::ToString_ABCD() const
{
	std::stringstream lStrStr;
	lStrStr << A() << " " << B() << " " << C() << " " << D();
    return lStrStr.str();
}
#pragma endregion

}
