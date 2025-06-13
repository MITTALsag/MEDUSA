#include "Matrix.h"
//#include "Polynom.h"
//#include "PolyMatrix.h"
#include "ToolBox.h"
#include <LzServices/LzLog.h>
#include <iomanip>


#ifdef USE_CLI
#include <math.h>
#endif


// @TODO LzMath::Matrix : templeter cette classe avec specialisations
// pour Polynom, Map, etc ?


namespace LzMath
{
#ifdef USE_CLI
using System::String;
#endif
using std::vector;


#pragma region "Construction & destruction"
//================================================================================
Matrix::Matrix( unsigned int pH/*=0*/, unsigned int pW/*=1*/ ) : mH( 0 ), mW( 0 )
{
    SetDims( pH, pW );
}

//================================================================================
Matrix::Matrix( unsigned int pH, unsigned int pW, float const * pRowByRowElts ) : mH( 0 ), mW( 0 )
{
    SetDims( pH, pW );

    // Set entries
    unsigned int e = 0;
    for( unsigned int i = 0 ; i < pH ; i++ )
    for( unsigned int j = 0 ; j < pW ; j++ )
        Elt( i, j ) = pRowByRowElts[e++];
}

//================================================================================
Matrix::Matrix( unsigned int pH, unsigned int pW, double const * pRowByRowElts ) : mH( 0 ), mW( 0 )
{
    SetDims( pH, pW );

    // Set entries
    unsigned int e = 0;
    for( unsigned int i = 0 ; i < pH ; i++ )
    for( unsigned int j = 0 ; j < pW ; j++ )
        Elt( i, j ) = pRowByRowElts[e++];
}

//================================================================================
Matrix::Matrix( unsigned int pH, unsigned int pW, const string & pRowByRowElts ) : mH( 0 ), mW( 0 )
{
#ifdef USE_CLI
    SetDims( pH, pW );

	String ^lCLI_RowByRowElts = gcnew String( pRowByRowElts.c_str() );

    // Split
    cli::array<String^> ^lItems = lCLI_RowByRowElts->Split( gcnew cli::array<wchar_t> { ' ', '\t' }, System::StringSplitOptions::RemoveEmptyEntries );

    // Check
    if( lItems->Length != pH * pW )
        LzLogException("", "Wrong number of initializers! Expected "<<(pH * pW)<<", found "<<lItems->Length<<" in string '"<<pRowByRowElts<<"'." );

    // Set entries
    unsigned int e = 0;
    for( unsigned int i = 0 ; i < pH ; i++ )
    for( unsigned int j = 0 ; j < pW ; j++ )
        Elt( i, j ) = TxServices::StrToDouble( lItems[e++] );
#else
    SetDims( pH, pW );

    // Split
    vector<string> lItems;
    LzServices::SplitString( pRowByRowElts, {" ", "\t"}, lItems );

    // Check
    if( lItems.size() != pH * pW )
        LzLogException("", "Wrong number of initializers! Expected "<<(pH * pW)<<", found "<<lItems.size()<<" in string '"<<pRowByRowElts<<"'." );

    // Set entries
    unsigned int e = 0;
    for( unsigned int i = 0 ; i < pH ; i++ )
    for( unsigned int j = 0 ; j < pW ; j++ )
        Elt( i, j ) = LzServices::StrToDouble( lItems[e++] );
#endif
}

//================================================================================
Matrix::Matrix( unsigned int pH, unsigned int pW, const std::initializer_list<double> & pRowByRowElts ) : mH( 0 ), mW( 0 )
{
    SetDims( pH, pW );

    // Check
    if( pRowByRowElts.size() != pH * pW )
        LzLogException("", "Wrong number of initializers! Expected "<<(pH * pW)<<", found "<<pRowByRowElts.size()<<"." );

    // Set entries
    unsigned int e = 0;
    for( unsigned int i = 0 ; i < pH ; i++ )
    for( unsigned int j = 0 ; j < pW ; j++ )
        Elt( i, j ) = pRowByRowElts.begin()[ e++ ];
}

//================================================================================
Matrix::Matrix( const Matrix & pM ) : mH( 0 ), mW( 0 )
{
    *this = pM;
}

//================================================================================
Matrix::~Matrix()
{
    Free();
}
#pragma endregion


#pragma region "Init"
//================================================================================
void Matrix::SetDims( unsigned int pH, unsigned int pW/*=1*/ )
{
    unsigned int lNewSize = pH * pW;

    // Need to reallocate space?
    if( lNewSize != mH * mW )
    {
        Free();

        // Non empty matrix?
        if( lNewSize )
            mMat = new double [lNewSize];
    }

    // Store new dimensions
    mW = pW;
    mH = pH;
}

//================================================================================
void Matrix::Free()
{
    if( mH * mW )
        delete mMat;

    mH = 0;
    mW = 0;
}

Matrix Matrix::Identity( unsigned int pH, unsigned int pW )
{
    Matrix lMat( pH, pW );
    lMat.LoadIdentity();

    return lMat;
}

#pragma endregion


#pragma region "Get & set"
//================================================================================
double & Matrix::Elt( unsigned int pI, unsigned int pJ/*=0*/ )
{
    // Check
    if( pI >= mH || pJ >= mW )
        LzLogException("","Accessing invalid element (i="<<pI<<", j="<<pJ<<")! Matrix is "<<mH<<" x "<<mW<<".");

//// Row-major
//return mMat[ pI*mW + pJ ];

// Col-major
return mMat[ pI + pJ*mH ];
}

//================================================================================
double Matrix::Elt( unsigned int pI, unsigned int pJ/*=0*/ ) const
{
    // Check
    if( pI >= mH || pJ >= mW )
        LzLogException("","Accessing invalid element (i="<<pI<<", j="<<pJ<<")! Matrix is "<<mH<<" x "<<mW<<".");

//// Row-major
//return mMat[ pI*mW + pJ ];

// Col-major
return mMat[ pI + pJ*mH ];
}

//================================================================================
const double * Matrix::ColMajorFrom( unsigned int pJ ) const
{
    // Check
    if( pJ >= mW )
        LzLogException("","Accessing invalid column (j="<<pJ<<")! Matrix is "<<mH<<" x "<<mW<<".");

	// Col-major
	return mMat + pJ*mH;
}
#pragma endregion



//********************************************* sore made

//********************************************* replace method separators !!

//================================================================================
//const double & CMatrix::Elt( int p_I, int p_J ) const
//{
//#ifdef _DEBUG
//  // Ok ?
//  if( p_I>=m_H || p_I<0
//   || p_J>=m_W || p_J<0 )
//  {
//      LOGERR("[CMatrix::Elt] : Accessing invalid element (i=%d,j=%d) !",p_I,p_J);
//      return ERR_ELT;
//  }
//#endif

//  return m_pMat[p_J][p_I];
//}

//================================================================================
//BOOL CMatrix::SubMatrixTo( int p_FirstRow, int p_LastRow, int p_FirstCol, int p_LastCol, CMatrix & p_To )
//{
//  // Ok ?
//  if( p_FirstRow>p_LastRow
//   || p_FirstCol>p_LastCol
//   || p_FirstRow<0
//   || p_FirstCol<0
//   || p_LastRow>=m_H
//   || p_LastCol>=m_W )
//  {
//      LOGERR("[CMatrix::SubMatrixTo] : Invalid range !");
//      return FALSE;
//  }

//  // Fill new matrix
//  p_To.SetDims( p_LastRow-p_FirstRow+1, p_LastCol-p_FirstCol+1 );
//  for( int i=0 ; i<p_To.Rows() ; i++ )
//  for( int j=0 ; j<p_To.Cols() ; j++ )
//  {
//      p_To.Elt(i,j) = Elt( p_FirstRow+i, p_FirstCol+j );
//  }

//  return TRUE;
//}

//================================================================================
//BOOL CMatrix::SubMatrixTo( int p_FirstRow, int p_LastRow, int p_Col, CVector & p_To )
//{
//  // Ok ?
//  if( p_FirstRow>p_LastRow
//   || p_FirstRow<0
//   || p_LastRow>=m_H
//   || p_Col<0
//   || p_Col>=m_W )
//  {
//      LOGERR("[CMatrix::SubMatrixTo] : Invalid range !");
//      return FALSE;
//  }

//  // Fill new matrix
//  p_To.SetDim( p_LastRow-p_FirstRow+1 );
//  for( int i=0 ; i<p_To.Rows() ; i++ )
//  {
//      p_To.Elt(i) = Elt( p_FirstRow+i, p_Col );
//  }

//  return TRUE;
//}

//================================================================================
//BOOL CMatrix::SubMatrixTo( int p_Col, CVector & p_To )
//{
//  return SubMatrixTo( 0, Rows()-1, p_Col, p_To );
//}

//================================================================================
//const double * CMatrix::GetColBuffer( int p_J ) const
//{
//  if( p_J<0 || p_J>m_W )
//  {
//      LOGERR("[CMatrix::GetColBuffer] : Invalid col index (%d) !",p_J);
//      return NULL;
//  }

//  return m_pMat[p_J];
//}


#pragma region "Operations"
//================================================================================
void Matrix::LoadIdentity()
{
    // Check
    if( mH != mW )
        LzLogException("","Matrix is not square! ("<<mH <<"x"<<mW<<")");

	// Load diagonal
    for( unsigned int i = 0 ; i < mH ; i++ )
    for( unsigned int j = 0 ; j < mW ; j++ )
        Elt( i, j ) = ( i == j ? 1 : 0 );
}

//================================================================================
void Matrix::LoadDiagonal( const Matrix & pDiagInCol, unsigned int pColIdx/*=0*/ )
{
    // Check
    if( mH != mW )
        LzLogException("","Matrix is not square! ("<<mH<<"x"<<mW<<")");

    // Check
    if( mH != pDiagInCol.Rows() )
        LzLogException("","Rows count mismatch! ("<<mH<<" != "<<pDiagInCol.Rows()<<")");

    // Check
    if( pColIdx >= pDiagInCol.Cols() )
        LzLogException("","Invalid column index! ("<<pColIdx<<" >= "<<pDiagInCol.Cols()<<")");

	// Load diagonal
    for( unsigned int i = 0 ; i < mH ; i++ )
    for( unsigned int j = 0 ; j < mW ; j++ )
        Elt( i, j ) = ( i == j ? pDiagInCol.Elt(i,pColIdx) : 0 );
}

//================================================================================
void Matrix::LoadValue( double pVal )
{
    for( unsigned int i = 0 ; i < mH ; i++ )
    for( unsigned int j = 0 ; j < mW ; j++ )
        Elt( i, j ) = pVal;
}

//================================================================================
void Matrix::Mult( const Matrix & pA, Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() || Cols() != pA.Rows() || !pA.Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Check
    if( &pRes == &pA || &pRes == this  )
        LzLogException("", "Illegal operation!" );

    pRes.SetDims( Rows(), pA.Cols() );

    for( unsigned int i = 0 ; i < pRes.Rows() ; i++ )
    for( unsigned int j = 0 ; j < pRes.Cols() ; j++ )
    {
        pRes.Elt( i, j ) = 0;

        for( unsigned int k = 0 ; k < Cols() ; k++ )
            pRes.Elt( i, j ) += Elt( i, k ) * pA.Elt( k, j );
    }
}

//================================================================================
void Matrix::SymmetricalMult( const Matrix & pA, Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() || Cols() != pA.Rows() || !pA.Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Check
    if( &pRes == &pA || &pRes == this  )
        LzLogException("", "Illegal operation!" );

    // Check
    if( Rows() != pA.Cols() )
        LzLogException("", "Invalid dimensions!" );

    pRes.SetDims( Rows(), pA.Cols() );

    for( unsigned int i = 0 ; i < pRes.Rows() ; i++ )
    for( unsigned int j = i ; j < pRes.Cols() ; j++ )
    {
        // Upper triangle
        pRes.Elt( i, j ) = 0;

        for( unsigned int k = 0 ; k < Cols() ; k++ )
            pRes.Elt( i, j ) += Elt( i, k ) * pA.Elt( k, j );

        // Lower triangle = symmetrical
        pRes.Elt( j, i ) = pRes.Elt( i, j );
    }
}

//================================================================================
void Matrix::SymmetricalMult_MtM( Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!" );

    if( &pRes == this  )
        LzLogException("", "Illegal operation!" );

    pRes.SetDims( Cols(), Cols() );

    for( unsigned int i = 0 ; i < pRes.Rows() ; i++ )
    for( unsigned int j = i ; j < pRes.Cols() ; j++ )
    {
        // Upper triangle
        pRes.Elt( i, j ) = 0;

        for( unsigned int k = 0 ; k < Rows() ; k++ )
            pRes.Elt( i, j ) += Elt( k, i ) * Elt( k, j );

        // Lower triangle = symmetrical
        pRes.Elt( j, i ) = pRes.Elt( i, j );
    }
}

//================================================================================
void Matrix::SymmetricalMult_MtM_upper_tri_only( Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!" );

    if( &pRes == this  )
        LzLogException("", "Illegal operation!" );

    pRes.SetDims( Cols(), Cols() );

    for( unsigned int i = 0 ; i < pRes.Rows() ; i++ )
    for( unsigned int j = i ; j < pRes.Cols() ; j++ )
    {
        // Upper triangle
        pRes.Elt( i, j ) = 0;

        for( unsigned int k = 0 ; k < Rows() ; k++ )
            pRes.Elt( i, j ) += Elt( k, i ) * Elt( k, j );
    }
}

//================================================================================
//void Matrix::Mult( const PolyMatrix & pA, PolyMatrix & pRes ) const
//{
//    // Check
//    if( !Rows() || !Cols() || Cols() != pA.Rows() || !pA.Cols() )
//        LzLogException("", "Invalid dimensions!" );
//
//    pRes.SetDims( Rows(), pA.Cols() );
//
//    for( unsigned int i = 0 ; i < pRes.Rows() ; i++ )
//    for( unsigned int j = 0 ; j < pRes.Cols() ; j++ )
//    {
//        pRes.Elt( i, j ) = 0;
//
//        for( unsigned int k = 0 ; k < Cols() ; k++ )
//            pRes.Elt( i, j ) = pRes.Elt( i, j ) + Elt( i, k ) * pA.Elt( k, j );
//    }
//}

//================================================================================
void Matrix::Mult( double pA, Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!" );

    pRes.SetDims( Rows(), Cols() );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pRes.Elt( i, j ) = pA * Elt( i, j );
}

//================================================================================
void Matrix::Mult( double pA )
{
    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        Elt( i, j ) *= pA;
}

//================================================================================
void Matrix::AddMult( double pA, const Matrix & pB, Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() || Rows() != pB.Rows() || Cols() != pB.Cols() )
        LzLogException("", "Invalid dimensions!" );

    pRes.SetDims( Rows(), Cols() );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pRes.Elt( i, j ) = Elt( i, j ) + pA * pB.Elt( i, j );
}

//================================================================================
void Matrix::Add( const Matrix & pA, Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() || Rows() != pA.Rows() || Cols() != pA.Cols() )
        LzLogException("", "Invalid dimensions!" );

    pRes.SetDims( Rows(), Cols() );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pRes.Elt( i, j ) = Elt( i, j ) + pA.Elt( i, j );
}

//================================================================================
void Matrix::Add( const Matrix & pA )
{
    // Check
    if( !Rows() || !Cols() || Rows() != pA.Rows() || Cols() != pA.Cols() )
        LzLogException("", "Invalid dimensions!" );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        Elt( i, j ) += pA.Elt( i, j );
}

//================================================================================
void Matrix::Sub( const Matrix & pA, Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() || Rows() != pA.Rows() || Cols() != pA.Cols() )
        LzLogException("", "Invalid dimensions!" );

    pRes.SetDims( Rows(), Cols() );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pRes.Elt( i, j ) = Elt( i, j ) - pA.Elt( i, j );
}

//================================================================================
//BOOL CMatrix::Pow( int p_Pow, CMatrix & p_Res ) const
//{
//  if( Cols() != Rows() )
//  {
//      LOGERR("[CMatrix::Pow] : Matrix not square !");
//      return FALSE;
//  }

//  CMatrix l_Tmp;
//  if( !p_Res.SetDims(Rows(),Cols())
//   || !l_Tmp.SetDims(Rows(),Cols()) )
//  {
//      LOGERR("[CMatrix::Pow] : Unable to set matrices !");
//      return FALSE;
//  }

//  p_Res.LoadIdentity();
//
//  while( p_Pow )
//  {
//      Mult(p_Res,l_Tmp);
//      p_Res = l_Tmp;
//      p_Pow--;
//  }
//
//  return TRUE;
//}

//================================================================================
void Matrix::SquareTranspose()
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Check
    if( Rows() != Cols() )
        LzLogException("", "Matrix is not square!" );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < i      ; j++ )
    {
        double lTmp = Elt( i, j );
        Elt( i, j ) = Elt( j, i );
        Elt( j, i ) = lTmp;
    }
}

//================================================================================
void Matrix::TransposeTo( Matrix & pRes ) const
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Check
    if( this == &pRes )
        LzLogException("", "Cannot transpose to self!" );

    pRes.SetDims( Cols(), Rows() );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pRes.Elt( j, i ) = Elt( i, j );
}

//================================================================================
bool Matrix::IsSymmetric( double pTol/*=-1*/ ) const
{
    // Check
    if( !Rows() || Rows() != Cols() )
        LzLogException("", "Invalid dimensions!" )

    // Check all elements
    for( unsigned int i = 0   ; i < Rows() - 1 ; i++ )
        for( unsigned int j = i + 1 ; j < Cols()   ; j++ )
        {
            if( !LzMath::ToolBox::IsZero( Elt( i, j ) - Elt( j, i ), pTol ) )
            {
//#ifdef USE_CLI
//                LzLogMsg("", "Matrix not symmetric: Elt(" + i + "," + j + ")=" + Elt( i, j ) + " != Elt(" + j + "," + i + ")=" + Elt( j, i ) + "." );
//#else
                LzLogMsg("", "Matrix not symmetric: Elt("
                          << i << "," << j << ")="
                          << Elt( i, j )
                          << " != Elt(" << j << "," << i
                          << ")=" << Elt( j, i ) << "." )
//#endif
                return false;
            }
        }

    return true;
}

//================================================================================
const Matrix & Matrix::operator=( const Matrix & pIn )
{
    SetDims( pIn.Rows(), pIn.Cols() );

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        Elt( i, j ) = pIn.Elt( i, j );

    return *this;
}

//================================================================================
double Matrix::DotProd( const Matrix & pA ) const
{
    // Check
    if( !Rows() || !Cols() || Rows() != pA.Rows() || Cols() != pA.Cols() )
        LzLogException("", "Invalid dimensions!")

    double lDotProd = 0;

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        lDotProd += Elt( i, j ) * pA.Elt( i, j );

    return lDotProd;
}

//================================================================================
double Matrix::MaxAbs() const
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!")

    double lMaxAbs = 0;

    for( unsigned int i = 0 ; i < Rows() ; i++ )
    for( unsigned int j = 0 ; j < Cols() ; j++ )
    {
        double lAbs = std::abs( Elt( i, j ) );

        if( lMaxAbs < lAbs )
            lMaxAbs = lAbs;
    }

    return lMaxAbs;
}

//================================================================================
double Matrix::Norm2() const
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!")

    double Norm2 = 0;

    for( unsigned int i = 0 ; i < mH ; i++ )
    for( unsigned int j = 0 ; j < mW ; j++ )
    {
        double lE = Elt( i, j );
        Norm2 += lE * lE;
    }

    return Norm2;
}

//================================================================================
double Matrix::M2x2_Det() const
{
    // Check
    if( Rows() != 2 || Cols() != 2 )
        LzLogException("", "Invalid dimensions!")

    // Read all 2x2 elements
    double a = Elt( 0, 0 );
    double b = Elt( 0, 1 );

    double c = Elt( 1, 0 );
    double d = Elt( 1, 1 );

    // Determinant
    return a * d - c * b;
}

//================================================================================
double Matrix::M3x3_Det() const
{
    // Check
    if( Rows() != 3 || Cols() != 3 )
        LzLogException("", "Invalid dimensions!")

    // Read all 3x3 elements
    double a = Elt( 0, 0 );
    double b = Elt( 0, 1 );
    double c = Elt( 0, 2 );

    double d = Elt( 1, 0 );
    double e = Elt( 1, 1 );
    double f = Elt( 1, 2 );

    double g = Elt( 2, 0 );
    double h = Elt( 2, 1 );
    double i = Elt( 2, 2 );

    // Determinant
    return a * ( e * i - h * f ) - b * ( d * i - g * f ) + c * ( d * h - g * e );
}

//================================================================================
void Matrix::M3x3_EigenPolynom( double & p3, double & p2, double & p1, double & p0 ) const
{
    // Check
    if( Rows() != 3 || Cols() != 3 )
        LzLogException("", "Invalid dimensions!")

    p3 = 1;

    p2 = -( Elt( 0, 0 ) + Elt( 1, 1 ) + Elt( 2, 2 ) );

    p1 = Elt( 0, 0 ) * Elt( 1, 1 ) + Elt( 0, 0 ) * Elt( 2, 2 ) + Elt( 1, 1 ) * Elt( 2, 2 ) - ( Elt( 0, 1 ) * Elt( 1, 0 ) + Elt( 0, 2 ) * Elt( 2, 0 ) + Elt( 1, 2 ) * Elt( 2, 1 ) );

    p0 = -M3x3_Det();
}

//================================================================================
double Matrix::M4x4_Det() const
{
    // Check
    if( Rows() != 4 || Cols() != 4 )
        LzLogException("", "Invalid dimensions!")

    // Read all 4x4 elements
    double A = Elt( 0, 0 );
    double B = Elt( 0, 1 );
    double C = Elt( 0, 2 );
    double D = Elt( 0, 3 );

    double E = Elt( 1, 0 );
    double F = Elt( 1, 1 );
    double G = Elt( 1, 2 );
    double H = Elt( 1, 3 );

    double I = Elt( 2, 0 );
    double J = Elt( 2, 1 );
    double K = Elt( 2, 2 );
    double L = Elt( 2, 3 );

    double M = Elt( 3, 0 );
    double N = Elt( 3, 1 );
    double O = Elt( 3, 2 );
    double P = Elt( 3, 3 );

    // Subcofactors
    double KP_OL = K * P - O * L;
    double JP_NL = J * P - N * L;
    double JO_NK = J * O - N * K;
    double IP_ML = I * P - M * L;
    double IO_MK = I * O - M * K;
    double IN_MJ = I * N - M * J;

    // Cofactors
    double CF_00 = +( F * ( KP_OL ) - G * ( JP_NL ) + H * ( JO_NK ) );
    double CF_01 = -( E * ( KP_OL ) - G * ( IP_ML ) + H * ( IO_MK ) );
    double CF_02 = +( E * ( JP_NL ) - F * ( IP_ML ) + H * ( IN_MJ ) );
    double CF_03 = -( E * ( JO_NK ) - F * ( IO_MK ) + G * ( IN_MJ ) );

    // Determinant
    return A * CF_00 + B * CF_01 + C * CF_02 + D * CF_03;
}

//================================================================================
void Matrix::M2x2_InverseTo( Matrix & pInv ) const
{
    // Check
    if( Rows() != 2 || Cols() != 2 )
        LzLogException("", "Invalid dimensions!")

    // Read all 2x2 elements
    double a = Elt( 0, 0 );
    double b = Elt( 0, 1 );

    double c = Elt( 1, 0 );
    double d = Elt( 1, 1 );

    // Determinant
    double lDet = a * d - c * b;
    if( LzMath::ToolBox::IsZero( lDet ) )
        LzLogException("", "Matrix is not inversible!")

    // Compute inverse as (transpose of cofactors matrix) / det
    pInv.SetDims( 2, 2 );
    pInv.Elt( 0, 0 ) = +d / lDet;
    pInv.Elt( 0, 1 ) = -b / lDet;
    pInv.Elt( 1, 0 ) = -c / lDet;
    pInv.Elt( 1, 1 ) = +a / lDet;
}

//================================================================================
bool Matrix::M2x2_InverseTo_NoExc( Matrix & pInv ) const
{
    // Check
    if( Rows() != 2 || Cols() != 2 )
	{
        LzLogErr("", "Invalid dimensions!")
		return false;
	}

    // Read all 2x2 elements
    double a = Elt( 0, 0 );
    double b = Elt( 0, 1 );

    double c = Elt( 1, 0 );
    double d = Elt( 1, 1 );

    // Determinant
    double lDet = a * d - c * b;
    if( LzMath::ToolBox::IsZero( lDet ) )
		return false;

    // Compute inverse as (transpose of cofactors matrix) / det
    pInv.SetDims( 2, 2 );
    pInv.Elt( 0, 0 ) = +d / lDet;
    pInv.Elt( 0, 1 ) = -b / lDet;
    pInv.Elt( 1, 0 ) = -c / lDet;
    pInv.Elt( 1, 1 ) = +a / lDet;

	return true;
}

//================================================================================
void Matrix::M3x3_InverseTo( Matrix & pInv ) const
{
    // Check
    if( Rows() != 3 || Cols() != 3 )
        LzLogException("", "Invalid dimensions!")

    // Read all 3x3 elements
    double a = Elt( 0, 0 );
    double b = Elt( 0, 1 );
    double c = Elt( 0, 2 );

    double d = Elt( 1, 0 );
    double e = Elt( 1, 1 );
    double f = Elt( 1, 2 );

    double g = Elt( 2, 0 );
    double h = Elt( 2, 1 );
    double i = Elt( 2, 2 );

    // Cofactors
    double CF[3][3];
    CF[0][0] = +( e * i - h * f );
    CF[0][1] = -( d * i - g * f );
    CF[0][2] = +( d * h - g * e );

    // Determinant
    double lDet = a * CF[0][0] + b * CF[0][1] + c * CF[0][2];
    if( LzMath::ToolBox::IsZero( lDet ) )
        LzLogException("", "Matrix is not inversible!")

    // Other cofactors
    CF[1][0] = -( b * i - h * c );
    CF[1][1] = +( a * i - g * c );
    CF[1][2] = -( a * h - g * b );

    CF[2][0] = +( b * f - e * c );
    CF[2][1] = -( a * f - d * c );
    CF[2][2] = +( a * e - d * b );

    // Compute inverse as (transpose of cofactors matrix) / det
    pInv.SetDims( 3, 3 );
    for( int i = 0 ; i < 3 ; i++ )
        for( int j = 0 ; j < 3 ; j++ )
            pInv.Elt( i, j ) = CF[j][i] / lDet;
}

//================================================================================
void Matrix::M4x4_InverseTo( Matrix & pInv ) const
{
    // Check
    if( Rows() != 4 || Cols() != 4 )
        LzLogException("", "Invalid dimensions!")

    // Read all 4x4 elements
    double A = Elt( 0, 0 );
    double B = Elt( 0, 1 );
    double C = Elt( 0, 2 );
    double D = Elt( 0, 3 );

    double E = Elt( 1, 0 );
    double F = Elt( 1, 1 );
    double G = Elt( 1, 2 );
    double H = Elt( 1, 3 );

    double I = Elt( 2, 0 );
    double J = Elt( 2, 1 );
    double K = Elt( 2, 2 );
    double L = Elt( 2, 3 );

    double M = Elt( 3, 0 );
    double N = Elt( 3, 1 );
    double O = Elt( 3, 2 );
    double P = Elt( 3, 3 );

    // Subcofactors
    double KP_OL = K * P - O * L;
    double JP_NL = J * P - N * L;
    double JO_NK = J * O - N * K;
    double IP_ML = I * P - M * L;
    double IO_MK = I * O - M * K;
    double IN_MJ = I * N - M * J;

    // Cofactors
    double CF[4][4];
    CF[0][0] = +( F * ( KP_OL ) - G * ( JP_NL ) + H * ( JO_NK ) );
    CF[0][1] = -( E * ( KP_OL ) - G * ( IP_ML ) + H * ( IO_MK ) );
    CF[0][2] = +( E * ( JP_NL ) - F * ( IP_ML ) + H * ( IN_MJ ) );
    CF[0][3] = -( E * ( JO_NK ) - F * ( IO_MK ) + G * ( IN_MJ ) );

    // Check determinant
    double lDet = A * CF[0][0] + B * CF[0][1] + C * CF[0][2] + D * CF[0][3];
    if( LzMath::ToolBox::IsZero( lDet ) )
        LzLogException("", "Matrix is not inversible!")

    // Other cofactors
    CF[1][0] = -( B * ( KP_OL ) - C * ( JP_NL ) + D * ( JO_NK ) );
    CF[1][1] = +( A * ( KP_OL ) - C * ( IP_ML ) + D * ( IO_MK ) );
    CF[1][2] = -( A * ( JP_NL ) - B * ( IP_ML ) + D * ( IN_MJ ) );
    CF[1][3] = +( A * ( JO_NK ) - B * ( IO_MK ) + C * ( IN_MJ ) );

    double GP_OH = G * P - O * H;
    double FP_NH = F * P - N * H;
    double FO_NG = F * O - N * G;
    double EP_MH = E * P - M * H;
    double EO_MG = E * O - M * G;
    double EN_MF = E * N - M * F;

    CF[2][0] = +( B * ( GP_OH ) - C * ( FP_NH ) + D * ( FO_NG ) );
    CF[2][1] = -( A * ( GP_OH ) - C * ( EP_MH ) + D * ( EO_MG ) );
    CF[2][2] = +( A * ( FP_NH ) - B * ( EP_MH ) + D * ( EN_MF ) );
    CF[2][3] = -( A * ( FO_NG ) - B * ( EO_MG ) + C * ( EN_MF ) );

    double GL_KH = G * L - K * H;
    double FL_JH = F * L - J * H;
    double FK_JG = F * K - J * G;
    double EL_IH = E * L - I * H;
    double EK_IG = E * K - I * G;
    double EJ_IF = E * J - I * F;

    CF[3][0] = -( B * ( GL_KH ) - C * ( FL_JH ) + D * ( FK_JG ) );
    CF[3][1] = +( A * ( GL_KH ) - C * ( EL_IH ) + D * ( EK_IG ) );
    CF[3][2] = -( A * ( FL_JH ) - B * ( EL_IH ) + D * ( EJ_IF ) );
    CF[3][3] = +( A * ( FK_JG ) - B * ( EK_IG ) + C * ( EJ_IF ) );

    // Compute inverse as (transpose of cofactors matrix) / det
    pInv.SetDims( 4, 4 );
    for( int i = 0 ; i < 4 ; i++ )
        for( int j = 0 ; j < 4 ; j++ )
            pInv.Elt( i, j ) = CF[j][i] / lDet;
}

//================================================================================
void Matrix::SubMatrixTo( unsigned int pI0, unsigned int pJ0, Matrix & pTo ) const
{
    // Check
    if( !pTo.Rows() || !pTo.Cols() || Rows() < pI0 + pTo.Rows() || Cols() < pJ0 + pTo.Cols() )
        LzLogException("", "Invalid dimensions!")

    for( unsigned int i = 0 ; i < pTo.Rows() ; i++ )
    for( unsigned int j = 0 ; j < pTo.Cols() ; j++ )
        pTo.Elt( i, j ) = Elt( pI0 + i, pJ0 + j );
}

//================================================================================
void Matrix::SubMatrixFrom( unsigned int pI0, unsigned int pJ0, const Matrix & pFrom )
{
    // Check
    if( !pFrom.Rows() || !pFrom.Cols() || Rows() < pI0 + pFrom.Rows() || Cols() < pJ0 + pFrom.Cols() )
        LzLogException("", "Invalid dimensions!")

    for( unsigned int i = 0 ; i < pFrom.Rows() ; i++ )
    for( unsigned int j = 0 ; j < pFrom.Cols() ; j++ )
        Elt( pI0 + i, pJ0 + j ) = pFrom.Elt( i, j );
}


#pragma region UNUSED
//////////////////////////////////////////////////////////////////////////////////
////
//// Permutations
////
//================================================================================
////BOOL CMatrix::LoadPermutation( const CContainer1<int> & p_Mapping )
////{
////    // Check dimensions
////    if( !SetDims(p_Mapping.Width(),p_Mapping.Width()) )
////    {
////        LOGERR("[CMatrix::LoadPermutation] : Dimensions error !");
////        return FALSE;
////    }
////
////    LoadValue(0);
////
////    // Fill
////    for( int i=0 ; i<m_H ; i++ )
////    {
////        int j = p_Mapping[i];
////        if( j<0 || j>=m_W )
////        {
////            LOGERR("[CMatrix::LoadPermutation] : Mapping error (i=%d,j=%d) !",i,j);
////            return FALSE;
////        }
////
////        Elt(i,j) = 1;
////    }
////
////    // Check
////    for( int j=0 ; j<m_W ; j++ )
////    {
////        int count = 0;
////        for( i=0 ; i<m_H ; i++ )
////        {
////            if( Elt(i,j) > 0.5 )
////                count++;
////        }
////
////        if( count != 1 )
////        {
////            LOGERR("[CMatrix::LoadPermutation] : Mapping error : col (j=%d) has %d 1s !",j,count);
////            return FALSE;
////        }
////    }
////
////    return TRUE;
////}

//////////////////////////////////////////////////////////////////////////////////
////
//// Rotations
////
//================================================================================
//BOOL CMatrix::LoadRotation( const CVector & pRoot, const CVector & pDir, double pDeg )
//{
//  // Check dims
//  if( pRoot.Rows()!=3 || pDir.Rows()!=3 )
//  {
//      LOGERR("[CMatrix::LoadRotation] : Invalid dims !");
//      return FALSE;
//  }

//  // Normalize dir
//  CVector Dir = pDir;
//  if( !Dir.Normalize() )
//  {
//      LOGERR("[CMatrix::LoadRotation] : Invalid rotation axis !");
//      return FALSE;
//  }

//  // Deg to rad
//  double Rad = pDeg * TxServices::PI / 180.0;

//  // Scal mat
//  CMatrix S;
//  Dir.ExtrProd( Dir, S );

//  // Cross mat
//  CMatrix W( 3,3, 0.0, -Dir.Elt(2), Dir.Elt(1),
//                  Dir.Elt(2), 0.0 ,-Dir.Elt(0),
//                  -Dir.Elt(1), Dir.Elt(0), 0.0 );
//  // Sin mat
//  CMatrix sinW;
//  W.Mult( sin(Rad), sinW );

//  // Cos mat
//  CMatrix cosWW;
//  W.Mult( W, cosWW );
//  cosWW.Mult( -cos(Rad), cosWW );

//  // Total cos mat + sin mat
//  CMatrix Sum;
//  S.Add( cosWW, Sum );
//  Sum.Add( sinW, Sum );

//  // This = pure rotation
//  SetDims( 4, 4 );
//  for( int i=0 ; i<3 ; i++ )
//  for( int j=0 ; j<3 ; j++ )
//  {
//      Elt(i,j) = Sum.Elt(i,j);
//  }
//  for( i=0 ; i<3 ; i++ )
//  {
//      Elt(i,3) = 0;
//      Elt(3,i) = 0;
//  }
//  Elt(3,3) = 1;

//  // Translation to / from origin
//  CMatrix T, U;
//  T.SetDims( 4, 4 );
//  U.SetDims( 4, 4 );
//  T.LoadIdentity();
//  U.LoadIdentity();
//  T.Elt(0,3) = pRoot.Elt(0);
//  T.Elt(1,3) = pRoot.Elt(1);
//  T.Elt(2,3) = pRoot.Elt(2);
//  U.Elt(0,3) = -pRoot.Elt(0);
//  U.Elt(1,3) = -pRoot.Elt(1);
//  U.Elt(2,3) = -pRoot.Elt(2);

//  // This = all together
//  CMatrix RU;
//  Mult( U, RU );
//  T.Mult( RU, *this );

//  return TRUE;
//}
//================================================================================
//BOOL CMatrix::SubsInAllRowsFor( int p_Var, CVector & p_F )
//{
//  // Check if matrix and vector are ok ?
//  if( !Rows() || Rows()!=Cols() || Rows()!=p_F.Rows() )
//  {
//      LOGERR("[CMatrix::SubsInAllRowsFor] : Wrong dims !");
//      return FALSE;
//  }

//  // Check if variable ok
//  if( p_Var<0 || p_Var>=Rows() )
//  {
//      LOGERR("[CMatrix::SubsInAllRowsFor] : Wrong variable !");
//      return FALSE;
//  }

//  // Divisor
//  double lPivot = Elt( p_Var, p_Var );
//  if( TxMathBox::IsZero(fabs(lPivot)) )
//  {
//      LOGERR("[CMatrix::SubsInAllRowsFor] : Could not substitute for variable %d ! Null coef.",p_Var);
//      return FALSE;
//  }

//  // Substitute in all rows
//  for( int i=0 ; i<Rows() ; i++ )
//  {
//      // Skip p_Var row
//      if( i == p_Var )
//          continue;

//      // Change constant in vector
//      p_F.Elt(i) = p_F.Elt(i) - p_F.Elt(p_Var)*Elt(i,p_Var)/lPivot;

//      // Change coefs in matrix
//      for( int j=0 ; j<Cols() ; j++ )
//      {
//          // Skip the col corresponding to p_Var (will be reset later)
//          if( j == p_Var )
//              continue;

//          Elt(i,j) = Elt(i,j) - Elt(i,p_Var)*Elt(p_Var,j)/lPivot;
//      }

//      // Reset coef for p_Var
//      Elt(i,p_Var) = 0;
//  }

//  return TRUE;
//}

//================================================================================
//BOOL CMatrix::KillRowCols( const CLst<int> & pRows, const CLst<int> & pCols, CMatrix & pToMat ) const
//{
//  // Check empty matrix
//  if( !Rows() || !Cols() )
//  {
//      LOGERR("[CMatrix::KillRowCols] : Empty matrix !");
//      return FALSE;
//  }

//  // Check re-entrance
//#pragma message("----- AJOUTER CHECKS DE REENTRANCE PARTOUT ? (#ifdef _DEBUG)")
//  if( &pToMat == this )
//  {
//      LOGERR("[CMatrix::KillRowCols] : Re-entering same matrix !");
//      return FALSE;
//  }

//  // Iterator
//  void * iPos;

//  // Check row idx ok
//  for( iPos=pRows.HeadPos() ; iPos ; pRows.Next(iPos) )
//  {
//      int lIdx = pRows.GetAt( iPos );

//      if( lIdx<0 || lIdx>=Rows() )
//      {
//          LOGERR("[CMatrix::KillRowCols] : Invalid row idx (%d) !",lIdx);
//          return FALSE;
//      }
//  }

//  // Check col idx ok
//  for( iPos=pCols.HeadPos() ; iPos ; pCols.Next(iPos) )
//  {
//      int lIdx = pCols.GetAt( iPos );

//      if( lIdx<0 || lIdx>=Cols() )
//      {
//          LOGERR("[CMatrix::KillRowCols] : Invalid col idx (%d) !",lIdx);
//          return FALSE;
//      }
//  }

//  // Reduced size
//  int lNewRows = Rows();
//  int lNewCols = Cols();

//  // Build bool vectors
//  static CVector lRowIn, lColIn;
//  lRowIn.SetDim( Rows() );
//  lColIn.SetDim( Cols() );
//  lRowIn.LoadValue( 1 );
//  lColIn.LoadValue( 1 );

//  // Fill row bool vector
//  for( iPos=pRows.HeadPos() ; iPos ; pRows.Next(iPos) )
//  {
//      int lIdx = pRows.GetAt( iPos );

//      // In case of repeated idx
//      if( lRowIn.Elt(lIdx) )
//          lNewRows--;

//      lRowIn.Elt(lIdx) = 0;
//  }

//  // Fill col bool vector
//  for( iPos=pCols.HeadPos() ; iPos ; pCols.Next(iPos) )
//  {
//      int lIdx = pCols.GetAt( iPos );

//      // In case of repeated idx
//      if( lColIn.Elt(lIdx) )
//          lNewCols--;

//      lColIn.Elt(lIdx) = 0;
//  }

//  // Size matrix
//  if( !pToMat.SetDims(lNewRows,lNewCols) )
//  {
//      LOGERR("[CMatrix::KillRowCols] : Unable to set mat dims !");
//      return FALSE;
//  }

//  // For all rows
//  int lNewI = 0;
//  for( int i=0 ; i<Rows() ; i++ )
//  {
//      // This row is in ?
//      if( !lRowIn.Elt(i) )
//          continue;

//      // For all cols
//      int lNewJ = 0;
//      for( int j=0 ; j<Cols() ; j++ )
//      {
//          // This col is in ?
//          if( !lColIn.Elt(j) )
//              continue;

//          // Copy element
//          pToMat.Elt(lNewI,lNewJ) = Elt(i,j);

//          // Next new col
//          lNewJ++;
//      }

//      // Next new row
//      lNewI++;
//  }

//  return TRUE;
//}
#pragma endregion
#pragma endregion


#pragma region "Load & save"
//================================================================================
void Matrix::SaveBin( const string & /*pFileName*/ ) const
{
#ifdef USE_CLI
	// Writer
	System::IO::BinaryWriter ^lWriter = nullptr;

	_TxLogTry
	lWriter = gcnew System::IO::BinaryWriter( System::IO::File::Open( gcnew String(pFileName.c_str()), System::IO::FileMode::Create ) );

	// Header
	const unsigned char lHeader[5] = { 'T', 'x', 'M', 'a', 't' };
	for( int h=0 ; h<5 ; h++ )
		lWriter->Write( (unsigned char)lHeader[h] );

	// Size
	lWriter->Write( (System::UInt32)Rows() );
	lWriter->Write( (System::UInt32)Cols() );

	// Write mat elements
	for( unsigned int i=0 ; i<mH ; i++ )
	for( unsigned int j=0 ; j<mW ; j++ )
		lWriter->Write( Elt(i,j) );

	_TxLogCatchAndThrow("Could not save matrix! File= "<<pFileName<<".")
	finally
	{
		if( lWriter )
			lWriter->Close();
	}
#else
    LzLogException("","*** TODO: Matrix::SaveBin( const std::string & pFname )")
#endif
}

//================================================================================
void Matrix::LoadBin( const std::string & /*pFileName*/ )
{
#ifdef USE_CLI
	// Secure perimeter
	TxServices::Detonator lDet( std::bind(&Matrix::Free, this) );

	// Reader
	System::IO::BinaryReader ^lReader = nullptr;
		
	// Result
	bool lOk = false;

	_TxLogTry

	// Free previous
	Free(); 

	// Open file
	System::IO::FileStream ^lStream = System::IO::File::Open( gcnew String(pFileName.c_str()), System::IO::FileMode::Open );
	lReader = gcnew System::IO::BinaryReader( lStream );

	// Read header
	const unsigned char lHeader[5] = { 'T', 'x', 'M', 'a', 't' };
	for( int h=0 ; h<5 ; h++ )
	{
		if( (char)lReader->ReadByte() != lHeader[h] )
            LzLogException("","Invalid matrix header!");
	}

	// Read size
	unsigned int lRows = lReader->ReadUInt32();
	unsigned int lCols = lReader->ReadUInt32();

	// Set size
	SetDims( lRows, lCols );

	// Read values
	for( unsigned int i=0 ; i<mH ; i++ )
	for( unsigned int j=0 ; j<mW ; j++ )
		Elt(i,j) = lReader->ReadDouble();

	// All is well
	lDet.Defuse();

	_TxLogCatchAndThrow("Could not load matrix! File= "<<pFileName<<".")
	finally
	{
		if( lReader )
			lReader->Close();
	}
#else
    LzLogException("","*** TODO: Matrix::LoadBin( const std::string & pFileName )")
#endif
}
#pragma endregion


#pragma region "Co-factors"
//================================================================================
List<int> Matrix::mLocalRows;
List<int> Matrix::mLocalCols;

//================================================================================
double Matrix::CF_Det( bool pEmptyMat/*=false*/ ) const
{
    if( !Rows() || Rows() != Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Init local cols and rows
    InitLocalRowsCols();

    return CF_RecDet( pEmptyMat );
}

//================================================================================
int Matrix::MaxZeroLocalRow() const
{
    int lMaxZeroCount = -1;
    int lLocalRow = 0;
    int lBestRow;

    BrowseList( iRowPos, mLocalRows )
    {
        int lGlobalRow = mLocalRows.GetAt( iRowPos );
        int lZeroCount = 0;

        BrowseList( iColPos, mLocalCols )
        {
            int lGlobalCol = mLocalCols.GetAt( iColPos );
            if( LzMath::ToolBox::IsZero( Elt( lGlobalRow, lGlobalCol ) ) )
                lZeroCount++;
        }

        if( lZeroCount > lMaxZeroCount )
        {
            lMaxZeroCount = lZeroCount;
            lBestRow = lLocalRow;
        }

        lLocalRow++;
    }

    return lBestRow;
}

//================================================================================
double Matrix::CF_RecDet( bool pEmptyMat/*=false*/ ) const
{
    // 0x0 ?
    if( LocalCols() == 0 )
        return 1;
    else
    // 1x1 ?
    if( LocalCols() == 1 )
        return Elt( mLocalRows.GetHead(), mLocalCols.GetHead() );
    else
    // 2x2 ?
    if( LocalCols() == 2 )
    {
        double l00 = Elt( mLocalRows.GetHead(), mLocalCols.GetHead() );
        double l11 = Elt( mLocalRows.GetTail(), mLocalCols.GetTail() );
        double l01 = Elt( mLocalRows.GetHead(), mLocalCols.GetTail() );
        double l10 = Elt( mLocalRows.GetTail(), mLocalCols.GetHead() );

        return l00 * l11 - l01 * l10;
    }

    // Decrease dim and recurse

    // Remove first local row
    int lLocalI;
    if( pEmptyMat )
        lLocalI = MaxZeroLocalRow();
    else
        lLocalI = 0;

    int lGlobalRow = RemoveLocalRow( lLocalI );

    // Develop along first row
    double lDet = 0;
    for( int iLocalJ = 0 ; iLocalJ < LocalCols() ; iLocalJ++ )
    {
        // Remove local col
        int lGlobalCol = RemoveLocalCol( iLocalJ );

        // Compute sub det only if element non null
        double lElt = Elt( lGlobalRow, lGlobalCol );
        if( !LzMath::ToolBox::IsZero( fabs( lElt ) ) )
            lDet += ( ( lLocalI + iLocalJ ) % 2 ? -1 : +1 ) * lElt * CF_RecDet( pEmptyMat );

        // Add global col for further recursion
        AddGlobalCol( lGlobalCol );
    }

    // Add global row for further recursion
    AddGlobalRow( lGlobalRow );

    return lDet;
}

//================================================================================
void Matrix::InitLocalRowsCols() const
{
    mLocalRows.DelAll();
    mLocalCols.DelAll();
    for( unsigned int i = 0 ; i < mH ; i++ )
    {
        mLocalRows.AddTail( i );
        mLocalCols.AddTail( i );
    }
}
int Matrix::RemoveLocalRow( int pLocalI ) const
{
    void * lPos = mLocalRows.FindIdx( pLocalI );
    int lGlobalRow = mLocalRows.GetAt( lPos );
    mLocalRows.DelAt( lPos );
    return lGlobalRow;
}
int Matrix::RemoveLocalCol( int pLocalJ ) const
{
    void * lPos = mLocalCols.FindIdx( pLocalJ );
    int lGlobalCol = mLocalCols.GetAt( lPos );
    mLocalCols.DelAt( lPos );
    return lGlobalCol;
}
void Matrix::AddGlobalRow( int pGlobalI ) const
{
    void * iPos = mLocalRows.HeadPos();
    while( iPos && mLocalRows.GetAt( iPos ) < pGlobalI )
        mLocalRows.Next( iPos );
    if( iPos )
        mLocalRows.AddBefore( iPos, pGlobalI );
    else
        mLocalRows.AddTail( pGlobalI );
}
void Matrix::AddGlobalCol( int pGlobalJ ) const
{
    void * iPos = mLocalCols.HeadPos();
    while( iPos && mLocalCols.GetAt( iPos ) < pGlobalJ )
        mLocalCols.Next( iPos );
    if( iPos )
        mLocalCols.AddBefore( iPos, pGlobalJ );
    else
        mLocalCols.AddTail( pGlobalJ );
}

//================================================================================
void Matrix::CF_InverseTo( Matrix & pInv, bool pEmptyMat/*=false*/ )
{
    // Check
    if( !Rows() || Rows() != Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Inversible ?
    double lDet = CF_Det( pEmptyMat );
    if( LzMath::ToolBox::IsZero( fabs( lDet ) ) )
        LzLogException("", "Matrix not inversible!" );

    // Dims
    pInv.SetDims( Rows(), Cols() );

    // Fill with transpose(co-factors)
    for( unsigned int i = 0 ; i < Rows() ; i++ )
        for( unsigned int j = 0 ; j < Cols() ; j++ )
        {
            // Remove
            RemoveLocalRow( i );
            RemoveLocalCol( j );

            // Compute co-factor and transpose
            pInv.Elt( j, i ) = ( ( i + j ) % 2 ? -1 : +1 ) * CF_RecDet( pEmptyMat ) / lDet;

            // Add
            AddGlobalRow( i );
            AddGlobalCol( j );
        }
}
#pragma endregion


#pragma region "LU decomposition"
//================================================================================
void Matrix::LU_ify( Matrix & pL, Matrix & pU, double pTol )
{
    if( !Rows() || Rows() != Cols() )
        LzLogException("", "Invalid dimensions!" );

    // Dimension and reset recipients
    pL.SetDims( Cols(), Cols() );
    pU.SetDims( Cols(), Cols() );
    pL.LoadValue( 0 );
    pU.LoadValue( 0 );

    // Init first line U
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pU.Elt( 0, j ) = Elt( 0, j );

    // Init diagonal L
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        pL.Elt( j, j ) = 1;

    // Fill columns(L) and lines(U)
    unsigned int iColRow = 0;
    while( iColRow < Cols() - 1 )
    {
        // Have solution ?
        double lK = pU.Elt( iColRow, iColRow );

        if( LzMath::ToolBox::IsZero( lK, pTol ) )
        {
//#ifdef USE_CLI
//            LzLogException("", "Unable to LU-ify! Found an excessively small diagonal element= " + lK + "." );
//#else
            LzLogException("", "Unable to LU-ify!"
                            << " Found an excessively small diagonal element= "
                            << lK << "." )
//#endif
        }

        // Fill col(L)
        for( unsigned int i = iColRow + 1 ; i < Cols() ; i++ )
        {
            // Calc L( row=i, col=iColRow )
            double lSum = 0;
            for( unsigned int k = 0 ; k < iColRow ; k++ )
                lSum += pL.Elt( i, k ) * pU.Elt( k, iColRow );

            pL.Elt( i, iColRow ) = ( Elt( i, iColRow ) - lSum ) / lK;
        }

        // Next row for U
        iColRow++;

        // Fill row(U)
        for( unsigned int j = iColRow ; j < Cols() ; j++ )
        {
            // Calc U( row=i_ColRow, col=j )
            double lSum = 0;
            for( unsigned int k = 0 ; k < iColRow ; k++ )
                lSum += pL.Elt( iColRow, k ) * pU.Elt( k, j );

            pU.Elt( iColRow, j ) = Elt( iColRow, j ) - lSum;
        }
    }
}

//================================================================================
double Matrix::LU_Det()
{
    Matrix L, U;
    LU_ify( L, U, -1 );

    return U.TriDet();
}

//================================================================================
double Matrix::TriDet()
{
    // Check
    if( !Rows() || Rows() != Cols() )
        LzLogException("", "Invalid dimensions!")

    double lDet = 1;
    for( unsigned int i = 0 ; i < Cols() ; i++ )
        lDet *= Elt( i, i );

    return lDet;
}

//================================================================================
void Matrix::L_InverseTo( Matrix & pInvL )
{
    // Invertible ?
    const double lTD = TriDet();
    if( LzMath::ToolBox::IsZero( lTD ) )
    {
//#ifdef USE_CLI
//        LzLogException("", "Matrix not invertible! (tri. det.= " + lTD + ")" );
//#else
        LzLogException("", "Matrix not invertible! (tri. det.= " << lTD << ")")
//#endif
    }

    // Fill-in inverse col after col
    pInvL.SetDims( Rows(), Cols() );
    for( unsigned int j = 0 ; j < Cols() ; j++ )
        for( unsigned int i = 0 ; i < Rows() ; i++ )
        {
            if( i < j )
                pInvL.Elt( i, j ) = 0;
            else if( i == j )
                pInvL.Elt( i, i ) = 1 / Elt( i, i );
            else
            {
                double lSum = 0;
                for( unsigned int k = j ; k < i ; k++ )
                    lSum += Elt( i, k ) * pInvL.Elt( k, j );

                pInvL.Elt( i, j ) = -lSum / Elt( i, i );
            }
        }
}

//================================================================================
void Matrix::U_InverseTo( Matrix & pInvU )
{
    // Invertible ?
    if( LzMath::ToolBox::IsZero( TriDet() ) )
        LzLogException("", "Matrix not invertible!")

    // Fill-in inverse backwards col after col
    pInvU.SetDims( Rows(), Cols() );
    for( int j = Cols() - 1 ; j >= 0 ; j-- )
        for( int i = Rows() - 1 ; i >= 0 ; i-- )
        {
            if( i > j )
                pInvU.Elt( i, j ) = 0;
            else if( i == j )
                pInvU.Elt( i, i ) = 1 / Elt( i, i );
            else
            {
                double lSum = 0;
                for( int k = Rows() - 1 ; k > i ; k-- )
                    lSum += Elt( i, k ) * pInvU.Elt( k, j );

                pInvU.Elt( i, j ) = -lSum / Elt( i, i );
            }
        }
}

//================================================================================
void Matrix::LU_InverseTo( Matrix & pInv )
{
    // Decomposition
    Matrix L, U;
    LU_ify( L, U, -1 );

    // Inverse each triangular mat
    Matrix InvL, InvU;
    L.L_InverseTo( InvL );
    U.U_InverseTo( InvU );

    // Put all together
    InvU.Mult( InvL, pInv );
}

//================================================================================
//BOOL CMatrix::LU_Solve( CVector & p_X, const CVector & p_Y )
//{
//  // This is a square non-empty matrix ?
//  if( !Rows() || Rows()!=Cols() )
//  {
//      LOGERR("[CMatrix::LU_Solve] : Empty or not square matrix !");
//      return FALSE;
//  }

//  // Y vector ok ?
//  if( Cols() != p_Y.Rows() )
//  {
//      LOGERR("[CMatrix::LU_Solve] : Dimensions mismatch !");
//      return FALSE;
//  }

//  CMatrix L,U;
//  if( !LU_ify(L,U) )
//  {
//      LOGERR("[CMatrix::LU_Solve] : Unable to LU-ify !");
//      return FALSE;
//  }

//  // Z <== Z / LZ=Y
//  CVector Z;
//  Z.SetDim( p_Y.Rows() );
//  for( int i=0 ; i<p_Y.Rows() ; i++ )
//  {
//      Z.Elt(i) = p_Y.Elt(i);
//
//      for( int j=i-1 ; j>=0 ; j-- )
//          Z.Elt(i) -= L.Elt(i,j)*Z.Elt(j);
//  }

//  // X <== X / UX=Z
//  p_X.SetDim( p_Y.Rows() );
//  for( i=p_Y.Rows()-1 ; i>=0 ; i-- )
//  {
//      p_X.Elt(i) = Z.Elt(i) / U.Elt(i,i);

//      for( int j=i+1 ; j<p_Y.Rows() ; j++ )
//          p_X.Elt(i) -= U.Elt(i,j)*p_X.Elt(j)/U.Elt(i,i);
//  }

//  return TRUE;
//}
#pragma endregion


#pragma region "Cholesky LLt: sym def pos matrix"
//================================================================================
void Matrix::LLt_ify( Matrix & p_L )
{
    // Check
    if( !Rows() || Rows() != Cols() )
        LzLogException("", "Invalid dimensions!")

    //-------------------------------------------------------------------
    // Only the lower part of (*this) matrix is accessed for computation
    //-------------------------------------------------------------------

    // Ok ?
    p_L.SetDims( Cols(), Cols() );

    // Reset
    p_L.LoadValue( 0 );

    //unsigned int i,j,k;

    // Fill col by col
    for( unsigned int j = 0 ; j < Cols() ; j++ )
    {
        // Ljj
        double Ljj_sum = Elt( j, j );
        for( unsigned int i = 0 ; i < j ; i++ )
            Ljj_sum -= p_L.Elt( j, i ) * p_L.Elt( j, i );

        if( Ljj_sum <= 0 )
        {
//#ifdef USE_CLI
//            LzLogException("", "Unable to LLt-ify! sum= " + Ljj_sum + "." );
//#else
            LzLogException("", "Unable to LLt-ify! sum= " << Ljj_sum << ".")
//#endif
        }


        p_L.Elt( j, j ) = sqrt( Ljj_sum );

        // Fill rows below Ljj
        for( unsigned int i = j + 1 ; i < Cols() ; i++ )
        {
            double Lij_sum = Elt( i, j );
            for( unsigned int k = 0 ; k < j ; k++ )
                Lij_sum -= p_L.Elt( j, k ) * p_L.Elt( i, k );

            p_L.Elt( i, j ) = Lij_sum / p_L.Elt( j, j );
        }
    }
}

//================================================================================
void Matrix::LLt_Solve( const Matrix & pL, Matrix & pX, const Matrix & pY )
{
    // @TODO CHECKER !!!! et modifier code pour autoriser pY.Cols()>1

    // Check
    if( !pL.Rows() || pL.Rows() != pL.Cols() || pY.Rows() != pL.Rows() || !pY.Cols() )
        LzLogException("", "Invalid dimensions!")

    //// L matrix ok ?
    //if( !p_L.mpMat || p_L.Rows()!=p_L.Cols() )
    //{
    //  LzLogErr("","Empty or not square L matrix!");
    //  return false;
    //}

    //// Y vector ok ?
    //if( p_Y.Rows()!=p_L.Cols() || p_Y.Cols()!=1 )
    //{
    //  LzLogErr("","Dimensions mismatch!");
    //  return false;
    //}

    // Check L
    for( unsigned int i = 0 ; i < pL.Rows() ; i++ )
    {
        if( LzMath::ToolBox::IsZero( pL.Elt( i, i ) ) )
            LzLogException("", "Unable to LLt solve!")
    }

    // Z <== Z / LZ=Y
    Matrix Z;
    Z.SetDims( pY.Rows(), 1 );
    for( unsigned int i = 0 ; i < Z.Rows() ; i++ )
    {
        Z.Elt( i, 0 ) = pY.Elt( i, 0 );

        for( unsigned int j = 0 ; j < i ; j++ )
            Z.Elt( i, 0 ) -= pL.Elt( i, j ) * Z.Elt( j, 0 );

        Z.Elt( i, 0 ) /= pL.Elt( i, i );
    }

    // X <== X / LtX=Z
    pX.SetDims( pY.Rows(), 1 );
    for( int i = pY.Rows() - 1 ; i >= 0 ; i-- )
    {
        pX.Elt( i, 0 ) = Z.Elt( i, 0 );

        for( unsigned int j = i + 1 ; j < pY.Rows() ; j++ )
            pX.Elt( i, 0 ) -= pL.Elt( j, i ) * pX.Elt( j, 0 );

        pX.Elt( i, 0 ) /= pL.Elt( i, i );
    }
}
#pragma endregion


#pragma region "Iterative solver: gradient and conjugate gradient"
//================================================================================
void Matrix::SolveSteepestDescent( Matrix & pU, const Matrix & pF, unsigned int /*pMaxIterations*/, double /*pErrMaxAbs*/ ) const
{
    // Assuming: *this = K
    Matrix lKU;

    //***************** CODE FROM SCHEWCHUK + commentaires: lien sur la publi dans la bib sur pc !

//************** checker code Grad_Solve vs. schewchuk pour voir s'il faut tout refaire
//************** conj grad: a priori oui...

    // Compute energy
    double lE;
    {
        Mult( pU, lKU );
        lE = 0.5 * pU.DotProd( lKU ) - pU.DotProd( pF );
    }

    // Compute gradient
    Matrix lGE;
    {
        Mult( pU, lKU );
        lKU.Sub( pF, lGE );
    }

//f(0) = E(0)
//f(1) = E(1)
//f'(0) = -lGE.DotProd(lGE);

    double df_dt;
    {
        df_dt = -lGE.DotProd( lGE );
    }


    //************ la meme chose en sparse matrix


    //************ et gradient conjugue : preconditionn� !!!!! (ca marchait pendant la these non? cf manuscrit)

}

//================================================================================
void Matrix::SolveConjGrad( Matrix & /*pX*/, const Matrix & /*pY*/, unsigned int /*pMaxIterations*/, double /*pTol*/ ) const
//BOOL CMatrix::Grad_Solve( CVector & p_X, const CVector & p_Y, double p_EMax ) const
{
    /*
    // Check dims
    if( Rows()!=Cols() || Rows()!=p_X.Rows() || Rows()!=p_Y.Rows() )
    LzLogException("","Dim mismatch!");

    // Value for current solution
    Matrix V;
    V.SetDims( Rows(), 1 );

    // Value at gradient
    Matrix VG;
    VG.SetDims( Rows(), 1 );

    // Gradient of error
    Matrix G;
    G.SetDims( Rows(), 1 );

    // Gradient descent
    unsigned int iStep = 0;
    while( iStep < pMaxIterations )
    {
        // Evaluate current solution
        Mult( p_X, V );


    //{
    //   CVector E;
    //   V.Sub( p_Y, E );
    //   LOGMSG("[%d] Error norm = %f",iStep,sqrt(E.Norm2()));
    //}

        // Compute gradient of error
        for( int j=0 ; j<Cols() ; j++ )
        {
            G.Elt(j) = 0;
            for( int i=0 ; i<Rows() ; i++ )
                G.Elt(j) += 2*Elt(i,j)*( V.Elt(i) - p_Y.Elt(i) );
        }

        // Compute gradient norm
        double lNormG2 = G.Norm2();

        // Reached minimum
        if( lNormG2 < pTol )
        {
            LzLogMsg("","Found solution in "+iStep+" step(s).");
            return;
        }

        // Compute min along gradient
        Mult( G, VG );
        double A = 0;
        double B = 0;
        for( int i=0 ; i<Rows() ; i++ )
        {
            A += VG.Elt(i) * VG.Elt(i);
            B += VG.Elt(i) * (V.Elt(i)-p_Y.Elt(i));
        }

        // Update current solution
        G.Mult( -B/A, G );
        p_X.Add( G, p_X );


    //p_X.Log("temp X");

        // Next step
        iStep++;
    }

    LzLogException("","Couldn't find solution after "+iStep+" step(s).");
    */
}

//================================================================================
//BOOL CMatrix::ConjGrad_Solve( CVector & p_X, const CVector & p_Y, double p_EMax ) const
//{
////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Implementation BOF !!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Implementation BOF !!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  // Check dims
//  if( Rows()!=Cols() || Rows()!=p_X.Rows() || Rows()!=p_Y.Rows() )
//  {
//      LOGERR("[CMatrix::Grad_Solve] : Dim mismatch !");
//      return FALSE;
//  }

//CGraphTools::StartPerfMeter();

//  // Residual
//  CVector Rk;
//  Rk.SetDim( Rows() );

//  // Conjugate direction
//  CVector Pk;
//  Pk.SetDim( Rows() );

//  // Init
//  CVector Tmp;
//  Mult( p_X, Tmp );
//  p_Y.Sub( Tmp, Rk );
//  double RkNorm2 = Rk.Norm2();

//  Pk = Rk;

//  // Loop
//  int iStep = 0;
//  while( iStep < 10000 )
//  {
//      // Alpha k
//      Mult( Pk, Tmp );
//      double lPk_APk;
//      Pk.ScalProd( Tmp, lPk_APk );
//      double Alpha = RkNorm2 / lPk_APk;

//      p_X.AddMult( Alpha, Pk, p_X );      // p_X = p_X + Alpha*Pk
//      Rk.AddMult( -Alpha, Tmp, Rk );      // Rkp1 = Rk - Alpha*APk

//      double Rkp1Norm2 = Rk.Norm2();
//      if( sqrt(Rkp1Norm2) < p_EMax )
//          break;

//      // Beta
//      double Beta = Rkp1Norm2 / RkNorm2;
//      RkNorm2 = Rkp1Norm2;

//      Rk.AddMult( Beta, Pk, Pk );         // Pkp1 = Rkp1 + Beta*Pk

//      iStep++;
//  }

//  // Log final result
//  double lErr = sqrt( Rk.Norm2() );
//  LOGMSG("[CMatrix::ConjGrad_Solve] : Error= %f, after %d iterations.",lErr,iStep);

//CGraphTools::LogPerfMeter("CMatrix::ConjGrad_Solve");
//  return TRUE;
////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Implementation BOF !!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Implementation BOF !!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//}
#pragma endregion


#pragma region Debug
//================================================================================
bool Matrix::IsEqual( const Matrix & pOther, double pTol/*=-1*/ )
{
    // Check
    if( Rows() != pOther.Rows() || Cols() != pOther.Cols() )
        LzLogException("", "Invalid dimensions!")

    for( unsigned int i = 0 ; i < Rows() ; i++ )
        for( unsigned int j = 0 ; j < Cols() ; j++ )
        {
            if( !ToolBox::IsZero( Elt( i, j ) - pOther.Elt( i, j ), pTol ) )
            {
//#ifdef USE_CLI
//                LzLogMsg("", "Matrices differ: this[" + i + "," + j + "]=" + Elt( i, j ) + ", other[" + i + "," + j + "]=" + pOther.Elt( i, j ) + "." );
//#else
                LzLogMsg("", "Matrices differ: this[" <<
                          i << "," << j << "]=" << Elt( i, j ) <<
                          ", other[" << i << "," << j << "]=" <<
                          pOther.Elt( i, j ) << ".")
//#endif
                return false;
            }
        }

    return true;
}

//================================================================================
//BOOL CMatrix::IsIdentity( double tol/*=-1*/ )
//{
//  if( !Rows() || Rows()!=Cols() )
//  {
//      LOGERR("[CMatrix::IsIdentity] : Matrix is empty or not square !");
//      return FALSE;
//  }

//  for( int i=0 ; i<Rows() ; i++ )
//  for( int j=0 ; j<Cols() ; j++ )
//  {
//      if( fabs(Elt(i,j)-(i==j?1:0)) >= tol )
//          return FALSE;
//  }

//  return TRUE;
//}

//================================================================================
bool Matrix::IsZero( double pTol/*=-1*/ )
{
	// Check
    if( !Rows() || !Cols() )
        LzLogException("","Invalid dimensions: "<<Rows()<<"x"<<Cols()<<"!")

    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
    {
        if( !LzMath::ToolBox::IsZero( Elt( i, j ), pTol ) )
        {
            LzLogMsg("","Matrix is not zero: Elt("<<i<<","<<j<<")="<<Elt(i,j)<<".")
			return false;
        }
    }

    return true;
}

//================================================================================
bool Matrix::IsIdentity( double pTol/*=-1*/ )
{
	// Check
    if( !Rows() || !Cols() )
        LzLogException("","Invalid dimensions: "<<Rows()<<"x"<<Cols()<<"!");

    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
    {
        if( i != j )
		{
            if( !LzMath::ToolBox::IsZero(Elt(i,j), pTol) )
			{
                LzLogMsg("","Matrix is not identity: Elt("<<i<<","<<j<<")="<<Elt( i, j )<<".")
				return false;
			}
		}
		else
		{
            if( !LzMath::ToolBox::IsZero(Elt(i,j) - 1, pTol) )
			{
                LzLogMsg("","Matrix is not identity: Elt("<<i<<","<<j<<")="<<Elt( i, j )<<".")
				return false;
			}
		}
    }

    return true;
}

//================================================================================
void Matrix::LoadRandom( bool pIntValues/*=true*/, double pCenter/*=0*/, double pAmplitude/*=100*/ )
{
    // Check
    if( !Rows() || !Cols() )
        LzLogException("", "Invalid dimensions!")

	// Create uniform distribution
	std::uniform_real_distribution<double> lDist( pCenter-0.5*pAmplitude, pCenter+0.5*pAmplitude );

    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
    {
        double lR = lDist(LzServices::RandomEngine());

        if( pIntValues )
            Elt( i, j ) = ( int )lR;
        else
            Elt( i, j ) = lR;
    }
}

//================================================================================
//void CMatrix::LoadTridiagonal()
//{
//  if( m_H != m_W )
//  {
//      LOGMSG("[CMatrix::LoadTridiagonal] : Matrix not square !");
//      return;
//  }
//
//  srand( (unsigned)time( NULL ) );
//
//  LoadValue(0);
//
//  for( int i=0 ; i<m_H ; i++ )
//  {
//      if( i-1 >= 0 )
//          Elt(i,i-1) = 5;
//
//      Elt(i,i) = 6;
//
//      if( i+1 < m_W )
//          Elt(i,i+1) = 5;
////
////        if( i-1 >= 0 )
////            Elt(i,i-1) = (int)( 20*(double)rand()/(double)RAND_MAX - 10 );
////
////        Elt(i,i) = (int)( 20*(double)rand()/(double)RAND_MAX - 10 );
////
////        if( i+1 < m_W )
////            Elt(i,i+1) = (int)( 20*(double)rand()/(double)RAND_MAX - 10 );
//  }
//}

//================================================================================
//void CMatrix::LoadRandomEmpty()
//{
//  srand( (unsigned)time( NULL ) );

//  for( int i=0 ; i<m_H ; i++ )
//  for( int j=0 ; j<m_W ; j++ )
//  {
//      if( (int)( 4*(double)rand()/(double)RAND_MAX ) == 0 )
//          Elt(i,j) = (int)( 10*(double)rand()/(double)RAND_MAX - 5 );
//      else
//          Elt(i,j) = 0;
//  }
//}

//================================================================================
void Matrix::Log( const std::string & pStd_Prefix/*=""*/, double pRoundOutTol/*=-1*/ ) const
{
 #ifdef USE_CLI
	String ^pPrefix = gcnew String( pStd_Prefix.c_str() );

   if( !Rows() || !Cols() )
    {
        // Header
        LzLogMsg("", pStd_Prefix << ": " << Rows() << "x" << Cols() << " = EMPTY MATRIX" );
    }
    else
    {
        // Header
        LzLogMsg("", pStd_Prefix << ": " << Rows() << "x" << Cols() );

        // Format
        int lLenH = LzMath::ToolBox::StringLength( Rows() - 1 );

        // Variable column width
        Matrix lColWidth( 1, Cols() );
        for( unsigned int j = 0 ; j < Cols() ; j++ )
        {
            lColWidth.Elt( 0, j ) = -1;

            for( unsigned int i = 0 ; i < Rows() ; i++ )
            {
                double E = Elt( i, j );

                // Round out zeros
                if( LzMath::ToolBox::IsZero( E, pRoundOutTol ) )
                    E = 0;

                // Round out ones
                if( LzMath::ToolBox::IsZero( E - 1, pRoundOutTol ) )
                    E = 1;

                int lLen = LzMath::ToolBox::StringLength( E );
                if( lColWidth.Elt( 0, j ) < lLen )
                    lColWidth.Elt( 0, j ) = lLen;
            }
        }

        // Log matrix
        for( unsigned int i = 0 ; i < mH ; i++ )
        {
			String ^lRow = String::Format( "{0}[{1," + lLenH + "}] = |", pPrefix, i );
            for( unsigned int j = 0 ; j < mW ; j++ )
            {
                double E = Elt( i, j );

                // Round out zeros
                if( LzMath::ToolBox::IsZero( E, pRoundOutTol ) )
                    E = 0;

                // Round out ones
                if( LzMath::ToolBox::IsZero( E - 1, pRoundOutTol ) )
                    E = 1;

                lRow += String::Format( " {0," + lColWidth.Elt( 0, j ) + "}", E );
            }

            lRow += " |";

            TxLogM( TxServices::StdStrFromCLI(lRow) );
        }
    }

    // Footer
    TxLogM( "" );
#else
   if( !Rows() || !Cols() )
    {
        // Header
        LzLogMsg("", pStd_Prefix << ": " << Rows() << "x" << Cols() << " = EMPTY MATRIX" );
    }
    else
    {
        // Header
        LzLogMsg("", pStd_Prefix << ": " << Rows() << "x" << Cols() );

        // Format
        int lLenH = LzMath::ToolBox::StringLength( Rows() - 1 );

        // Variable column width
        Matrix lColWidth( 1, Cols() );
        for( unsigned int j=0 ; j<Cols() ; j++ )
        {
            lColWidth.Elt( 0, j ) = -1;

            for( unsigned int i=0 ; i<Rows() ; i++ )
            {
                double E = Elt( i, j );

                // Round out zeros
                if( LzMath::ToolBox::IsZero( E, pRoundOutTol ) )
                    E = 0;

                // Round out ones
                if( LzMath::ToolBox::IsZero( E - 1, pRoundOutTol ) )
                    E = 1;

                int lLen = LzMath::ToolBox::StringLength( E );
                if( lColWidth.Elt( 0, j ) < lLen )
                    lColWidth.Elt( 0, j ) = lLen;
            }
        }

        // Log matrix
        for( unsigned int i = 0 ; i < mH ; i++ )
        {
			std::stringstream lSS;
			lSS << pStd_Prefix << "[" << std::setw(lLenH) << i << "] = |";

            for( unsigned int j = 0 ; j < mW ; j++ )
            {
                double E = Elt( i, j );

                // Round out zeros
                if( LzMath::ToolBox::IsZero(E, pRoundOutTol) )
                    E = 0;

                // Round out ones
                if( LzMath::ToolBox::IsZero(E-1, pRoundOutTol) )
                    E = 1;

				lSS << " " << std::setw(lColWidth.Elt(0, j)) << E;
            }

			lSS << " |";

            LzLogM("", lSS.str())
        }
    }

    // Footer
    LzLogM("", "")
#endif
}

//================================================================================
std::string Matrix::M3x3_ToString() const
{
	// Check
	if( Rows()!=3 || Cols()!=3 )
        LzLogException("","Matrix is not 3x3!")

	std::stringstream lStream;
    lStream << std::setprecision(15);
    //
	lStream << Elt(0,0) << " " << Elt(0,1) << " " << Elt(0,2);
	lStream << Elt(1,0) << " " << Elt(1,1) << " " << Elt(1,2);
	lStream << Elt(2,0) << " " << Elt(2,1) << " " << Elt(2,2);

	return lStream.str();
}

//================================================================================
std::string Matrix::M4x4_ToString() const
{
	// Check
	if( Rows()!=4 || Cols()!=4 )
        LzLogException("","Matrix is not 4x4!")

	std::stringstream lStream;
    lStream << std::setprecision(15);
    //
    lStream << Elt(0,0) << " " << Elt(0,1) << " " << Elt(0,2) << " " << Elt(0,3) << " ";
	lStream << Elt(1,0) << " " << Elt(1,1) << " " << Elt(1,2) << " " << Elt(1,3) << " ";
	lStream << Elt(2,0) << " " << Elt(2,1) << " " << Elt(2,2) << " " << Elt(2,3) << " ";
	lStream << Elt(3,0) << " " << Elt(3,1) << " " << Elt(3,2) << " " << Elt(3,3);

	return lStream.str();
}

//================================================================================
void Matrix::M4x4_FromString( const string & /*pStd_From*/ )
{
#ifdef USE_CLI
	String ^pFrom = gcnew String( pStd_From.c_str() );

	// Parse
	array<wchar_t> ^lSpaceTab = { ' ', '\t' };
	array<String^> ^lTokens = pFrom->Split( lSpaceTab, System::StringSplitOptions::RemoveEmptyEntries );

	// Check
	if( lTokens->Length != 16 )
        LzLogException("","Could not read 16 elements from string '"<<pStd_From<<"'! (found "<<lTokens->Length<<" token(s))");

	// Set (new) dims
	SetDims( 4 ,4 );


	// Set values
	Elt(0,0) = TxServices::StrToDouble( lTokens[ 0] );
	Elt(0,1) = TxServices::StrToDouble( lTokens[ 1] );
	Elt(0,2) = TxServices::StrToDouble( lTokens[ 2] );
	Elt(0,3) = TxServices::StrToDouble( lTokens[ 3] );

	Elt(1,0) = TxServices::StrToDouble( lTokens[ 4] );
	Elt(1,1) = TxServices::StrToDouble( lTokens[ 5] );
	Elt(1,2) = TxServices::StrToDouble( lTokens[ 6] );
	Elt(1,3) = TxServices::StrToDouble( lTokens[ 7] );

	Elt(2,0) = TxServices::StrToDouble( lTokens[ 8] );
	Elt(2,1) = TxServices::StrToDouble( lTokens[ 9] );
	Elt(2,2) = TxServices::StrToDouble( lTokens[10] );
	Elt(2,3) = TxServices::StrToDouble( lTokens[11] );

	Elt(3,0) = TxServices::StrToDouble( lTokens[12] );
	Elt(3,1) = TxServices::StrToDouble( lTokens[13] );
	Elt(3,2) = TxServices::StrToDouble( lTokens[14] );
	Elt(3,3) = TxServices::StrToDouble( lTokens[15] );
#else
    LzLogException("","TODO: Matrix::M4x4_FromString( const std::string & pStd_From )")
#endif
}

//================================================================================
//void CMatrix::NonNullCellStatPerRow( int & Max, int & Min, int & Avg, int & Tot, double tol/*=-1*/ )
//{
//  Max = 0;
//  Min = Cols();
//  Tot = 0;

//  for( int i=0 ; i<Rows() ; i++ )
//  {
//      int Row = 0;

//      for( int j=0 ; j<Cols() ; j++ )
//          if( fabs(Elt(i,j)) >= tol )
//              Row++;

//      if( Row > Max )
//          Max = Row;

//      if( Row < Min )
//          Min = Row;

//      Tot += Row;
//  }

//  Avg = Tot / Rows();
//}

//================================================================================
//void CMatrix::LogBinToRAWImage( CString p_FnamePrefix/*=""*/ )
//{
//  if( !Rows() )
//  {
//      LOGMSG("[CMatrix::LogBinToRAWImage] : %s [.] = empty matrix - no file written",p_FnamePrefix);
//  }
//  else
//  {
//      CArray2<BYTE> Pic(Rows(),Cols());
//      Pic.SetValue( 255 );

//      for( int I=0 ; I<Rows() ; I++ )
//      for( int J=0 ; J<Cols() ; J++ )
//      {
//          if( fabs( Elt(I,J) ) > 1e-6 )
//              Pic[I][J] = 0;
//      }


//#pragma message("----- LogBinToRAWImage : fwrite directement !!! pas de CArray2")


//      // Save
//      CString FName;
//      FName.Format("%s_Mat_%dx%d.raw",p_FnamePrefix,Rows(),Cols());

//      if( CGraphTools::SaveRaw(FName,Pic.GetBuffer(),Cols(),Rows(),FALSE/*H flip*/) )
//          LOGMSG("[CMatrix::LogBinToRAWImage] : Writing matrix to file <%s>",FName);
//      else
//          LOGERR("[CMatrix::LogBinToRAWImage] : Could not write to file <%s> !",FName);
//  }
//}
#pragma endregion
}
