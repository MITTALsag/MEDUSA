#pragma once

#include <LzServices/List.h>
#include <string>


namespace LzMath
{
using LzServices::List;
using std::string;

//class PolyMatrix;


class Matrix
{
#pragma region "Construction/destruction"
public:
    // /!\ Warning: use 'explicit' to avoidthe following possible usage: "Matrix M = 1;"
    explicit Matrix( unsigned int pH = 0, unsigned int pW = 1 );
    Matrix( unsigned int pH, unsigned int pW, float const * pRowByRowElts );
    Matrix( unsigned int pH, unsigned int pW, double const * pRowByRowElts );
//#ifdef USE_CLI
    //Matrix( unsigned int pH, unsigned int pW, String ^pRowByRowElts );
    Matrix( unsigned int pH, unsigned int pW, const string & pRowByRowElts );
    Matrix( unsigned int pH, unsigned int pW, const std::initializer_list<double> & pRowByRowElts );
//#endif
    Matrix( const Matrix & pM );
    virtual ~Matrix();
    #pragma endregion


#pragma region "Init"
public:
    virtual void SetDims( unsigned int pH/*Rows*/, unsigned int pW = 1/*Cols*/ );
    virtual void Free();
    static Matrix Identity( unsigned int pH, unsigned int pW );
#pragma endregion


#pragma region "Get/set"
public:
    unsigned int Rows() const { return mH; }
    unsigned int Cols() const { return mW; }

    // 'Elt' methods implement internal storage order. Matrix data should never be accessed by reading mMat directly!
    double & Elt( unsigned int pI, unsigned int pJ = 0 );
    double Elt( unsigned int pI, unsigned int pJ = 0 ) const;
//const double * RowByRow_Data() const { return mMat; }
const double * ColMajorFrom( unsigned int pJ ) const;
#pragma endregion


#pragma region "Operations"
public:
    void LoadIdentity();
	void LoadDiagonal( const Matrix & pDiagInCol, unsigned int pColIdx=0 );
    void LoadValue( double pVal );
    void Mult( const Matrix & pA, Matrix & pRes ) const;
    void SymmetricalMult( const Matrix & pA, Matrix & pRes ) const;
    void SymmetricalMult_MtM( Matrix & pRes ) const;
    void SymmetricalMult_MtM_upper_tri_only( Matrix & pRes ) const;
//    void Mult( const PolyMatrix & pA, PolyMatrix & pRes ) const;
    void Mult( double pA, Matrix & pRes ) const;
    void Mult( double pA );
    void AddMult( double pA, const Matrix & pB, Matrix & pRes ) const;
    void Add( const Matrix & pA, Matrix & pRes ) const;
    void Add( const Matrix & pA );
    void Sub( const Matrix & pA, Matrix & pRes ) const;
    //  BOOL Pow( int pPow, CMatrix & pRes ) const;
    void SquareTranspose();
    void TransposeTo( Matrix & pRes ) const;
    bool IsSymmetric( double pTol = -1 ) const;
    virtual const Matrix & operator=( const Matrix & pIn );
    //  double MaxFabs() const;
    double DotProd( const Matrix & pA ) const;
    double MaxAbs() const;
    double Norm2() const;
    double M2x2_Det() const;
    double M3x3_Det() const;
    void M3x3_EigenPolynom( double & p3, double & p2, double & p1, double & p0 ) const;
    double M4x4_Det() const;
    void M2x2_InverseTo( Matrix & pInv ) const;
    bool M2x2_InverseTo_NoExc( Matrix & pInv ) const;
    void M3x3_InverseTo( Matrix & pInv ) const;
    void M4x4_InverseTo( Matrix & pInv ) const;
    void SubMatrixTo( unsigned int pI0, unsigned int pJ0, Matrix & pTo ) const;
    void SubMatrixFrom( unsigned int pI0, unsigned int pJ0, const Matrix & pFrom );
    //  // Permutations
    ////    BOOL LoadPermutation( const CContainer1<int> & pMapping );
    //
    //  // Rotations
    //  BOOL LoadRotation( const CVector & pRoot, const CVector & pDir, double pDeg );
    //  // Gaussian elimination: makes all equations pVar-independent (except row #pVar)
    //  BOOL SubsInAllRowsFor( int pVar, CVector & pF );
    //
    //  // Produces a matrix with reduced dimension by deleting listed rows and cols
    //  BOOL KillRowCols( const List<int> & pRows, const List<int> & pCols, CMatrix & pToMat ) const;
#pragma endregion


#pragma region "Load/save"
public:
      void SaveBin( const string & pFileName ) const;
      void LoadBin( const string & pFileName );
#pragma endregion


#pragma region "Co-factors"
public:
    double CF_Det( /*virer cette option : par defaut true!*/ bool pEmptyMat = false ) const;
    void CF_InverseTo( Matrix & pInv, bool pEmptyMat = false );
protected:
    int MaxZeroLocalRow() const;
    double CF_RecDet( bool pEmptyMat = false ) const;
    void InitLocalRowsCols() const;
    static List<int> mLocalRows;
    static List<int> mLocalCols;
    int LocalCols() const
    {
        return static_cast<int>( mLocalCols.Count() );
    }
    int RemoveLocalRow( int pLocalI ) const;
    int RemoveLocalCol( int pLocalJ ) const;
    void AddGlobalRow( int pGlobalI ) const;
    void AddGlobalCol( int pGlobalJ ) const;
#pragma endregion


#pragma region "LU decomposition"
public:
    void LU_ify( Matrix & pL, Matrix & pU, double pTol );
    double LU_Det();
    double TriDet();
    void L_InverseTo( Matrix & pInvL );
    void U_InverseTo( Matrix & pInvU );
    void LU_InverseTo( Matrix & pInv );
    //  BOOL LU_Solve( CVector & pX, const CVector & pY );
#pragma endregion


#pragma region "Cholesky LLt: sym def pos matrix"
public:
    void LLt_ify( Matrix & pL );
    ////===>    BOOL LLt_Solve( CVector & pX, const CVector & pY ) const;
    static void LLt_Solve( const Matrix & pL, Matrix & pX, const Matrix & pY );
#pragma endregion


#pragma region "Iterative solver: gradient and conjugate gradient"
public:
    void SolveSteepestDescent( Matrix & pU, const Matrix & pF, unsigned int pMaxIterations, double pErrMaxAbs ) const;
    void SolveConjGrad( Matrix & pX, const Matrix & pY, unsigned int pMaxIterations, double pErrMaxAbs ) const;
    void SolvePrecondConjGrads( Matrix & pU, const Matrix & pF, unsigned int pMaxIterations, double pErrMaxAbs ) const;
#pragma endregion


#pragma region "Debug"
public:
    bool IsEqual( const Matrix & pOther, double pTol = -1 );
    //  BOOL IsIdentity( double pTol=-1 );
    bool IsZero( double pTol = -1 );
    bool IsIdentity( double pTol = -1 );
    void LoadRandom( bool pIntValues = true, double pCenter = 0, double pAmplitude = 100 );
//  void LoadTridiagonal();
//  void LoadRandomEmpty();
    string M3x3_ToString() const;
    string M4x4_ToString() const;
    void M4x4_FromString( const string & pStd_From/*String ^pFrom*/ );
    void Log( const string & pPrefix=""/*String ^pPrefix = ""*/, double pRoundOutTol = -1 ) const;
//  void NonNullCellStatPerRow( int & Max, int & Min, int & Avg, int & Tot, double pTol=-1 );
//  void LogBinToRAWImage( CString pFnamePrefix="" );
#pragma endregion


protected:
    unsigned int mH;
    unsigned int mW;
    double * mMat;
};
}
