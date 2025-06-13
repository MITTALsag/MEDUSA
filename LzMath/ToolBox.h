#pragma once

//#include "Complex.h"
#include "StringOperations.h"
#include <LzGeom/Point3D.h>
#include <LzGeom/Vector3D.h>
#include <LzGeom/RigidTr3D.h>
#include <LzGeom/Tree3D.h>
#include <LzServices/Vector.h>
#include <LzServices/LzLog.h>
#include <vector>


namespace LzMath
{
class Matrix;


#pragma region "Constants"
static const double PI = 3.141592653589793238462643383279502884197;
static const double DPI = 2.0 * PI;
static const double LzEpsilon = 1e-12;
static const double RAD_2_DEG = 180.0 / PI;
static const double DEG_2_RAD = PI / 180.0;
#pragma endregion


namespace ToolBox
{
using LzGeom::Point3D;
using LzGeom::Vector3D;
using LzGeom::RigidTr3D;
using LzGeom::Tree3D;
using LzMath::Matrix;
using LzServices::Vector;

using std::vector;
using std::string;
using std::wstring;


#pragma region "Lost and found"
// Fitting a cartesian ellipsoid to a set of points
// Least squares sphere
void FitCartesianEllipsoid( const vector<Point3D> & pPts, Point3D & pCenter, double & pRadX, double & pRadY, double & pRadZ );
void UnitTest_FitCartesianEllipsoid();

// Iterative dichotomic solving of a monotonic function
double DichotSolve( double pMin, double pMax, std::function<double(double)> pFunc, double pY, double pYTol=1e-6, double pXTol=1e-6 );

// Perimeter of an ellpsoid, implementing Ramanujan's approximation
double EllipsoidPerimeter( double a, double b );

// Find an ellipsoid arc given the geometric constraints.
// /§\ pA or pB can be negative if XDir is not oriented inwards to the ellipse.
void ComputeEllArc( const Point3D & pXPnt, const Vector3D & pXDir_In, const Point3D & pTPnt, const Vector3D & pTDir,
					Point3D & pCenter, Vector3D & pXDir, Vector3D & pYDir, double & pA, double & pB, double & pRad );

// Find all intersections between an edge and a hyper ellipsoid
// The ell is centered in 0 and aligned along cartesian axes, so A and B must be properly rotated and translated
// before calling the function
void FindHyperEdgeInterEll( const Point3D & pA, const Point3D & pB, const Vector<double> & pABCs, List<Point3D> & pInters );
void FindHyperEdgeInterEll( const Matrix & pA, const Matrix & pB, const Vector<double> & pABCs, List<Matrix> & pInters );
void UnitTest_FindHyperEdgeInterEll( double & pMaxErr );

// Find nearest point on a hyper ellipsoid, given a metric tensor.
//
// The ellipsoid is centered in 0 and aligned along cartesian axes, so A and B must be properly rotated and translated
// before calling the function.
//
// If ppNear or ppFar are nullptr, then they are not computed. At least one must be non-nullptr.
//
// If ppMetric is set to nullptr then the default Euclidean metric will be used.
//
// Source: ALGORITHMS OF PROJECTION OF A POINT ONTO AN ELLIPSOID, Yu. N. Kiseliov, Lithuanian Mathematical Journal, Vol. 34, No. 2, 1994
void ProjectOnHyperEll( const Matrix & pA, const Vector<double> & pABCs, Matrix * ppNear, Matrix * ppFar, Matrix * ppMetric,
                        size_t pMaxIter=100, double pMaxErr=1e-10 );
void UnitTest_ProjectOnHyperEll( double & pMaxErr );

// System::Double::Epsilon: smallest such that x+DBL_EPSILON != x
//http://en.wikipedia.org/wiki/Machine_epsilon
// Better way of doing this? ==> http://floating-point-gui.de/errors/comparison/
//static const double TxEpsilon = 1e-12; //System::Double::Epsilon;
//************ EPSILON = 1e-16 : A TESTER sur cas concrets (distance entre 2 points 3D etc...) sinon TxEpsilon

// If pTol<0, uses default system tolerance TxEpsilon
bool IsZero( double pX, double pTol = -1 );


// Sigmoid transition function S(0)=0, S(1)=1
extern inline double SigmoidPolyRank3( double t );


#if defined(USING_QT) && !defined(NO_OPENGL)
// Generic color scale
void SetScaleColor( double pMin, double pMax, double pVal );
void SetScaleColor( double pMin, double pMax, double pVal, double pMyAlpha );
#endif


// Quick sort. KEYs are sorted in an increasing order (pKeys[i] <= pKeys[i+1] <= ... <= pKeys[j])
template <class KEY> bool IsSorted( KEY * pKeys, int pCount );
template <class KEY> void QuickSort( KEY * pKeys, int pCount )
{
    QuickSort( pKeys, 0, pCount-1, (char*)nullptr );
}
template <class KEY> void QuickSort( Vector<KEY> & pKeys )
{
    QuickSort( pKeys.Buffer(), 0, pKeys.Size()-1, (char*)nullptr );
}

// Should be hidden by implementation
template <class KEY, class ITEM> void QuickSort( KEY * pKeys, int i, int j, ITEM * pItems );
//template <class KEY, class ITEM> void QuickSort(KEY * pKeys, int i, int j, int pAxis, ITEM * pItems);
template <class KEY, class ITEM> void QuickSort( KEY * pKeys, int i, int j, int pAxis, ITEM * pItems );
    
template <class KEY, class ITEM> void QuickSort( Vector<KEY> & pKeys, Vector<ITEM> & pItems );
template <class KEY, class ITEM> void QuickSort( KEY * pKeys, int pCount, ITEM * pItems )
{
    QuickSort( pKeys, 0, pCount - 1, pItems );
}

// Min max in an array
template <class T> void GetMinMax( T * pArray, size_t pCount, T & pMin, T & pMax );

//**** USE MeanMaxStdDev
//// Mean max RMS of a set of errors
//void MeanMaxRMS( const vector<double> & pErr, double & pMean, double & pMax, double & pRMS );
//void MeanMaxRMS( const Vector<double> & pErr, double & pMean, double & pMax, double & pRMS );
//**** USE MeanMaxStdDev

// Mean max standard deviation of a set of errors
void MeanMaxStdDev( const vector<double> & pErr, double & pMean, double & pMax, double & pStdDev, size_t * ppMaxIdx=nullptr );
void MeanMaxStdDev( const Vector<double> & pErr, double & pMean, double & pMax, double & pStdDev, size_t * ppMaxIdx=nullptr );
void MeanMaxStdDev( const List<double> & pErr, double & pMean, double & pMax, double & pStdDev, size_t * ppMaxIdx=nullptr );

// Mass center
Point3D Centroid( const vector<Point3D> & pPoints );
Point3D Centroid( const Vector<Point3D> & pPoints );
Point3D Centroid( const List<Point3D> & pPoints );
Point3D WeightedCentroid( const Vector<Point3D> & pPoints, const Vector<double> & pWeights );

// Finds loops of nearest vertices between two point clouds. Explores all attractors starting from all point in first dataset.
class FindAttractorsParams
{
public:
    FindAttractorsParams( bool pUseTree=false,
						  bool pLogTree=true,
						  double pMinDiameter=0.001,
                          size_t pMaxLeafCount=10,
                          size_t pMaxTreeDepth=10 )
     : mUseTree(pUseTree), mLogTree(pLogTree), mMinDiameter(pMinDiameter), mMaxLeafCount(pMaxLeafCount), mMaxTreeDepth(pMaxTreeDepth)
    {}

    bool mUseTree;
    bool mLogTree;
    double mMinDiameter;
    size_t mMaxLeafCount;
    size_t mMaxTreeDepth;
};
//void FindAttractors( const vector<Point3D> & pA, const vector<Point3D> & pB, List< List<size_t> > & pAtoB, const FindAttractorsParams & pFAP, const Tree3D * ppTree_A=nullptr, const Tree3D * ppTree_B=nullptr );
//void FindAttractors( const Vector<Point3D> & pA, const Vector<Point3D> & pB, List< List<size_t> > & pAtoB, const FindAttractorsParams & pFAP, const Tree3D * ppTree_A=nullptr, const Tree3D * ppTree_B=nullptr );


#pragma region "Polygons"
// pBoundaryIsIn = DEPRECATED
// Returns true if pPt lies inside the triangle
// or on the boundary of the triangle (Behavior: pBoundaryIsIn = true)
// Reason: to determin whether a point is *on* a boundary requires a tolerance, and makes things complicated... (we like simplicity)
//
// Warning: this function throws an exception if the passed triangle is flat (withing IsZero tolerance)!
bool PointIsInTriangle( const Point3D & pPt, const Point3D & pA, const Point3D & pB, const Point3D & pC/*, bool pBoundaryIsIn*/ );

// Returns true if pPt lies inside the polygon pPoly
//
// pPoly : Closed polygon.
//          No need to repeat first point at the end of the vector.
//          Ex.: triangle is defined by 3 vertices = A, B, C (no need to repeat A here)
//
//************* deprecated ==> use new version
//bool PointIsInPolygon( const Point3D & pPt, const vector<Point3D> & pPoly, bool pBoundaryIsIn );
//bool PointIsInPolygon( const Point3D & pPt, const Vector<Point3D> & pPoly, bool pBoundaryIsIn );
enum class PointAndPolygon { In, Out, Surface };
PointAndPolygon PointIsInPolygon( const Point3D & pPt, const std::vector<Point3D> & pPoly );
PointAndPolygon PointIsInPolygon( const Point3D & pPt, const Vector<Point3D> & pPoly );
//************* deprecated ==> use new version
bool SegmentCutsConvexQuad( const Point3D & pA, const Point3D & pB, const vector<Point3D> & pQuad );

// Polygon = closed outline = circular list of points = (P0, P1, P2, ..., PN, P0)
// S = 0.5 * ( X1*Y2 - X2*Y1 + X2*Y3 - X3*Y2 + ... + Xn-1*Yn - Xn*Yn-1 + Xn*Y1 - X1*Yn)
double PolygonArea( const List<Point3D> & pPoly );

//// Returns true if the split is possible
//bool SplitPolygon( const List<Point3D> & pPoly, int pLink1, int pLink2, List<Point3D> & pP1, List<Point3D> & pP2 );
#pragma endregion


//// Size of an array
//template <typename T,size_t S> inline size_t arraysize( const T (&v)[S] ) { return S; }
//template <size_t N> void func(const int (&a)[N]) {}

// Logical rotation in 32 bit
// Uses size_ts "The C specification does not specify if the sign bit is shifted over or not. It is implementation dependent."
unsigned int RotL_32bit( unsigned int pV, int pShift );
unsigned int RotR_32bit( unsigned int pV, int pShift );


#pragma endregion


#pragma region "Piecewise linear function"
//******* Dans Vector: utiliser size_t et pas unsigned int pour indexer des bidules...
//******* Dans Vector: utiliser size_t et pas unsigned int pour indexer des bidules...
//******* Dans Vector: utiliser size_t et pas unsigned int pour indexer des bidules...
//******* Dans Vector: utiliser size_t et pas unsigned int pour indexer des bidules...
//******* Dans Vector: utiliser size_t et pas unsigned int pour indexer des bidules...

//******* Recuperer code de T3D pour lire une piecewise linear func from file
//******* Recuperer code de T3D pour lire une piecewise linear func from file
//******* Recuperer code de T3D pour lire une piecewise linear func from file
//******* Recuperer code de T3D pour lire une piecewise linear func from file
//******* Recuperer code de T3D pour lire une piecewise linear func from file

class PiecewiseLinearFunction
{
public:
    PiecewiseLinearFunction();
    PiecewiseLinearFunction( const Vector<double> & pXs, const Vector<double> & pYs );
    virtual ~PiecewiseLinearFunction();
protected:
    virtual void Free();
    void SetSamples( const Vector<double> & pXs, const Vector<double> & pYs );

public:
    double EvalAt( double pX, bool pClampValue=false ) const;
    const Vector<double> & Xs() const
    {
        return mXs;
    }
    void Scale( double pScale );

public:
    void Log() const;
    std::string ToString() const;
    void FromString( const std::string & pFrom );

protected:
    Vector<double> mXs;
    Vector<double> mYs;
    Vector<double> mAs;
};
#pragma endregion


#pragma region "Eigen values"
#pragma region "Cubic roots"
// Solving a 3rd grade polynomial: a*Z^3 + b*Z^2 + c*Z + d = 0
// Source: http://en.wikipedia.org/wiki/Cubic_function
//void CubicPolyRoots( Vector<Complex> & pRoots, double a, double b, double c, double d );

// Returns the maximal modulus of the 3 residuals
//double CheckCubicPolyRoots( const Vector<Complex> & pRoots, double a, double b, double c, double d );

// Returns an accurate real root
double AccurateRealRoot( double pInitRoot, double pMaxErr, size_t pMaxIter, double a, double b, double c, double d );
//double AccurateRealRoot_Dichot( double pInitRoot, double pMaxErr, size_t pMaxIter, double a, double b, double c, double d );
#pragma endregion


// Jacobi Transformations of a Symmetric Matrix cf. Numerical Recipes in C++ (2007, pp.570-576)
// A = V D Vt
// Where:
//        . D is a Nx1 vector containing eigenvalues
//        . V is a NxN matrix
//
// Columns of matrix V = normalized eigenvectors = (same order as eigenvalues i.e. decreasing order in D)
void EigenSort( Matrix & d, Matrix * v=nullptr );
void EigenJacobi( const Matrix & pA, Matrix & pDiag, Matrix & pVecs );

// Polar decomposition of any matrix into a sym def pos matrix pSym and an orthogonal matrix
// A = Sym Orth
// or
// A = Orth Sym
void PolarDecompose_SymOrth( const Matrix & pA, Matrix & pSym, Matrix & pOrth );
void PolarDecompose_OrthSym( const Matrix & pA, Matrix & pOrth, Matrix & pSym );

// PCA analysis of a points cloud
//
// pV[i] = principal axis of the points cloud = ortho-normalized vectors
void PCA( const List<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] );
void PCA( const Vector<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] );
void PCA( const vector<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] );
void WeightedPCA__KindOfBuggy( const Vector<Point3D> & pCloud, const Vector<double> & pWeights, Point3D & pMean, Vector3D pV[3] );
void WeightedPCA( const List<Point3D> & pCloud, const List<double> & pWeights, Point3D & pMean, Vector3D pV[3] );
void WeightedPCA( const Vector<Point3D> & pCloud, const Vector<double> & pWeights, Point3D & pMean, Vector3D pV[3] );
void PCAMinMax( const vector<Point3D> & pCloud, const Point3D & pO, const Vector3D pV[3], Vector<double> & pMin, Vector<double> & pMax );
void PCAMinMax( const List<Point3D> & pCloud, const Point3D & pO, const Vector3D pV[3], Vector<double> & pMin, Vector<double> & pMax );
void PCAMinMax( const Vector<Point3D> & pCloud, const Point3D & pO, const Vector3D pV[3], Vector<double> & pMin, Vector<double> & pMax );
//
void InertialPointCloudPCA( const vector<Point3D> & pCloud, const vector<double> & pWeights, Point3D & pMean, Vector3D pV[3] );

//*************** voir numerical recipes : difference entre singular values, eigen values - diagonalisation et decomposition SVD
//*************** si matrice carree alors decompo svd = diagonalisation (si diagonalisation unique - a priori oui)

// SVD
// /!\ The eigen values are not sorted at the end of the procedure.
// Singular Value Decomposition of a matrix : A = U W Vt, where W is a column vector containing the singular values
//
// A: mxn, where m >= n
//
// U: mxn
// W: nx1, containing the diagonal of the corresponding nxn matrix
// V: nxn
//
// If A is squared, sym def pos then U == V and W contains the Eigenvalues of A
// that need to be sorted using EigenSort
void SVD( const Matrix & pA, Matrix & pU, Matrix & pW, Matrix & pV );
//--------------------------------------------------------------------------------
//  Code from Numerical recipies.
//  a (mxn) matrix, is modified and returns the matrix U of the decomposition
//  m number of rows in the matrix a
//  n number of columns in the matrix a
//  w (n) vector, returns the singular values in the diagonal matrix W of the decomposition
//  v (nxn) matrix, returns the matrix V (not Vt...) of the decomposition
//  /!\ all indices range in [1,m] and [1,n]
//--------------------------------------------------------------------------------
void SVD( double **a, int m, int n, double *w, double **v );
#pragma endregion


#pragma region "Pressure conversion"
enum class PaUnit : int { Raw = 0, Pa, kPa, Ncm2, kgcm2, gcm2, mmHg, psi };
double ConvertPa( PaUnit pFrom, PaUnit pTo, double pVal );
std::string PaUnitSuffix( PaUnit pUnit );
#pragma endregion
   

#pragma region "Gradient descent"
class GradResult
{
public:
    enum class Status { Converged, StillRunning, EnergyRising, NotComputed };

    GradResult( Status pStatus=Status::NotComputed, size_t pSteps=0, double pGradNorm=-1, double pEnergy=-1 )
        : mStatus(pStatus), mSteps(pSteps), mGradNorm(pGradNorm), mEnergy(pEnergy)
    {}

    string toString() const
    {
        string lRes;
        switch( mStatus )
        {
        case Status::Converged: lRes = "Converged."; break;
        case Status::StillRunning: lRes = "Still running."; break;
        case Status::EnergyRising: lRes = "Energy rising."; break;
        case Status::NotComputed: lRes = "--Not computed--."; break;
        default: LzLogException("", "Unexpected status!")
        }
        return lRes + " Steps= " + std::to_string(mSteps)
                    + ", grad norm= " + std::to_string(mGradNorm)
                    + ", energy= " + std::to_string(mEnergy)
                    + ".";
    }

public:
    Status mStatus;
    size_t mSteps;
    double mGradNorm;
    double mEnergy;
};
//
class GradDomain
{
public:
    // Init
    GradDomain( size_t pDim );

    // Check if dimension is constrained
    bool IsBound( size_t pDim ) const;

    // Checking bounds
    bool NearMin( const Vector<double> & pX, size_t pDim ) const;
    bool NearMax( const Vector<double> & pX, size_t pDim ) const;

    // Check if point is out of bounds
    bool IsOutOfBounds( const Vector<double> & pX ) const;

    // Clamping
    // Assumtion: A is within bounds; B might not be
    // If B is within bounds, returns 1.0
    // If not, returns a squeeze factor in ]0.0, 1.0[
    double Clamp( const Vector<double> & pA, const Vector<double> & pB ) const;

public:
    // Data
    Vector<double> mMins;
    Vector<double> mMaxs;
    double mMinMaxMargin;
};
//
enum class GradLog { All, AllPara, None };
//
GradResult GradDescent( GradLog pLogMode, Vector<double> & pX, std::function<double(const Vector<double> &)> pEnergy,
                        size_t pMaxSteps, double pMinGradNorm, const GradDomain & pDom, double pFiniteDiffStep=1e-5 );
#pragma endregion


#pragma region "Conjugate gradient"
/*
 * An Introduction to the Conjugate Gradient Method Without the Agonizing Pain, Edition 1 1/4
 * Jonathan Richard Shewchuk
 * August 4, 1994
 *
 */
//tester avec systeme lineaire
//tester avec systeme lineaire
//tester avec systeme lineaire
//tester avec systeme lineaire
//tester avec systeme lineaire
//tester avec systeme lineaire

//*** cf doc NLOpt : trust region, COBYLA, constrained conjugate gradient... ???
//
GradResult ConjugateGradient( GradLog pLogMode, Vector<double> & pX, std::function<double(const Vector<double> &)> pEnergy,
                              size_t pMaxSteps, double pMinGradNorm, const GradDomain & pDom, double pFiniteDiffStep=1e-5 );/*,
                              double & pGradErr,
                              double pIncrementdX,*/
                              /*CG_GetGradient ^pGetGradient,*/ /*, CG_PostStep ^pPostStep );*/
// ******************* a virer !!!
//void CG_GetGradient_FinDiff( const Matrix & pX, Matrix & pG )
//{
//  /* compute gradient using finite differences */
//}
// ******************* a virer !!!
#pragma endregion


#pragma region "Registration"
/**
 *  @brief Produce all alias matrices R'(X', Y', Z') to R(X, Y, Z)
 */
void GetCartesianPermutations( vector<RigidTr3D> & pTo,
                               const List<int> & pForbiddenXs={},
                               const List<int> & pForbiddenYs={}/*,
                               const List<int> & pForbiddenZs={}*/ );

// Rigid registration: Arun

// Returns T such as: sum_i( | T*P1_i - P2_i |^2 ) --> least squares min
RigidTr3D Arun( const Vector<Point3D> & pP1, const Vector<Point3D> & pP2 );
RigidTr3D Arun( const vector<Point3D> & pP1, const vector<Point3D> & pP2 );

// Compute error
void LogArunError( const RigidTr3D & pAlibi_P1_to_P2, const Vector<Point3D> & pP1, const Vector<Point3D> & pP2 );

// ICP: Iterative Arun until all Euler rotations are < max and norm of translation is < max
//
//*** NB: points pairing is brute force! for loop across all destination points... not the best implementation.
//
RigidTr3D ArunLoop( const vector<Point3D> & pP1, const vector<Point3D> & pP2, double pMaxEulerRotDeg=0.1, double pMaxTransNorm=0.1, size_t pMaxSteps=300 );
RigidTr3D ArunLoop( const Vector<Point3D> & pP1, const Vector<Point3D> & pP2, double pMaxEulerRotDeg=0.1, double pMaxTransNorm=0.1, size_t pMaxSteps=300 );
RigidTr3D ArunLoop( const vector<Point3D> & pP1, const Tree3D & pDstTree, double pMaxEulerRotDeg=0.1, double pMaxTransNorm=0.1, size_t pMaxSteps=300 );
RigidTr3D ArunLoop( const Vector<Point3D> & pP1, const Tree3D & pDstTree, double pMaxEulerRotDeg=0.1, double pMaxTransNorm=0.1, size_t pMaxSteps=300 );

// Affine registration: scaling along principal axis
Matrix AffinePCAFitPairedPoints( vector<Point3D> & pSrc, const vector<Point3D> & pDst, double pTolToIdentity );
Matrix AffinePCAFitPairedPoints( Vector<Point3D> & pSrc, const Vector<Point3D> & pDst, double pTolToIdentity );

// Performs:
// 1) rigid Arun registration between the paired point sets
// 2) an affine optimal scaling along the principal axis of the source data set
//
// Source points are modified by the operation, destination points are not.
//
// Returns true if BOTH the rigid Arun AND the affine deformations are identity within the specified tolerance.
//
// pSrcV are normalized eigen vectors produced by PCA( const Vector<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] )
//
// Does not recompute PCA axes at each iteration but merely updates them
// Optional parameter ppLockedScales lists the directions (0, 1, or 2) which will have scales fixed to 1.0 during the iteration
bool AffinePCAFitPairedPoints_Iteration( Vector<Point3D> & pSrc, const Vector<Point3D> & pDst, RigidTr3D & pArun, Matrix & pAffine, double pTolToIdentity,
                                         Point3D & pSrcG, Vector3D pSrcV[3], List<size_t> * ppLockedScales=nullptr );

//// ICP iteration ==> Arun
//RigidTr3D ICP_Iteration( const Vector<Point3D> & pSrc, const Vector<Point3D> & pDst );

// Affine scale de la bbox sur base PCA des vecteurs (to nearest)
// Affine scale de la bbox sur base PCA des vecteurs (to nearest)
// Affine scale de la bbox sur base PCA des vecteurs (to nearest)
#pragma endregion


#pragma region "Pivot calibration"
void BuildPivot6D( const Vector<RigidTr3D> & pTrRefA2RefB, Point3D & pPivotInRefA, Point3D & pPivotInRefB, double pVar[3] );
#pragma endregion
}
}


#pragma region "Templates implementation"
//================================================================================
template <class KEY> bool LzMath::ToolBox::IsSorted( KEY * pKeys, int pCount )
{
    for( int i = 1 ; i < pCount ; i++ )
    {
        if( pKeys[i - 1] > pKeys[i] )
        {
            LzLogM("", "(Key[" << ( i - 1 ) << "]=" << pKeys[i - 1] << ") > (Key[" << i << "]=" << pKeys[i] << ")")
            return false;
        }
    }

    return true;
}

//================================================================================
template <class KEY, class ITEM> int Partition( KEY * pKeys, int p, int r, int pAxis, ITEM * pItems )
{
    double x = pKeys[p].mV[pAxis];

    int i = p - 1;
    int j = r + 1;
    //==> i <= j

    while( true )
    {
        do
        {
            i++;
        }
        while( pKeys[i].mV[pAxis] < x );   //==> pKeys[i].mV[pAxis] >= x

        do
        {
            j--;
        }
        while( pKeys[j].mV[pAxis] > x );   //==> pKeys[j].mV[pAxis] <= x

        if( i < j )
        {
            // Swap keys
            KEY tmpKey = pKeys[i];
            pKeys[i] = pKeys[j];
            pKeys[j] = tmpKey;

            // Swap items (if any)
            if ( pItems )
            {
                ITEM tmpItem = pItems[i];
                pItems[i] = pItems[j];
                pItems[j] = tmpItem;
            }
        }
        else
            return j;
    }
}

//================================================================================
template <class KEY, class ITEM> int Partition( KEY * pKeys, int p, int r, ITEM * pItems )
{
    KEY x = pKeys[p];

    int i = p - 1;
    int j = r + 1;
    //==> i <= j

    while( true )
    {
        do
        {
            i++;
        }
        while ( pKeys[i] < x );   //==> pKeys[i] >= x

        do
        {
            j--;
        }
        while ( pKeys[j] > x );   //==> pKeys[j] <= x

        if( i < j )
        {
            // Swap keys
            KEY tmpKey = pKeys[i];
            pKeys[i] = pKeys[j];
            pKeys[j] = tmpKey;

            // Swap items (if any)
            if( pItems )
            {
                ITEM tmpItem = pItems[i];
                pItems[i] = pItems[j];
                pItems[j] = tmpItem;
            }
        }
        else
            return j;
    }
}

//================================================================================
class IJ
{
public:
    IJ( int pI, int pJ ) : mI( pI ), mJ( pJ ) {}

    int mI;
    int mJ;
};

//================================================================================
template <class KEY, class ITEM> void LzMath::ToolBox::QuickSort( KEY * pKeys, int i, int j, ITEM * pItems )
{
    // Stacks
    List<IJ> lStack;

    // Push
    lStack.AddTail( IJ( i, j ) );

    // Sort
    while( lStack.Count() )
    {
        // Pop
        const IJ & lPopIJ = lStack.GetHead();
        int i = lPopIJ.mI;
        int j = lPopIJ.mJ;
        lStack.DelHead();

        if( i < j )
        {
            int q = Partition( pKeys, i, j, pItems );

            // Push: QuickSort( pKeys, i,   q, pItems )
            lStack.AddTail( IJ( i, q ) );

            // Push: QuickSort( pKeys, q+1, j, pItems )
            lStack.AddTail( IJ( q + 1, j ) );
        }
    }
}

// Recursivity ==> stack overflow
//================================================================================
//template <class KEY, class ITEM> void LzMath::ToolBox::QuickSort( KEY * pKeys, int i, int j, ITEM * pItems )
//{
//  if( i < j )
//  {
//      int q = Partition( pKeys, i, j, pItems );
//
//      QuickSort( pKeys, i,   q, pItems );
//      QuickSort( pKeys, q+1, j, pItems );
//  }
//}
// Recursivity ==> stack overflow

//================================================================================
template <class KEY, class ITEM> void LzMath::ToolBox::QuickSort( KEY * pKeys, int i, int j, int pAxis, ITEM * pItems )
{
    // Stacks
    List<IJ> lStack;

    // Push
    lStack.AddTail( IJ( i, j ) );

    // Sort
    while( lStack.Count() )
    {
        // Pop
        const IJ & lPopIJ = lStack.GetHead();
        int i = lPopIJ.mI;
        int j = lPopIJ.mJ;
        lStack.DelHead();

        if( i < j )
        {
            int q = Partition( pKeys, i, j, pAxis, pItems );

            // Push: QuickSort( pKeys, i,   q, pAxis, pItems )
            lStack.AddTail( IJ( i, q ) );

            // Push: QuickSort( pKeys, q+1, j, pAxis, pItems )
            lStack.AddTail( IJ( q + 1, j ) );
        }
    }
}

// Recursivity ==> stack overflow
//================================================================================
//template <class KEY, class ITEM> void LzMath::ToolBox::QuickSort( KEY * pKeys, int i, int j, int pAxis, ITEM * pItems )
//{
//  if( i < j )
//  {
//      int q = Partition( pKeys, i, j, pAxis, pItems );
//
//      QuickSort( pKeys, i,   q, pAxis, pItems );
//      QuickSort( pKeys, q+1, j, pAxis, pItems );
//  }
//}

//================================================================================
template <class KEY, class ITEM> void LzMath::ToolBox::QuickSort( Vector<KEY> & pKeys, Vector<ITEM> & pItems )
{
    // Check
    if( pKeys.Size() != pItems.Size() )
        LzLogException("",  "Keys and items do not have the same size! (" << pKeys.Size() << " != " << pItems.Size() << ")")

    // Do the job
    QuickSort( pKeys.Buffer(), pKeys.Size(), pItems.Buffer() );
}

//================================================================================
template <class T> void LzMath::ToolBox::GetMinMax( T * pArray, size_t pCount, T & pMin, T & pMax )
{
	// Check
	if( !pCount )
        LzLogException("",  "Empty array!")

	pMin = pMax = pArray[0];
    for( size_t i = 1 ; i < pCount ; i++ )
    {
        double t = pArray[i];
        if( pMin > t ) pMin = t;
        if( pMax < t ) pMax = t;
    }
}    
#pragma endregion
