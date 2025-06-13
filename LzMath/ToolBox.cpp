#include "ToolBox.h"
#include "Matrix.h"
#include <LzGeom/Line3D.h>
#include <LzGeom/Plane3D.h>
#include <LzServices/LzLog.h>
#include <limits>
#if defined(USING_QT) && !defined(NO_OPENGL)
    #include <QtOpenGL> // Need to include this before including glu.h otherwise MSVC does not compile
    #include <GL/gl.h>
#endif
#include <stdlib.h>
#include <cmath>


//*************** conflicting with std::numeric_limits<double>::max()
#undef max
//*************** conflicting with std::numeric_limits<double>::max()

namespace LzMath
{
namespace ToolBox
{
using LzGeom::Line3D;
using LzGeom::Plane3D;


#pragma region "Lost and found"
//================================================================================
double DichotSolve( double pMin, double pMax, std::function<double(double)> pFunc, double pY, double pYTol/*=1e-6*/, double pXTol/*=1e-6*/ )
{
    // Check
    if( pMin >= pMax )
        LzLogException("", "Invalid bounds! Min= "<<pMin<<", max= "<<pMax<<".")

    // Determine function trend
    bool lFunc_increases = pFunc(pMax) > pFunc(pMin);

//LzLogM("", lFunc_increases ? "INCREASING" : "DECREASING" );

    // Solve pFunc( lX ) = pY
    while( true )
    {
        // Evaluate guess
        double lX = 0.5*(pMin + pMax);
        double lFunc_of_X = pFunc( lX );

//LzLogM("", "X= "<<lX<<", lFunc_of_X= "<<lFunc_of_X);

        // Found solution?
        if( fabs(lFunc_of_X - pY) <= pYTol )
            return lX;

        // Interval collapse
        if( pMax - pMin <= pXTol )
            LzLogException("", "Search interval collapsed below tolerance and no solution was found! Min= "<<pMin<<", max= "<<pMax<<".")

        // Move bounds
        if( pY < lFunc_of_X )
        {
            if( lFunc_increases )
//{
//LzLogM("", "Moving MAX down");
                pMax = lX;
//}
            else
                pMin = lX;
        }
        else
        {
            if( lFunc_increases )
//{
//LzLogM("", "Moving MIN up");
                pMin = lX;
//}
            else
                pMax = lX;
        }

    }
}

//================================================================================
double EllipsoidPerimeter( double a, double b )
{
	// Ramanujan
    return PI*( 3*(a + b) - sqrt((3*a + b)*(a + 3*b)) );
}

//================================================================================
void ComputeEllArc( const Point3D & pXPnt, const Vector3D & pXDir_In, const Point3D & pTPnt, const Vector3D & pTDir,
					Point3D & pCenter, Vector3D & pXDir, Vector3D & pYDir, double & pA, double & pB, double & pRad )
{
	// Set local transform
	RigidTr3D lAlias_Ell2W;
	lAlias_Ell2W.SetOrthonormal( pXDir_In, pTPnt-pXPnt, pXPnt );

	// Correct X and Y dir
	pXDir = lAlias_Ell2W * Vector3D(1, 0, 0);
	pYDir = lAlias_Ell2W * Vector3D(0, 1, 0);

	// Compute slope in local ref
	const double lTan_X = pTDir * pXDir;
	const double lTan_Y = pTDir * pYDir;

	// Check
	if( fabs(lTan_X) < 1e-6 )
        LzLogException("", "Cannot compute local slope at tangent point!");

	// Slope
	const double p = lTan_Y / lTan_X;

	// Compute tan point in Ref Ell
	const RigidTr3D lAlias_W2lEl = lAlias_Ell2W.Inverse();
	const Point3D lTPnt_Ell = lAlias_W2lEl * pTPnt;
	const double x = lTPnt_Ell.X();
	const double y = lTPnt_Ell.Y();

	// Check
	// y must be > 2px
	if( y <= 2*p*x )
        LzLogException("", "Cannot compute ellipsoid arc! y must be > 2px; x= "<<x<<", p= "<<p<<"; 2*p*x= "<<2*p*x<<", y= "<<y<<".");

	// Ell_Y scale factor
	const double k = sqrt( x*x / (y*y - 2*p*x*y) );

	// Ellipsoid axes
	pA = x + k*k*p*y;
	pB = pA / k;

	// Ellipsoid center
	pCenter = pXPnt + pA*pXDir;

	// Angular position of tan point
    // eq. to: pRad = LzMath::PI - atan2( y/pB, (x - pA)/pA );
	pRad = atan2( y/pB, (pA - x)/pA );
}

//================================================================================
void FindHyperEdgeInterEll( const Matrix & pA, const Matrix & pB, const Vector<double> & pABCs, List<Matrix> & pInters )
{
	// Determine in which dimension we are working
    const size_t DIM = pABCs.Size();

	// Check
	if( DIM == 0 )
        LzLogException("", "Invalid dimension! Expected at least 1, found "<<pABCs.Size()<<".")

	// Check
	if( pA.Rows()!=DIM || pA.Cols()!=1 )
        LzLogException("", "Invalid hyperpoint A! Expected "<<DIM<<"x1, found "<<pA.Rows()<<"x"<<pA.Cols()<<".")

	// Check
	if( pB.Rows()!=DIM || pB.Cols()!=1 )
        LzLogException("", "Invalid hyperpoint B! Expected "<<DIM<<"x1, found "<<pB.Rows()<<"x"<<pB.Cols()<<".")

	// Clean previous
	pInters.DelAll();

	//
	// Inters = A + t(B - A). Now we must find t. t is solution of a quadratic eq: there can be 0, 1 or 2 inters.
	//

	//-----------------------------------
	// Compute terms
	//-----------------------------------

	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	{
		// Squared axes
		Vector<double> lSqABCs( DIM );
		Vector<double> lB_minus_A( DIM );
        for( size_t d=0 ; d<DIM ; d++ )
		{
			lSqABCs[d] = pABCs[d] * pABCs[d];
			lB_minus_A[d] = pB.Elt(d) - pA.Elt(d);
		}

		// a
        for( size_t d=0 ; d<DIM ; d++ )
			a += lB_minus_A[d] * lB_minus_A[d] / lSqABCs[d];

		// b
        for( size_t d=0 ; d<DIM ; d++ )
			b += pA.Elt(d) * lB_minus_A[d] / lSqABCs[d];
		b *= 2;

		// c
        for( size_t d=0 ; d<DIM ; d++ )
			c += pA.Elt(d) * pA.Elt(d) / lSqABCs[d];
		c -= 1.0;
	}

	//-----------------------------------
	// Solve a*t*t + b*t + c = 0
	//-----------------------------------

	// Determinant
	const double Det = b*b - 4*a*c;

	// No solution
	if( Det < 0 )
		return;
	
	if( IsZero(Det) )
	{
		// Root
		double t = -b / (2.0 * a);

		// Keep only if within bounds
		if( t>=0.0 && t<=1.0 )
		{
			// Find inter
			Matrix I( DIM );
            for( size_t d=0 ; d<DIM ; d++ )
				I.Elt(d) = pA.Elt(d) + t*(pB.Elt(d) - pA.Elt(d));

			// Stash
			pInters.AddTail( I );
		}
	}
	else
	{
		// Roots
		double rootDet = sqrt( Det );
		double t[2] =
		{
			(-b - rootDet) / (2.0 * a),
			(-b + rootDet) / (2.0 * a)
		};

		// Find inters
		for( int i=0 ; i<2 ; i++ )
		{
			// Skip if out of bounds
			if( t[i]<0.0 || t[i]>1.0 )
				continue;

			// Find inter
			Matrix I( DIM );
            for( size_t d=0 ; d<DIM ; d++ )
				I.Elt(d) = pA.Elt(d) + t[i]*(pB.Elt(d) - pA.Elt(d));

			// Stash
			pInters.AddTail( I );

//Matrix I0( DIM ), I1( DIM );
//for( size_t d=0 ; d<DIM ; d++ )
//{
//	I0.Elt(d) = pA.Elt(d) + t0*(pB.Elt(d) - pA.Elt(d));
//	I1.Elt(d) = pA.Elt(d) + t1*(pB.Elt(d) - pA.Elt(d));
//}

//// Stash
//pInters.AddTail( I0 );
//pInters.AddTail( I1 );
		}
	}
}

//================================================================================
void FindHyperEdgeInterEll( const Point3D & pA, const Point3D & pB, const Vector<double> & pABCs, List<Point3D> & pInters )
{
	// Convert to matrices
	Matrix lA( 3 ), lB( 3 );
	for( int i=0 ; i<3 ; i++ )
	{
		lA.Elt(i) = pA.mV[i];
		lB.Elt(i) = pB.mV[i];
	}

	// Find inters
	List<Matrix> lInters;
	FindHyperEdgeInterEll( lA, lB, pABCs, lInters );

	// Convert to 3D points
	pInters.DelAll();
	BrowseList( iI, lInters )
	{
		const Matrix & lInter = lInters.GetAt( iI );
		pInters.AddTail( Point3D(lInter.Elt(0), lInter.Elt(1), lInter.Elt(2)) );
	}
}

//================================================================================
void UnitTest_FindHyperEdgeInterEll( double & pMaxErr )
{
	// Pick a random dimension
//const size_t DIM = std::uniform_int_distribution<size_t>(1, 20)( LzServices::RandomEngine() );
const size_t DIM = std::uniform_int_distribution<size_t>(1, 10)( LzServices::RandomEngine() );

    LzLogNode("", "Testing LzMath::ToolBox::FindHyperEdgeInterEll in dimension "<<DIM<<".");

	// Set-up data holders
	std::uniform_real_distribution<double> lRnd(-1000, +1000);
	Vector<double> lABCs( DIM );
	Matrix lA( DIM ), lB( DIM );
    size_t test_count = 0;
	do
	{
		// Set values
        for( size_t d=0 ; d<DIM ; d++ )
		{
            lABCs[d] = 1.0/*to avoid degenerate ellipsoids*/ + fabs( lRnd(LzServices::RandomEngine()) ) * 0.5;

            lA.Elt(d) = lRnd(LzServices::RandomEngine());
            lB.Elt(d) = lRnd(LzServices::RandomEngine());
		}

		// Compute inters
		List<Matrix> lInters;
		FindHyperEdgeInterEll( lA, lB, lABCs, lInters );

		// Check if any
		if( lInters.Count() )
		{
			// Check all
			BrowseList( iI, lInters )
			{
				// Compute error
				double lErr = 0.0;
				{
					const Matrix & lInt = lInters.GetAt( iI );
                    for( size_t d=0 ; d<DIM ; d++ )
						lErr += lInt.Elt(d) * lInt.Elt(d) / (lABCs[d] * lABCs[d]);
					lErr -= 1.0;
				}

				// Log
                LzLogM("", "Error "<<test_count<<"= "<<lErr);

				// Update 
				if( pMaxErr < lErr )
				{
					pMaxErr = lErr;

                    LzLogN("", "New max error");
                    LzLogM("", lABCs.ToString0("ABCs "));
					lA.Log("A");
					lB.Log("B");
				}
			}
			
			// Done
			test_count++;
		}
	}
	while( test_count < 10 );
}

//================================================================================
void ProjectOnHyperEll( const Matrix & pA, const Vector<double> & pABCs, Matrix * ppNear, Matrix * ppFar, Matrix * ppMetric,
                        size_t pMaxIter/*=100*/, double pMaxErr/*=1e-10*/ )
{
	// Determine in which dimension we are working
    const size_t DIM = pABCs.Size();

	// Check
	if( DIM == 0 )
        LzLogException("", "Invalid dimension! Expected at least 1, found "<<pABCs.Size()<<".");

	// Check
	if( pA.Rows()!=DIM || pA.Cols()!=1 )
        LzLogException("", "Invalid hyperpoint A! Expected "<<DIM<<"x1, found "<<pA.Rows()<<"x"<<pA.Cols()<<".");

	// Check
	if( ppNear==nullptr && ppFar==nullptr )
        LzLogException("", "At least one of near and far must be specified! Both are nullptr.");

	//------------------------------------
	// Fork between Metric or No Metric
	//------------------------------------

	if( ppMetric )
	{
		//-------------------
		// Yes metric
		//-------------------

		// Check
		if( ppMetric->Rows()!=DIM || ppMetric->Cols()!=DIM )
            LzLogException("", "Invalid metric! Expected "<<DIM<<"x"<<DIM<<", found "<<ppMetric->Rows()<<"x"<<ppMetric->Cols()<<".");

		
//************** transformer le pb en projection sur ellipsoid dans un espace avec metrique Euclidienne par defaut
//************** transformer le pb en projection sur ellipsoid dans un espace avec metrique Euclidienne par defaut
//************** transformer le pb en projection sur ellipsoid dans un espace avec metrique Euclidienne par defaut
//************** transformer le pb en projection sur ellipsoid dans un espace avec metrique Euclidienne par defaut
//************** transformer le pb en projection sur ellipsoid dans un espace avec metrique Euclidienne par defaut

		
		/*
		// Find Lagrange multiplier Mu
		double lMu = 0.0;
        size_t iIter = 0;
		double lErr;
		do
		{

	//toudou


	lErr = 0.0;


			// Next iteration
			iIter++;
		}
		while( iIter<pMaxIter && lErr>pMaxErr );

		// Final check
		if( lErr > pMaxErr )
            LzLogException("", "Still have error= "<<lErr<<" >= max err= "<<pMaxErr<<", after "<<iIter<<" iteration (s)!");

		//-----------------------
		// Compute solution
		//-----------------------

		Matrix MuQ_plus_R( DIM, DIM );
		MuQ_plus_R.LoadValue( 0.0 );

		// Init diagonal
        for( size_t i=0 ; i<DIM ; i++ )
			MuQ_plus_R.Elt(i, i) = lMu*pABCs[i] + (ppMetric ? ppMetric->Elt(i, i) : 1.0);
	
		// Add metric
		if( ppMetric )
		{
            for( size_t i=0   ; i<DIM ; i++ )
            for( size_t j=i+1 ; j<DIM ; j++ )
				MuQ_plus_R.Elt(i, j) = MuQ_plus_R.Elt(j, i) = ppMetric->Elt(i, j );
		}

		// Inverse
		Matrix MuQ_p_R_inv;
		MuQ_plus_R.LU_InverseTo( MuQ_p_R_inv );

		// Solve
		if( ppMetric )
		{
			Matrix RA;
			ppMetric->Mult( pA, RA );
			MuQ_p_R_inv.Mult( RA, pProj );
		}
		else
			MuQ_p_R_inv.Mult( pA, pProj );
*/
	}
	else
	{
		//-------------------
		// Metric = default
		//-------------------

		// Solution(s)= ( x_i ), outer point= ( o_i ), lambda= Lagrangian, half-axes= { a_i }
		// NB: Solution(s) = near point, and far point
		//
		// x_i( 1 + lambda / a_i^2 ) = o_i
		//
		// Find roots of G(lambda) := sum_i{ a_i^2 o_i^2 / (a_i^2 + lambda)^2 } - 1
		//
		// Derivative: G'(lambda) = sum_i{ -a_i^2 o_i^2 / (a_i^2 + lambda)^3 }
		//

		// Find Lagrange multiplier L
//		double lMu = 0.0;
        size_t iIter = 0;
		double lErr;
		do
		{

	//toudou

			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )
			//renvoyer near et far !!! (avec option a la Eigen : ComputeNear, ComputeFar, ComputeNearAndFar )

	lErr = 0.0;


			// Next iteration
			iIter++;
		}
		while( iIter<pMaxIter && lErr>pMaxErr );

		// Final check
		if( lErr > pMaxErr )
            LzLogException("", "Still have error= "<<lErr<<" >= max err= "<<pMaxErr<<", after "<<iIter<<" iteration (s)!");

		//-----------------------
		// Compute solution
		//-----------------------

//Matrix MuQ_plus_R( DIM, DIM );
//MuQ_plus_R.LoadValue( 0.0 );

//// Init diagonal
//for( size_t i=0 ; i<DIM ; i++ )
//	MuQ_plus_R.Elt(i, i) = lMu*pABCs[i] + (ppMetric ? ppMetric->Elt(i, i) : 1.0);
	
//// Add metric
//if( ppMetric )
//{
//	for( size_t i=0   ; i<DIM ; i++ )
//	for( size_t j=i+1 ; j<DIM ; j++ )
//		MuQ_plus_R.Elt(i, j) = MuQ_plus_R.Elt(j, i) = ppMetric->Elt(i, j );
//}

//// Inverse
//Matrix MuQ_p_R_inv;
//MuQ_plus_R.LU_InverseTo( MuQ_p_R_inv );

//// Solve
//if( ppMetric )
//{
//	Matrix RA;
//	ppMetric->Mult( pA, RA );
//	MuQ_p_R_inv.Mult( RA, pProj );
//}
//else
//	MuQ_p_R_inv.Mult( pA, pProj );


	}
}

//================================================================================
void UnitTest_ProjectOnHyperEll( double & /*pMaxErr*/ )
{
	//***** TODO **************************************************
	//*****
	//***** debug premiere version : calculer point le plus proche en discretisant violemment l'ellipsoide (ca permettra de valider la projection avec metrique,
	//***** avant d'implementer la "bonne" projection avec methode de Newton et criteres pour trouver premiere singularite (near) et deuxieme singularite (far))
	//*****
	//***** TODO **************************************************



	// Case 1: circles (known solution)


	// Case 2: hyper ellipsoids (validation through intensive sampling)
	// Pick a dimension
	// Pick ellipsoid size
	// Pick a bunch of test points
	// Check by sampling the hyperellipsoid


	// Case 3: hyper ellipsoids with non Euclidean metrics

	// same shit


		//toudou

	Matrix A( 3 );
	Vector<double> lABCs = { 12, 34, 56 };
	Matrix P;

	ProjectOnHyperEll( A, lABCs, /*ppNear=*/&P, /*ppFar=*/nullptr, /*ppMetric=*/nullptr );

}

//================================================================================
bool IsZero( double pX, double pTol/*=-1*/ )
{
    if( pTol < 0 )
        return std::abs( pX ) < LzEpsilon;
    else
        return std::abs( pX ) < pTol;

// Better way of doing this? ==> http://floating-point-gui.de/errors/comparison/
}

//================================================================================
//double Deg2Rad( double pDeg )
//{
//    return LzMath::PI * pDeg / 180.0;
//}
//
//================================================================================
//double Rad2Deg( double pRad )
//{
//    return 180.0 * pRad / LzMath::PI;
//}

//================================================================================
double SigmoidPolyRank3( double t )
{
    if( t < 0 )
        return 0;
    else if( t > 1 )
        return 1;
    else
        return t * t * ( 3 - 2 * t );
}

#if defined(USING_QT) && !defined(NO_OPENGL)
//================================================================================
void SetScaleColor( double pMin, double pMax, double pVal )
{
    // Scale colors
    static const double COLS[][3] = { {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0} };
    static const int M = sizeof( COLS ) / sizeof( COLS[0] ) - 1;

    if( pMin < pMax )
    {
        // Normalize pVal
        pVal = ( pVal - pMin ) / ( pMax - pMin );

        // Ensures that pVal in [0,1]
        if( pVal < 0 )
            pVal = 0;
        if( pVal > 1 )
            pVal = 1;

        // Find right interval
        int I_low = ( int )( pVal * M );
        if( I_low == M )
        {
            // Max value color
            glColor4d( COLS[M][0], COLS[M][1], COLS[M][2], 1.0 );
        }
        else
        {
            // Interpolate color
            double a = pVal * M - I_low;
            double b = 1 - a;

            double R1 = COLS[I_low][0];
            double G1 = COLS[I_low][1];
            double B1 = COLS[I_low][2];

            double R2 = COLS[I_low + 1][0];
            double G2 = COLS[I_low + 1][1];
            double B2 = COLS[I_low + 1][2];

            glColor4d( b*R1 + a*R2, b*G1 + a*G2, b*B1 + a*B2, pVal * pVal/**pVal*/ );
        }
    }
    else
    {
        // Default 0 color
        glColor4d( COLS[0][0], COLS[0][1], COLS[0][2], 0.0 );
    }
}

//================================================================================
void SetScaleColor( double pMin, double pMax, double pVal, double pMyAlpha )
{
    // Scale colors
    static const double COLS[][3] = { {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0} };
    static const int M = sizeof( COLS ) / sizeof( COLS[0] ) - 1;

    if( pMin < pMax )
    {
        // Normalize pVal
        pVal = ( pVal - pMin ) / ( pMax - pMin );

        // Ensures that pVal in [0,1]
        if( pVal < 0 )
            pVal = 0;
        if( pVal > 1 )
            pVal = 1;

        // Find right interval
        int I_low = ( int )( pVal * M );
        if( I_low == M )
        {
            // Max value color
            glColor4d( COLS[M][0], COLS[M][1], COLS[M][2], pMyAlpha );
        }
        else
        {
            // Interpolate color
            double a = pVal * M - I_low;
            double b = 1 - a;

            double R1 = COLS[I_low][0];
            double G1 = COLS[I_low][1];
            double B1 = COLS[I_low][2];

            double R2 = COLS[I_low + 1][0];
            double G2 = COLS[I_low + 1][1];
            double B2 = COLS[I_low + 1][2];

            glColor4d( b*R1 + a*R2, b*G1 + a*G2, b*B1 + a*B2, pMyAlpha );
        }
    }
    else
    {
        // Default 0 color
        glColor4d( COLS[0][0], COLS[0][1], COLS[0][2], pMyAlpha );
    }
}
#endif

//================================================================================
//void MeanMaxRMS( const vector<double> & pErr, double & pMean, double & pMax, double & pRMS )
//{
//    pMean = pMax = pRMS = 0;
//
//    if( !pErr.size() )
//        return;
//
//    pMax = -std::numeric_limits<double>::max();
//    for( size_t e = 0 ; e < pErr.size() ; e++ )
//    {
//        double E = pErr[e];
//
//        pMean += E;
//        pRMS += E * E;
//        if( E > pMax ) pMax = E;
//    }
//
//    pMean /= pErr.size();
//    pRMS = std::sqrt( pRMS / pErr.size() );
//}

//================================================================================
//void MeanMaxRMS( const Vector<double> & pErr, double & pMean, double & pMax, double & pRMS )
//{
//    vector<double> lErr( pErr.Size() );
//    for( size_t i = 0 ; i < pErr.Size() ; i++ )
//        lErr[i] = pErr[i];
//
//    MeanMaxRMS( lErr, pMean, pMax, pRMS );
//}

//================================================================================
void MeanMaxStdDev( const vector<double> & pErr,
                    double & pMean, double & pMax, double & pStdDev, size_t * ppMaxIdx/*=nullptr*/ )
{
    // Check
    if( pErr.size() == 0 )
        LzLogException("", "Cannot compute statistics on an empty dataset!")

    // Reset std deviation
    pStdDev = 0;

	// Special case: 1 entry
    if( pErr.size() == 1 )
	{
		pMean = pMax = pErr[0];
        if( ppMaxIdx ) *ppMaxIdx = 0;
		return;
	}

    // Reset mean and max
    pMean = pMax = 0;

    // Mean and max
    pMax = -std::numeric_limits<double>::max();
    for( size_t e=0 ; e<pErr.size() ; e++ )
    {
        // Read value
        double E = pErr[e];

        // Update mean
        pMean += E;

        // Update max
        if( pMax < E )
        {
            pMax = E;

            // Return index of max, if needed
            if( ppMaxIdx ) *ppMaxIdx = e;
        }
    }
    pMean /= pErr.size();

    // Std dev
    for( size_t e=0 ; e<pErr.size() ; e++ )
    {
        double EmeanE = pErr[e] - pMean;
        pStdDev += EmeanE * EmeanE;
    }

    // --> Nope: pStdDev = std::sqrt( pStdDev / (pErr.size() - 1) );
    // --> Warning: sqrt( unbiased est. of variance ) = BIASED estimator of std dev!
    //
    // Approximates the unbiased std dev estimator, see https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
    // The remaining bias is relatively small: for n=3 it is equal to 1.3%, and for n=9 the bias is already 0.1%
    pStdDev = std::sqrt( pStdDev / (pErr.size() - 1.5) );
}

//================================================================================
void MeanMaxStdDev( const Vector<double> & pErr, double & pMean, double & pMax, double & pStdDev, size_t * ppMaxIdx/*=nullptr*/ )
{
    // Convert
    vector<double> lErr( pErr.Size() );
    for( size_t i = 0 ; i < pErr.Size() ; i++ )
        lErr[i] = pErr[i];

    // Compute
    MeanMaxStdDev( lErr, pMean, pMax, pStdDev, ppMaxIdx );
}

//================================================================================
void MeanMaxStdDev( const List<double> & pErr, double & pMean, double & pMax, double & pStdDev, size_t * ppMaxIdx/*=nullptr*/ )
{
    // Convert
    vector<double> lErr;
    pErr.ToVector( lErr );

    // Compute
    MeanMaxStdDev( lErr, pMean, pMax, pStdDev, ppMaxIdx );
}

//================================================================================
Point3D Centroid( const vector<Point3D> & pPoints )
{
    // Check
    if( !pPoints.size() )
        LzLogException("",  "Cannot compute centroid: empty set of points!" );

    // Compute
    Point3D lCentroid( 0, 0, 0 );
    for( const Point3D & iPt : pPoints )
    {
        lCentroid.X() += iPt.X();
        lCentroid.Y() += iPt.Y();
        lCentroid.Z() += iPt.Z();
    }
    lCentroid.X() /= pPoints.size();
    lCentroid.Y() /= pPoints.size();
    lCentroid.Z() /= pPoints.size();

    return lCentroid;
}

//================================================================================
Point3D Centroid( const List<Point3D> & pPoints )
{
    if( !pPoints.Count() )
        LzLogException("",  "Cannot compute centroid: empty set of points!" );

    Point3D lCentroid( 0, 0, 0 );
    BrowseList( iP, pPoints )
    {
		const Point3D lP = pPoints.GetAt( iP );

        lCentroid.X() += lP.X();
        lCentroid.Y() += lP.Y();
        lCentroid.Z() += lP.Z();
    }
    lCentroid.X() /= pPoints.Count();
    lCentroid.Y() /= pPoints.Count();
    lCentroid.Z() /= pPoints.Count();

    return lCentroid;
}

//================================================================================
Point3D Centroid( const Vector<Point3D> & pPoints )
{
//vector<Point3D> lPoints( pPoints.Size() );
//for( size_t i = 0 ; i < pPoints.Size() ; i++ )
//    lPoints[i] = pPoints[i];
//
//return Centroid( lPoints );

    const size_t lSize = pPoints.Size();
	if( !lSize )
        LzLogException("",  "Cannot compute centroid: empty set of points!" );

    Point3D lCentroid( 0, 0, 0 );
    for( size_t i=0 ; i<lSize ; i++ )
    {
		const Point3D & lPt = pPoints[i];
        lCentroid.X() += lPt.X();
        lCentroid.Y() += lPt.Y();
        lCentroid.Z() += lPt.Z();
    }
    lCentroid.X() /= lSize;
    lCentroid.Y() /= lSize;
    lCentroid.Z() /= lSize;

    return lCentroid;
}

//================================================================================
Point3D WeightedCentroid( const Vector<Point3D> & pPoints, const Vector<double> & pWeights )
{
    // Check
    if( pPoints.Size() != pWeights.Size() )
        LzLogException("",  "Cannot compute wighted centroid: vector size mismatch! (" << pPoints.Size() << " != " << pWeights.Size() << ")" );

    Point3D lCentroid;
    double lTotWeight = 0;

    for( size_t i=0 ; i<pPoints.Size() ; i++ )
    {
        const Point3D & lP = pPoints[i];
        const double lW = pWeights[i];

        lTotWeight += lW;

        lCentroid.X() += lW * lP.X();
        lCentroid.Y() += lW * lP.Y();
        lCentroid.Z() += lW * lP.Z();
    }

    // Check
    if( IsZero(lTotWeight) )
        LzLogException("",  "Could not normalize weighted centroid by weight= " << lTotWeight << "!" );

    // Normalize
    lCentroid /= lTotWeight;

    return lCentroid;
}

//================================================================================
//void FindAttractors( const vector<Point3D> & pA, const vector<Point3D> & pB, List< List<size_t> > & pAtoB, const FindAttractorsParams & pFAP,
//					 const Tree3D * ppTree_A/*=nullptr*/, const Tree3D * ppTree_B/*=nullptr*/ )
//{
//    Vector<Point3D> lA((size_t)pA.size() ), lB((size_t)pB.size() );
//
//    for( size_t a = 0 ; a < lA.Size() ; a++ )
//        lA[a] = pA[a];
//
//    for( size_t b = 0 ; b < lB.Size() ; b++ )
//        lB[b] = pB[b];
//
//    FindAttractors( lA, lB, pAtoB, pFAP, ppTree_A, ppTree_B );
//}

//================================================================================
//void FindAttractors( const Vector<Point3D> & pA, const Vector<Point3D> & pB, List< List<size_t> > & pAtoB, const FindAttractorsParams & pFAP,
//					 const Tree3D * ppTree_A/*=nullptr*/, const Tree3D * ppTree_B/*=nullptr*/ )
//{
//#pragma region "Initialize nearest tables"
//    Vector<size_t> lNearest_AtoB( pA.Size() );
//    Vector<size_t> lNearest_BtoA( pB.Size() );
//    {
//		// Detonators: in case we have ownership of pointers
//		LzServices::Detonator lDetTree_A;
//		LzServices::Detonator lDetTree_B;
//
//		Tree3D * lpTree_A;
//		Tree3D * lpTree_B;
//
//		if( pFAP.mUseTree )
//        {
//            //------------------------
//            // Initialize 3D-trees
//            //------------------------
//
//			// A
//			if( ppTree_A )
//			{
//				// Casting away constness here: I swear I won't mess with your tree...
//				lpTree_A = (Tree3D *)ppTree_A;
//			}
//			else
//			{
//				// Create new
//				lpTree_A = new Tree3D;
//				lDetTree_A.Set( [lpTree_A](){ delete lpTree_A; } );
//				//
//				lpTree_A->Create( pA, pFAP.mLogTree, pFAP.mMinDiameter, pFAP.mMaxLeafCount, pFAP.mMaxTreeDepth );
//			}
//
//			// B
//			if( ppTree_B )
//			{
//				// Casting away constness here: I swear I won't mess with your tree...
//				lpTree_B = (Tree3D *)ppTree_B;
//			}
//			else
//			{
//				// Create new
//				lpTree_B = new Tree3D;
//				lDetTree_B.Set( [lpTree_B](){ delete lpTree_B; } );
//				//
//				lpTree_B->Create( pB, pFAP.mLogTree, pFAP.mMinDiameter, pFAP.mMaxLeafCount, pFAP.mMaxTreeDepth );
//			}
//        }
//
//		// Tree aliases
//		const TxGeom::Tree3D & lTree_A = *lpTree_A;
//		const TxGeom::Tree3D & lTree_B = *lpTree_B;
//
//        // A to B
//        for( size_t a = 0 ; a < pA.Size() ; a++ )
//        {
//            const Point3D & lV = pA[a];
//
//            if( pFAP.mUseTree )
//            {
//                //------------------------
//                // 3D-tree
//                //------------------------
//
//                lNearest_AtoB[a] = lTree_B.GetNearestVertex( lV, TxGeom::Tree3D::NearestMode::Accurate );
//            }
//            else
//            {
//                //------------------------
//                // Slow...
//                //------------------------
//
//                double lMinD = +std::numeric_limits<double>::max();
//                for( size_t b = 0 ; b < pB.Size() ; b++ )
//                {
//                    double lD = lV.DistanceTo( pB[b] );
//                    if( lMinD > lD )
//                    {
//                        lMinD = lD;
//                        lNearest_AtoB[a] = b;
//                    }
//                }
//            }
//        }
//
//        // B to A
//        for( size_t b = 0 ; b < pB.Size() ; b++ )
//        {
//            const Point3D & lV = pB[b];
//
//            if( pFAP.mUseTree )
//            {
//                //------------------------
//                // 3D-tree
//                //------------------------
//
//                lNearest_BtoA[b] = lTree_A.GetNearestVertex( lV, TxGeom::Tree3D::NearestMode::Accurate );
//            }
//            else
//            {
//                //------------------------
//                // Slow...
//                //------------------------
//
//                double lMinD = +std::numeric_limits<double>::max();
//                for( size_t a = 0 ; a < pA.Size() ; a++ )
//                {
//                    double lD = lV.DistanceTo( pA[a] );
//                    if( lMinD > lD )
//                    {
//                        lMinD = lD;
//                        lNearest_BtoA[b] = a;
//                    }
//                }
//            }
//        }
//    }
//#pragma endregion
//
//
//    // Initialize loops
//    Vector<bool> lBeenThere_A, lBeenThere_B;
//    lBeenThere_A.Assign( pA.Size(), false );
//    lBeenThere_B.Assign( pB.Size(), false );
//    pAtoB.DelAll();
//
//
//#pragma region "Find all loops starting from A"
//    for( size_t a = 0 ; a < pA.Size() ; a++ )
//    {
//        // Been there?
//        if( lBeenThere_A[a] )
//            continue;
//
//        // Start new loop
//        List<size_t> lLoop;
//        lLoop.AddTail( a );
//        while( true )
//        {
//            // Pop index
//            const size_t lA0 = lLoop.GetTail();
//
//            //..........................
//            // Read mapping A ==> B
//            //..........................
//
//            const size_t lB0 = lNearest_AtoB[lA0];
//            {
//                // Loop already found?
//                if( lBeenThere_B[lB0] )
//                    break;
//
//                // Look for B-index in this loop
//                bool lIsB = false;
//                void * lPosB0 = nullptr;
//                BrowseList( iPos, lLoop )
//                {
//                    if( lIsB && lLoop.GetAt( iPos ) == lB0 )
//                    {
//                        lPosB0 = iPos;
//                        break;
//                    }
//
//                    lIsB = !lIsB;
//                }
//
//                // Found B-index?
//                if( lPosB0 )
//                {
//                    // Rotate list to start with A-index
//                    List<size_t> lFinalLoop;
//                    void * iPos = lLoop.NextPos( lPosB0 );
//                    while( iPos )
//                        lFinalLoop.AddTail( lLoop.GetAtAndNext( iPos ) );
//                    lFinalLoop.AddTail( lB0 );
//
//                    // Stash
//                    pAtoB.AddTail( lFinalLoop );
//                    break;
//                }
//                else
//                {
//                    // Grow list
//                    lLoop.AddTail( lB0 );
//                }
//            }
//
//            //..........................
//            // Read mapping B ==> A
//            //..........................
//
//            const size_t lA1 = lNearest_BtoA[lB0];
//            {
//                // Loop already found?
//                if( lBeenThere_A[lA1] )
//                    break;
//
//                // Look for A-index in this loop
//                bool lIsA = true;
//                void * lPosA1 = nullptr;
//                BrowseList( iPos, lLoop )
//                {
//                    if( lIsA && lLoop.GetAt( iPos ) == lA1 )
//                    {
//                        lPosA1 = iPos;
//                        break;
//                    }
//
//                    lIsA = !lIsA;
//                }
//
//                // Found B-index?
//                if( lPosA1 )
//                {
//                    // Clip list
//                    List<size_t> lFinalLoop;
//                    void * iPos = lPosA1;
//                    while( iPos )
//                        lFinalLoop.AddTail( lLoop.GetAtAndNext( iPos ) );
//
//                    // Stash
//                    pAtoB.AddTail( lFinalLoop );
//                    break;
//                }
//                else
//                {
//                    // Grow list
//                    lLoop.AddTail( lA1 );
//                }
//            }
//        }
//
//        // Mark all vertices as processed
//        bool lIsA = true;
//        BrowseList( iPos, lLoop )
//        {
//            size_t lIdx = lLoop.GetAt( iPos );
//            if( lIsA )
//                lBeenThere_A[lIdx] = true;
//            else
//                lBeenThere_B[lIdx] = true;
//
//            lIsA = !lIsA;
//        }
//    }
//#pragma endregion
//
//
//    // No need to look for loops starting from B!
//    //
//    // Proof:
//    //      Let b be a point in B. If there is a loop containing b then there is an a in A that projects onto b.
//    //      Which means that this loop has already been found as all points in A have been processed.
//    // \blacksquare
//}


#pragma region "Polygons"
//================================================================================
bool PointIsInTriangle( const Point3D & pPt, const Point3D & pA, const Point3D & pB, const Point3D & pC/*, bool pBoundaryIsIn*/ )
{
    // Triangle edges
    const Vector3D lAB = pB - pA;
    const Vector3D lAC = pC - pA;

    // Triangle surface
    const Vector3D lTri = lAB ^ lAC;
    size_t lNonNullCoord = lTri.MaxAbsValIdx(); // Always returns 0, 1 or 2

    // Triangle non null coord
    const double lABAC = lTri.mV[ lNonNullCoord ];

    // Check
    if( IsZero( lABAC ) )
        LzLogException("", "Found flat triangle! Tri= "<<lTri.ToString())

    // Point position
    const Vector3D lAP = pPt - pA;

    // Local coordinate x
    double x = ( lAP ^ lAC ).mV[lNonNullCoord] / ( +lABAC );

    // Local coordinate y
    double y = ( lAP ^ lAB ).mV[lNonNullCoord] / ( -lABAC );

    if( x < 0 || y < 0 || x + y > 1 )
        return false;

    return true;
}

//================================================================================
//  bool PointIsInPolygon( const Point3D & pPt, const vector<Point3D> & pPoly, bool pBoundaryIsIn )
//  {
//      // Find a normal vector to the polygon
//      Point3D pG;
//      Vector3D lEigV[3];
//      PCA( pPoly, pG, lEigV );
//
//      // Compute point position
//      double lAngleSum = 0;
//      for( size_t v=0 ; v<pPoly.size() ; v++ )
//      {
//          // Next index
//          size_t w = (v + 1) % pPoly.size();
//
//          // Two points
//          const Point3D & pV = pPoly[v];
//          const Point3D & lW = pPoly[w];
//
//          // Two vectors
//          Vector3D lVPt = pV - pPt;
//          Vector3D lWPt = lW - pPt;
//
//          // Check if point lies on one of the two vertices v or w
//          if( IsZero(lVPt.Norm()) || IsZero(lWPt.Norm()) )
//          {
//              // Point on boundary
//              return pBoundaryIsIn;
//          }
//
//          // Compute angle
//          double lAngle = lVPt.SignedAngleTo( lWPt, lEigV[2] );
//
//          // Check if point lies on edge [v,w] iff. angle close to +Ï€ or -Ï€, and distance to line is 0
//          if( (std::abs(lAngle-LzMath::PI)<1e-6 || std::abs(lAngle+LzMath::PI)<1e-6) && IsZero( Line3D(pV,lW).DistanceTo(pPt) ) )
//          {
//              // Point on boundary
//              return pBoundaryIsIn;
//          }
//
//          // Point inside or outside
//          lAngleSum += lAngle;
//      }
//
//      // Absolute value
//      double lAbsAngleSum = std::abs( lAngleSum );
//
//      // Out
//      if( lAbsAngleSum < 1e-6 ) // Any tolerance lower than Ï€/2
//          return false;
//      else
//      if( lAbsAngleSum > 2*LzMath::PI - 1e-6 )
//          return true;
//      else
//          return pBoundaryIsIn;
//
//#pragma region "Test PointIsInPolygon"
////vector<Point3D> lPoly;
////lPoly.push_back( Point3D(0,0,0) );
////lPoly.push_back( Point3D(10,0,0) );
////lPoly.push_back( Point3D(5,5,0) );
////lPoly.push_back( Point3D(0,5,0) );
////
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(2,2,0),lPoly,true) ? "IN" : "OUT"));
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(2,2,0),lPoly,false) ? "IN" : "OUT"));
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(2,-2,0),lPoly,true) ? "IN" : "OUT"));
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(2,-2,0),lPoly,false) ? "IN" : "OUT"));
////LzLogM("", "");
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(7.5,2.5,0),lPoly,true) ? "IN" : "OUT"));
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(7.5,2.5,0),lPoly,false) ? "IN" : "OUT"));
////LzLogM("", "");
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(1,0,0),lPoly,true) ? "IN" : "OUT"));
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(1,0,0),lPoly,false) ? "IN" : "OUT"));
////LzLogM("", "");
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(10,0,0),lPoly,true) ? "IN" : "OUT"));
////LzLogM("", "Test= "+(LzMath::ToolBox::PointIsInPolygon(Point3D(10,0,0),lPoly,false) ? "IN" : "OUT"));
//#pragma endregion
//  }

//================================================================================
//bool PointIsInPolygon( const Point3D & pPt, const Vector<Point3D> & pPoly, bool pBoundaryIsIn )
//{
//  vector<Point3D> lPoly( pPoly.Size() );
//  for( size_t p=0 ; p<pPoly.Size() ; p++ )
//      lPoly[p] = pPoly[p];
//
//  return PointIsInPolygon( pPt, lPoly, pBoundaryIsIn );
//}

//================================================================================
PointAndPolygon PointIsInPolygon( const Point3D & pPt, const vector<Point3D> & pPoly )
{
    Vector<Point3D> lPoly( (size_t)pPoly.size() );
    for( size_t p = 0 ; p < pPoly.size() ; p++ )
        lPoly[p] = pPoly[p];

    return PointIsInPolygon( pPt, lPoly );
}

//================================================================================
PointAndPolygon PointIsInPolygon( const Point3D & pPt, const Vector<Point3D> & pPoly )
{
    // Find a normal vector to the polygon
    Point3D pG;
    Vector3D lEigV[3];
    PCA( pPoly, pG, lEigV );

    // Compute point position
    double lAngleSum = 0;
    for( size_t v = 0 ; v < pPoly.Size() ; v++ )
    {
        // Next index
        size_t w = ( v + 1 ) % pPoly.Size();

        // Two points
        const Point3D & lV = pPoly[v];
        const Point3D & lW = pPoly[w];

        // Two vectors
        Vector3D lPtV = lV - pPt;
        Vector3D lPtW = lW - pPt;

        // Check if point lies near one of the two vertices v or w
        if( IsZero( lPtV.Norm(), 1e-6 ) || IsZero( lPtW.Norm(), 1e-6 ) )
        {
            // Point on boundary
            return PointAndPolygon::Surface;
        }

        // Compute angle
        double lAngle = lPtV.SignedAngleTo( lPtW, lEigV[2] );

        // Check if point lies on edge [v,w] iff. angle close to +Pi or -Pi
        if( std::abs(lAngle - LzMath::PI)<1e-6 || std::abs(lAngle + LzMath::PI)<1e-6 )
        {
            // Point on boundary
            return PointAndPolygon::Surface;
        }

        // Point inside or outside
        lAngleSum += lAngle;
    }

    // Absolute value
    double lAbsAngleSum = std::abs( lAngleSum );

    // Out
    if( lAbsAngleSum < 1e-6 ) // Any tolerance lower than Pi/2 will do
        return PointAndPolygon::Out;
    else if( lAbsAngleSum > 2*LzMath::PI - 1e-6 )
        return PointAndPolygon::In;
    else
	{
#if 0
        // Coul Raoul: default to Surface
        return PointAndPolygon::Surface;
#else
        // Maniac mode: check all configurations
        if( fabs(lAbsAngleSum - LzMath::PI) < 1e-6 )
            return PointAndPolygon::Surface;
        else
            LzLogException("",  "Unexpected situation! Point might be out of the polygon plane (angle sum.= " << lAbsAngleSum << ")." );
#endif
	}

#pragma region "Test PointIsInPolygon"
#if 0
    {
        LzLogNode("",  "DEBUG: test LzMath::ToolBox::PointIsInPolygon" );
        using namespace LzMath::ToolBox;

        LzServices::Vector<Point3D> lPoly;
        lPoly.PushBack( Point3D( 0, 0, 0 ) );
        lPoly.PushBack( Point3D( 10, 0, 0 ) );
        lPoly.PushBack( Point3D( 5, 5, 0 ) );
        lPoly.PushBack( Point3D( 0, 5, 0 ) );

#define TOSTR(A) (A==PointAndPolygon::Surface?"Surf":(A==PointAndPolygon::In?"In":"Out"))
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 2, 2, 0 ), lPoly ) ) );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 2, -2, 0 ), lPoly ) ) );
        LzLogM("",  "" );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 7.5, 2.5, 0 ), lPoly ) ) );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 1, 0, 0 ), lPoly ) ) );
        LzLogM("",  "" );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 0, 0, 0 ), lPoly ) ) );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 5, 5, 0 ), lPoly ) ) );
        LzLogM("",  "" );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 10, 10, 0 ), lPoly ) ) );
        LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( -10, 10, 0 ), lPoly ) ) );
        LzLogM("",  "" );
        try
        {
            LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( 2, 2, 5 ), lPoly ) ) );
            LzLogM("",  "Test= " + TOSTR( LzMath::ToolBox::PointIsInPolygon( Point3D( -10, 10, -10 ), lPoly ) ) );
        }
        catch( ... )
        {
            LzLogErr("",  "Exception: point out of polygon!" );
        }
#undef TOSTR
    }
#endif
#pragma endregion
}

//================================================================================
bool SegmentCutsConvexQuad( const Point3D & pA, const Point3D & pB, const vector<Point3D> & pQuad )
{
    // Check
    if( pQuad.size() != 4 )
        LzLogException("",  "Not a quad! (found " << pQuad.size() << " vertices)" );

    // Test if one or both of the extremities are in the polygon, including boundary
    if( PointIsInPolygon( pA, pQuad ) != PointAndPolygon::Out || PointIsInPolygon( pB, pQuad ) != PointAndPolygon::Out )
        return true;

    // Compute utils
    Point3D lQuadCenter = Centroid( pQuad );
    Vector3D lNor = Plane3D( pQuad ).Normal();

    // Both extremities out of quad: test angles
    Vector3D lAO = lQuadCenter - pA;
    Vector3D lBO = lQuadCenter - pB;

    double lMinA = +std::numeric_limits<double>::max();
    double lMinB = +std::numeric_limits<double>::max();
    double lMaxA = -std::numeric_limits<double>::max();
    double lMaxB = -std::numeric_limits<double>::max();

    for( int q = 0 ; q < 4 ; q++ )
    {
        double lAngleA = lAO.SignedAngleTo( pQuad[q] - pA, lNor );
        if( lMinA > lAngleA )
            lMinA = lAngleA;
        if( lMaxA < lAngleA )
            lMaxA = lAngleA;

        double lAngleB = lBO.SignedAngleTo( pQuad[q] - pB, lNor );
        if( lMinB > lAngleB )
            lMinB = lAngleB;
        if( lMaxB < lAngleB )
            lMaxB = lAngleB;
    }

    double lValueA;
    double lValueB;
    try
    {
        lValueA = lAO.SignedAngleTo( pB - pA, lNor );
        lValueB = lBO.SignedAngleTo( pA - pB, lNor );
    }
    catch( ... )
    {
        // Exception thrown if pB-pA is a null vector
        return false;
    }

    return lValueA > lMinA && lValueA < lMaxA && lValueB > lMinB && lValueB < lMaxB;
}

//================================================================================
double PolygonArea( const List<Point3D> & pPoly )
{
    // 1 triangle = 3 points = 4 elements in the list
    if( pPoly.Count() < 4 )
        return 0;

    Vector3D l_Surf( 0, 0, 0 );

    void * i_Pos = pPoly.HeadPos();
    const Point3D l_Root = pPoly.GetAtAndNext( i_Pos );
    Point3D l_P0   = pPoly.GetAtAndNext( i_Pos );
    while( i_Pos )
    {
        Point3D l_P1 = pPoly.GetAtAndNext( i_Pos );

        Vector3D l_V0 = l_P0 - l_Root;
        Vector3D l_V1 = l_P1 - l_Root;

        l_Surf += l_V0 ^ l_V1;

        l_P0 = l_P1;
    }

    return 0.5 * l_Surf.Norm();
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//bool SplitPolygon( const List<Point3D> & pPoly, int pLink1, int pLink2, List<Point3D> & pP1, List<Point3D> & pP2 )
//{
//  // 1 triangle = 3 points = 4 elements in the list
//  if( pPoly.Count() < 4 )
//      return false;
//
//  // Lien degenere
//  if( pLink1 == pLink2 )
//      return false;

//  // Sort links
//  int lLink1 = pLink1<pLink2 ? pLink1 : pLink2 ;
//  int lLink2 = pLink1>pLink2 ? pLink1 : pLink2 ;

//  // Vertices count (last == first)
//  const int lNodeCount = pPoly.Count() - 1;

//  // Impossible link
//  if( lLink1<0 || lLink2>=lNodeCount )
//      return false;

////    // Test si le lien suggere existe deja : les 2 liens se suivent
////    if( l_Link1==l_Link2-1 || l_Link1==(l_Link2+1)%l_NodeCount )
////        return FALSE;

////    // Verifier que le lien ne croise pas le polygone
////    CPoint3D l_A = p_Poly.GetAt( p_Poly.FindIndex(l_Link1) );
////    CPoint3D l_B = p_Poly.GetAt( p_Poly.FindIndex(l_Link2) );
////    CPoint3D l_Tmp;
////    POSITION i_Pos = p_Poly.GetHeadPosition() ;
////    CPoint3D l_I = p_Poly.GetNext(i_Pos);
////    while( i_Pos )
////    {
////        CPoint3D l_J = p_Poly.GetNext(i_Pos);

////        ECut l_Cut = SegCutSeg(l_A,l_B,l_I,l_J,l_Tmp);
////        if( l_Cut==SCS_REAL_CUT || l_Cut==SCS_CUT_I || l_Cut==SCS_CUT_J )
////            return FALSE;

////        l_I = l_J;
////    }

////    // Construction des polygones
////    p_P1.RemoveAll();
////    p_P2.RemoveAll();

////    i_Pos = p_Poly.GetHeadPosition();
////    for( int i=0 ; i<l_NodeCount ; i++ )
////    {
////        CPoint3D l_Pt = p_Poly.GetNext(i_Pos);

////        // P1
////        if( i<=l_Link1 || i>=l_Link2 )
////            p_P1.AddTail(l_Pt);

////        // P2
////        if( i>=l_Link1 && i<=l_Link2 )
////            p_P2.AddTail(l_Pt);
////    }

////    // Fermeture des polygones
////    p_P1.AddTail( p_P1.GetHead() );
////    p_P2.AddTail( p_P2.GetHead() );

//  // Test des surfaces
//  double lSTot = PolygonArea(pPoly);
//  double lS1   = PolygonArea(pP1);
//  double lS2   = PolygonArea(pP2);

//  return lS1<lSTot && lS2<lSTot;
//}

//================================================================================
//  bool MeshPolygon( const List<Point3D> & pPoly, List<Point3D> & pMesh )
//  {
////---------------------------------------------------------------------------------
///*! Triangule le polygone
// * Mesh = (T0, T1, ... ) = (P00,P01,P02 , P10,P12,P13 , ... )
// *///------------------------------------------------------------------------------
//      // 1 triangle = 3 points = 4 elements in the list
//      if( pPoly.Count() < 4 )
//          return false;
//
//      // 1 triangle
//      if( pPoly.Count() == 4 )
//      {
//  //      POSITION i_Pos = p_Poly.GetHeadPosition();
//  //      p_Mesh.AddTail( p_Poly.GetNext(i_Pos) );
//  //      p_Mesh.AddTail( p_Poly.GetNext(i_Pos) );
//  //      p_Mesh.AddTail( p_Poly.GetNext(i_Pos) );
//          return true;
//      }
//      else
//      {
//  //      // Recherche de split
//  //      CList<CPoint3D,CPoint3D&> l_Split1;
//  //      CList<CPoint3D,CPoint3D&> l_Split2;
//
//  //      double l_ST = PolygonSurface(p_Poly);
//  //      double l_BestScore;
//  //      int l_BestLink1 = -1;
//  //      int l_BestLink2;
//
//  //      for( int i_Link1=0 ; i_Link1<p_Poly.GetCount() ; i_Link1++ )
//  //      for( int i_Link2=0 ; i_Link2<p_Poly.GetCount() ; i_Link2++ )
//  //      {
//  //          if( SplitPolygon(p_Poly,i_Link1,i_Link2,l_Split1,l_Split2) )
//  //          {
//  //              double l_S1 = PolygonSurface(l_Split1);
//  //              double l_S2 = PolygonSurface(l_Split2);
//  //              double l_Score = (l_S1-l_ST/2)*(l_S2-l_ST/2);
//  //              if( l_BestLink1<0 || l_Score>l_BestScore )
//  //              {
//  //                  l_BestScore = l_Score;
//  //                  l_BestLink1 = i_Link1;
//  //                  l_BestLink2 = i_Link2;
//
//  //                  if( l_BestScore == 0 )
//  //                      break;
//  //              }
//  //          }
//  //      }
//
//  //      // Found a split ?
//  //      if( l_BestLink1 == -1 )
//  //      {
//              return false;
//  //      }
//  //      else
//  //      {
//  //          // Perform best split
//  //          SplitPolygon(p_Poly,l_BestLink1,l_BestLink2,l_Split1,l_Split2);
//  //          return MeshPolygon(l_Split1,p_Mesh) && MeshPolygon(l_Split2,p_Mesh);
//  //      }
//      }
//  }
#pragma endregion


//================================================================================
unsigned int RotL_32bit( unsigned int pV, int pShift )
{
    int lS = pShift >= 0 ? pShift % 32 : -( ( -pShift ) % 32 ) ;
    return ( pV << lS ) | ( pV >> ( 32 - lS ) );
}

//================================================================================
unsigned int RotR_32bit( unsigned int pV, int pShift )
{
    int lS = pShift >= 0 ? pShift % 32 : -( ( -pShift ) % 32 ) ;
    return ( pV >> lS ) | ( pV << ( 32 - lS ) );
}

#pragma endregion


#pragma region "Piecewise linear function"
//=============================================================================
PiecewiseLinearFunction::PiecewiseLinearFunction()
{
    vector<double> lXs( 2 ), lYs( 2 );
    lXs[0] = lYs[0] = 0;
    lXs[1] = lYs[1] = 1;

    SetSamples( lXs, lYs );
}

//=============================================================================
PiecewiseLinearFunction::PiecewiseLinearFunction( const Vector<double> & pXs, const Vector<double> & pYs )
{
    SetSamples( pXs, pYs );
}

//=============================================================================
PiecewiseLinearFunction::~PiecewiseLinearFunction()
{
    Free();
}

//=============================================================================
void PiecewiseLinearFunction::Free()
{
    mXs.Free();
    mYs.Free();
    mAs.Free();
}

//=============================================================================
void PiecewiseLinearFunction::SetSamples( const Vector<double> & pXs, const Vector<double> & pYs )
{
    // Check
    if( pXs.Size() != pYs.Size() || pXs.Size() < 2 )
        LzLogException("",  "Unable to define function! Invalid sample points count: " << pXs.Size() << " Xs and " << pYs.Size() << " Ys." );

#if 1
    // Check that the data is sorted
    for( size_t x=0 ; x+1<pXs.Size() ; x++ )
    {
        // Check sample spacing + margin
        static const double sMarginX = 1e-5;
        if( pXs[x+1] < pXs[x] + sMarginX  )
        {
            LzLogException("", "Error in X samples: "+std::to_string(pXs[x+1])+" < ("+std::to_string(pXs[x+1])+
                               " + "+std::to_string(sMarginX)+") !")
        }
    }
#else
    // Sort
    Vector<double> lXs = pXs;
    Vector<double> lYs = pYs;
    QuickSort( lXs.Buffer(), lXs.Size(), lYs.Buffer() );

    // Check for redundant points
    for( size_t i = 0 ; i < lXs.Size() - 1 ; i++ )
    {
        if( IsZero( lXs[i + 1] - lXs[i] ) )
            LzLogException("",  "Unable to define function! Redundant sample point: X= " << lXs[i] << "." );
    }
#endif

    // Free previous
    Free();

    // Commit
    mXs = pXs;
    mYs = pYs;

    // Compute slopes
    mAs.Resize( mXs.Size() - 1 );
    for( size_t i = 0 ; i < mAs.Size() ; i++ )
        mAs[i] = ( mYs[i + 1] - mYs[i] ) / ( mXs[i + 1] - mXs[i] );

//// Log
//LzLogMsg("", "Set "+lXs.Size()+" sample points in piecewise linear function.");
}

//=============================================================================
double PiecewiseLinearFunction::EvalAt( double pX, bool pClampValue/*=false*/ ) const
{
    // No need for prior Check a PLF is always valid (default function: F(x) = x)

    // Check if clamping needed
    if( pClampValue )
    {
        if( pX <= mXs[0] )
            return mYs[0];

        if( pX >= mXs.GetBack() )
            return mYs.GetBack();
    }

    // Find interval
    size_t i=1;
    for( ; i<mXs.Size()-1 ; i++ )
    {
        if( pX <= mXs[i] )
            break;
    }

    // Interval [i-1, i]
    return mYs[i - 1] + mAs[i - 1] * ( pX - mXs[i - 1] );
}

//=============================================================================
void PiecewiseLinearFunction::Scale( double pScale )
{
    // Copy Xs
    Vector<double> lXs = mXs;

    // Apply scale factor to Ys
    Vector<double> lYs( mYs.Size() );
    for( size_t s=0 ; s<lYs.Size() ; s++ )
        lYs[s] = mYs[s] * pScale;

    // Commit
    SetSamples( lXs, lYs );
}

//=============================================================================
void PiecewiseLinearFunction::Log() const
{
//    LzLogMsg("",  "Piecewise linear function has " << mXs.Size() << " sample points." );

    LzLogN("", "Piecewise linear function")
    for( size_t i=0 ; i<mXs.Size() ; i++ )
      LzLogM("", "F( "<<mXs[i]<<" )= "<<mYs[i])
}

//=============================================================================
std::string PiecewiseLinearFunction::ToString() const
{
	std::stringstream lStrStr;
	lStrStr << "{ ";
    for( size_t i=0 ; i<mXs.Size() ; i++ )
        lStrStr << mXs[i] << " " << mYs[i] << " ";

    lStrStr << "}";

    return lStrStr.str();

//String ^lStr = "{ ";
//for( size_t i = 0 ; i < mXs.Size() ; i++ )
//    lStr += "" + mXs[i] + " " + mYs[i] + " ";
//
//return lStr + "}";
}

//=============================================================================
void PiecewiseLinearFunction::FromString( const std::string & /*pFrom_Std*/ )
{
#if 0
	String ^pFrom = gcnew String( pFrom_Std.c_str() );

    cli::array<wchar_t> ^lSpaceTabNewLine = { ' ', '\t', '\n', '\r\n' };
    cli::array<String^> ^lItems = pFrom->Split( lSpaceTabNewLine, System::StringSplitOptions::RemoveEmptyEntries );

    // Check length
    if( lItems->Length < 4 || lItems->Length % 2 )
        LzLogException("",  "Cannot parse function! Found " << lItems->Length << " items in string." );

    // Read
    Vector<double> lXs( lItems->Length / 2 );
    Vector<double> lYs( lItems->Length / 2 );
    for( size_t i = 0 ; i < lXs.Size() ; i++ )
    {
        lXs[i] = LzServices::StrToDouble( lItems[2 * i + 0] );
        lYs[i] = LzServices::StrToDouble( lItems[2 * i + 1] );
    }

    // Commit
    SetSamples( lXs, lYs );
#else
    LzLogException("", "*** TODO PiecewiseLinearFunction::FromString( const std::string & pFrom )")
#endif
}
#pragma endregion


#pragma region "Eigen values"
#pragma region "Cubic roots"
//================================================================================
//void CubicPolyRoots( Vector<Complex> & pRoots, double a, double b, double c, double d )
//{
//    ////////////////////////////////////////////////     Implementer la version NRC !!!!
//    ////////////////////////////////////////////////     Implementer la version NRC !!!!
//    ////////////////////////////////////////////////     Implementer la version NRC !!!!
//    ////////////////////////////////////////////////     Implementer la version NRC !!!!
//    ////////////////////////////////////////////////     Implementer la version NRC !!!!
//    ////////////////////////////////////////////////     Implementer la version NRC !!!!
//
//    // Source: http://en.wikipedia.org/wiki/Cubic_function
//
//    // Check
//    if( IsZero( a ) )
//        LzLogException("",  "Found null 3rd degree monomial! (" << a << " x Z^3)" );
//
//    // Resize
//    pRoots.Resize( 3 );
//
//    //----------------------------------------------------------------------------------------------
//    //
//    // Normalization ==> catastrophic results.....
//    //
//    ///-- Eigen values for E
//    //| ( 211.112628241819 ) + i( 0 )
//    //| ( 211.112628241819 ) + i( 0 )
//    //| ( 41.1385329814317 ) + i( 0 )
//    //|
//    //| Err= 443613.09774871
//    //|
//    //| a= 1, b= -463.36378946507, c= 59002.1471223113, d= -1269083.56424396
//    //\--
//    //
//    //
//    ///-- Eigen Jacobi
//    //| 237.02568010392
//    //| 199.5
//    //| 26.8381093611496
//    //\--
//    //
//    //
//    // For large a, b, c and/or d, the normalization reduces by a tiny bit the residual error.
//    //
//    // This code is left here only because it doesn't take up too much time... but could
//    // as well be removed.
//    //
//    // Normalize parameters
//    ////////////////////////double lMaxAbs = 0;
//    ////////////////////////if( lMaxAbs < abs(a) ) lMaxAbs = abs(a);
//    ////////////////////////if( lMaxAbs < abs(b) ) lMaxAbs = abs(b);
//    ////////////////////////if( lMaxAbs < abs(c) ) lMaxAbs = abs(c);
//    ////////////////////////if( lMaxAbs < abs(d) ) lMaxAbs = abs(d);
//    ////////////////////////
//    ////////////////////////// Normalize
//    ////////////////////////a /= lMaxAbs;
//    ////////////////////////b /= lMaxAbs;
//    ////////////////////////c /= lMaxAbs;
//    ////////////////////////d /= lMaxAbs;
//    //
//    //----------------------------------------------------------------------------------------------
//
//    // Discriminants
//    const double D = 18 * a * b * c * d - 4 * b * b * b * d + b * b * c * c - 4 * a * c * c * c - 27 * a * a * d * d;
//    const double D0 = b * b - 3 * a * c;
//
////static const double sTol = 1e-12;
//    static const double sTol = 1e-20;
//
//    // Specific cases
//    if( !IsZero( D, sTol ) && IsZero( D0, sTol ) )
//    {
//        const double D1 = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d;
//        Complex C = Complex( D1, 0 ).Pow( 1.0 / 3.0 );
//
//        // Roots of unit
//        Complex u[3];
//        u[0] = Complex( 1, 0 );
//        u[1] = Complex( -0.5, +sqrt( 3 ) / 2 );
//        u[2] = Complex( -0.5, -sqrt( 3 ) / 2 );
//
//        // Roots
//        for( int k = 0 ; k < 3 ; k++ )
//            pRoots[k] = ( -1 / ( 3 * a ) ) * ( b + u[k] * C + D0 / ( u[k] * C ) );
//    }
//    else if( IsZero( D, sTol ) && IsZero( D0, sTol ) )
//    {
//        // Roots
//        pRoots[0] = pRoots[1] = pRoots[2] = -b / ( 3 * a );
//    }
//    else if( IsZero( D, sTol ) && !IsZero( D0, sTol ) )
//    {
//        // Roots
//        pRoots[0] = pRoots[1] = ( 9 * a * d - b * c ) / ( 2 * D0 );
//        pRoots[2] = ( 4 * a * b * c - 9 * a * a * d - b * b * b ) / ( a * D0 );
//    }
//    else
//    {
//        const double D1 = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d;
//        const double D12_m_4D03 = -27 * a * a * D;
//        Complex C;
//        {
//            Complex R1 = D1 + Complex::Sqrt( D12_m_4D03 );
//            Complex R2 = D1 - Complex::Sqrt( D12_m_4D03 );
//            const Complex & R = R1.Modulus() > R2.Modulus() ? R1 : R2 ;
//
//            C = ( 0.5 * R ).Pow( 1.0 / 3.0 );
//        }
//
//        // Roots of unit
//        Complex u[3];
//        u[0] = Complex( 1, 0 );
//        u[1] = Complex( -0.5, +sqrt( 3 ) / 2 );
//        u[2] = Complex( -0.5, -sqrt( 3 ) / 2 );
//
//        // Roots
//        for( int k = 0 ; k < 3 ; k++ )
//            pRoots[k] = ( -1 / ( 3 * a ) ) * ( b + u[k] * C + D0 / ( u[k] * C ) );
//    }
//}

//================================================================================
//double CheckCubicPolyRoots( const Vector<Complex> & pRoots, double a, double b, double c, double d )
//{
//    // Check
//    if( pRoots.Size() != 3 )
//        LzLogException("",  "Needs exactly 3 cubic roots to check! (" << pRoots.Size() << " != 3)" );
//
//    // Check
//    double lMaxModulus = 0;
//    for( int k = 0 ; k < 3 ; k++ )
//    {
//        // The complex
//        const Complex & xk = pRoots[k];
//
//        // Check
//        Complex lResidu = a * xk * xk * xk + b * xk * xk + c * xk + d;
//
//        double lMod = lResidu.Modulus();
//        if( lMaxModulus < lMod )
//            lMaxModulus = lMod;
//    }
//
//    return lMaxModulus;
//}

//================================================================================
double AccurateRealRoot( double pInitRoot, double pMaxErr, size_t pMaxIter /*** UTILE ***/, double a, double b, double c, double d )
{
#if 1
    const double lInitResid = a * pInitRoot * pInitRoot * pInitRoot + b * pInitRoot * pInitRoot + c * pInitRoot + d;
    const double lInitSlope = 3 * a * pInitRoot * pInitRoot + 2 * b * pInitRoot + c;

    // Find positive/negative approximations
    double lRoot_Pos, lRoot_Neg;


    if( lInitResid > +pMaxErr )
    {
//Loop_PosRoot:
        lRoot_Pos = pInitRoot;

        lRoot_Neg = lRoot_Pos - 2 * lInitResid / lInitSlope;

        double lNewResid = a * lRoot_Neg * lRoot_Neg * lRoot_Neg + b * lRoot_Neg * lRoot_Neg + c * lRoot_Neg + d;
        LzLogM("",  "init resid POS" );
        if( lNewResid > 0 )
        {
            LzLogM("",  "new resid POS" );
//***** CHECK
//***** CHECK
//***** CHECK
            return AccurateRealRoot( lRoot_Neg, pMaxErr, pMaxIter, a, b, c, d );
//***** CHECK
//***** CHECK
//***** CHECK
        }
    }
    else if( lInitResid < -pMaxErr )
    {
        lRoot_Neg = pInitRoot;

        lRoot_Pos = lRoot_Neg - 2 * lInitResid / lInitSlope;

        double lNewResid = a * lRoot_Pos * lRoot_Pos * lRoot_Pos + b * lRoot_Pos * lRoot_Pos + c * lRoot_Pos + d;
        LzLogM("",  "init resid NEG" );

        if( lNewResid < 0 )
        {
            LzLogM("",  "new resid NEG" );
//***** CHECK
//***** CHECK
//***** CHECK
            return AccurateRealRoot( lRoot_Pos, pMaxErr, pMaxIter, a, b, c, d );
//***** CHECK
//***** CHECK
//***** CHECK
        }
    }
    else
        return pInitRoot;

    {
        LzLogN("",  "" );
        LzLogM("",  "Neg root= " << lRoot_Neg << ", residu= " << ( a * lRoot_Neg * lRoot_Neg * lRoot_Neg + b * lRoot_Neg * lRoot_Neg + c * lRoot_Neg + d ) << "." );
        LzLogM("",  "Pos root= " << lRoot_Pos << ", residu= " << ( a * lRoot_Pos * lRoot_Pos * lRoot_Pos + b * lRoot_Pos * lRoot_Pos + c * lRoot_Pos + d ) << "." );


    }

    // Dichotomy
//    const double lInitDelta = abs( lRoot_Pos - lRoot_Neg );

    double lNew_Root = 0.5 * ( lRoot_Pos + lRoot_Neg ); ///********************** boucle while pas jolie .... a corriger......
    while( abs( lRoot_Pos - lRoot_Neg ) > pMaxErr )
    {
//LzLogN("", "");
//LzLogM("", "Pos root= "+lRoot_Pos );
//LzLogM("", "Neg root= "+lRoot_Neg );
//LzLogM("", "");
//LzLogM("", "delta= "+abs(lRoot_Pos-lRoot_Neg)+", max err= "+pMaxErr );
//LzLogM("", "delta ratio= "+(lInitDelta/abs(lRoot_Pos-lRoot_Neg)) );

        lNew_Root = 0.5 * ( lRoot_Pos + lRoot_Neg );
        double lNewResid = a * lNew_Root * lNew_Root * lNew_Root + b * lNew_Root * lNew_Root + c * lNew_Root + d;

//LzLogM("", "");
//LzLogM("", "New root= "+lNew_Root );
//LzLogM("", "New resid= "+lNewResid);


        if( lNewResid > 0 )
            lRoot_Pos = lNew_Root;
        else if( lNewResid < 0 )
            lRoot_Neg = lNew_Root;
        else
            return lNew_Root;

    }

    return lNew_Root;
#else
    double lRoot = pInitRoot;
    size_t lIter = 0;
    double lErr = abs( a * lRoot * lRoot * lRoot + b * lRoot * lRoot + c * lRoot + d );

    // Dichotomy
    double lPosRoot;
    double lNegRoot;

    while( lErr > pMaxErr )
    {
        // Compute slope
        double p = 3 * a * lRoot * lRoot + 2 * b * lRoot + c;

        // Compute residual
        double r = a * lRoot * lRoot * lRoot + b * lRoot * lRoot + c * lRoot + d;

        // Displace root
        double lNewRoot = lRoot - r / p;

        lRoot = lNewRoot;

        // Update error
        lErr = abs( a * lRoot * lRoot * lRoot + b * lRoot * lRoot + c * lRoot + d );


//LzLogM("", "r= "+r+", p= "+p+", Err= "+lErr);


        // Next iteration
        lIter++;

        // Check
        if( lIter > pMaxIter )
//LzLogException("", "Maximal number of iterations ("+pMaxIter+") reached! Err= "+lErr+".");
            break;
    }

    return lRoot;
#endif
}

//================================================================================
//double AccurateRealRoot_Dichot( double pInitRoot, double pMaxErr, size_t pMaxIter, double a, double b, double c, double d )
//{


//}
#pragma endregion


//================================================================================
//Given the eigenvalues d[0..n-1] and (optionally) the eigenvectors v[0..n-1][0..n-1] as determined
//by Jacobi (11.1) or tqli (11.4), this routine sorts the eigenvalues into descending
//order and rearranges the columns of v correspondingly. The method is straight insertion.
//static void EigenSort( /*VecDoub_IO &d*/ double *&d, int pSize_d, /*MatDoub_IO*/ double ***v=0 )
void EigenSort( Matrix & d, /*int pSize_d,*/ /*MatDoub_IO*/ Matrix * v/*=nullptr*/ )
{
    int k;
    int n = d.Rows();
    for( int i = 0 ; i < n - 1 ; i++ )
    {
        double p = d.Elt( k = i, 0 );
        for( int j = i ; j < n ; j++ )
            if ( d.Elt( j, 0 ) >= p )
                p = d.Elt( k = j, 0 );

        if( k != i )
        {
            d.Elt( k, 0 ) = d.Elt( i, 0 );
            d.Elt( i, 0 ) = p;
            if( v )
            {
                for( int j = 0 ; j < n ; j++ )
                {
                    p = v->Elt( j, i );         //p=(*v)[j][i]; careful with NRC index order!
                    v->Elt( j, i ) = v->Elt( j, k ); //(*v)[j][i]=(*v)[j][k];
                    v->Elt( j, k ) = p;         //(*v)[j][k]=p;
                }
            }
        }
    }
}

//================================================================================
//static void rot( /*MatDoub_IO &*/double **&a, double/*Doub*/ s, double/*Doub*/ tau,
//               int/*Int*/ i, int/*Int*/ j, int/*Int*/ k, int/*Int*/ l )
static void rot( Matrix & a, double s, double tau, int i, int j, int k, int l )
{
    double g = a.Elt( i, j ); //a[i][j];
    double h = a.Elt( k, l ); //a[k][l];
    a.Elt( i, j ) = g - s * ( h + g * tau ); //a[i][j] = g - s*(h+g*tau);
    a.Elt( k, l ) = h + s * ( g - h * tau ); //a[k][l] = h + s*(g-h*tau);
}

//================================================================================
void EigenJacobi( const Matrix & pA, Matrix & pDiag, Matrix & pVecs )
{
    if( !pA.Rows() || pA.Rows() != pA.Cols() )
        LzLogException("",  "Non-square matrix!" );

    int n = pA.Rows();
    int nrot = 0;
    double EPS = 1e-12;

    Matrix a( n, n ), v( n, n ), d( n, 1 );

    double *b = new double[n];
    double *z = new double[n];

    int i, j, ip, iq;
    double tresh, theta, tau, t, sm, s, h, g, c;

    // Initialize a
    a = pA;

    // Initialize v to the identity matrix.
    v.LoadIdentity();

    // Initialize b and d to the diagonal of a.
    for( ip = 0 ; ip < n ; ip++ )
    {
        b[ip] = d.Elt( ip, 0 ) = a.Elt( ip, ip ); //b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0; //z[ip] = 0.0; //This vector will accumulate terms of the form tapq as in equation (11.1.14).
    }

    for( i = 1 ; i <= 50 ; i++ )
    {
        sm = 0.0;

        //Sum magnitude of off-diagonal elements.
        for( ip = 0 ; ip < n - 1 ; ip++ )
        {
            for( iq = ip + 1 ; iq < n ; iq++ )
                sm += fabs( a.Elt( ip, iq ) ); //fabs(a[ip][iq]);
        }

        // The normal return, which relies on quadratic convergence to machine underflow.
        if( IsZero( sm, EPS ) ) // == 0.0)
        {
            EigenSort( d, &v );
            delete [] b;
            delete [] z;

            pDiag = d; //***** ou renvoyer matrice diagonale ??
            pVecs = v;

            // Finished ok
            return;
        }

        if( i < 4 )
            // On the first three sweeps...
            tresh = 0.2 * sm / ( n * n );
        else
            // ...thereafter.
            tresh = 0.0;

        for( ip = 0 ; ip < n - 1 ; ip++ )
        {
            for( iq = ip + 1 ; iq < n ; iq++ )
            {
                g = 100.0 * fabs( a.Elt( ip, iq ) ); //  [ip][iq]);

                //After four sweeps, skip the rotation if the off-diagonal element is small.
                if( i > 4 && g <= EPS * fabs( d.Elt( ip, 0 )/*[ip]*/ ) && g <= EPS * fabs( d.Elt( iq, 0 )/*[iq]*/ ) )
                    a.Elt( ip, iq ) = 0; //a[ip][iq]=0.0;
                else if( fabs( a.Elt( ip, iq )/*[ip][iq]*/ ) > tresh )
                {
                    h = d.Elt( iq, 0 )/*[iq]*/ - d.Elt( ip, 0 )/*[ip]*/;
                    if ( g <= EPS * fabs( h ) )
                        t = ( a.Elt( ip, iq )/*[ip][iq]*/ ) / h; // t = 1/(2*Theta)
                    else
                    {
                        theta = 0.5 * h / ( a.Elt( ip, iq )/*[ip][iq]*/ ); // Equation (11.1.10).
                        t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta ) );
                        if ( theta < 0.0 )
                            t = -t;
                    }
                    c = 1.0 / sqrt( 1 + t * t );
                    s = t * c;
                    tau = s / ( 1.0 + c );
                    h = t * a.Elt( ip, iq ); //[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d.Elt( ip, 0 )/*[ip]*/ -= h;
                    d.Elt( iq, 0 )/*[iq]*/ += h;
                    a.Elt( ip, iq )/*[ip][iq]*/ = 0.0;

                    // Case of rotations 0 <= j < p.
                    for ( j = 0; j < ip; j++ )
                        rot( a, s, tau, j, ip, j, iq );

                    // Case of rotations p < j <q.
                    for ( j = ip + 1; j < iq; j++ )
                        rot( a, s, tau, ip, j, j, iq );

                    // Case of rotations q < j <n.
                    for ( j = iq + 1; j < n; j++ )
                        rot( a, s, tau, ip, j, iq, j );

                    for ( j = 0; j < n; j++ )
                        rot( v, s, tau, j, ip, j, iq );

                    ++nrot;
                }
            }
        }

        for ( ip = 0; ip < n; ip++ )
        {
            b[ip] += z[ip];
            d.Elt( ip, 0 ) = b[ip]; // Update d with the sum of tapq,
            z[ip] = 0.0;    // and reinitialize z.
        }
    }

    // Failed
    delete [] b;
    delete [] z;
    LzLogException("",  "Too many iterations!" );
}

#if 0
{
    // CHECK CODE
    Matrix A( 2, 2 );
    A.Elt( 0, 0 ) = +5.0;
    A.Elt( 0, 1 ) = -2.0;
    A.Elt( 1, 0 ) = -2.0;
    A.Elt( 1, 1 ) = +2.0;

    Matrix D, V;
    LzMath::ToolBox::EigenJacobi( A, D, V );

    D.Log();
    V.Log();

    //Matrix:
    //  +5.000000e+00  -2.000000e+00
    //  -2.000000e+00  +2.000000e+00
    //
    //Eigenvalues:
    //  Number 1: +6.000000e+00
    //  Number 2: +1.000000e+00
    //
    //Eigenvectors:
    //  Number 1: +8.944272e-01 -4.472136e-01
    //  Number 2: +4.472136e-01 +8.944272e-01
    //
    //Number of Jacobi rotations = 1
}
#endif

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
void PolarDecompose_OrthSym( const Matrix & pA, Matrix & pOrth, Matrix & pSym )
{
	// Check
	if( !pA.Rows() || pA.Rows()!=pA.Cols() )
        LzLogException("", "Matrix A is empty or not square!");

	// Assemble sym def pos mat
	Matrix At, AtA;
	pA.TransposeTo( At );
	At.Mult( pA, AtA );

	// AtA = V D Vt
	Matrix AtA_D, AtA_V, AtA_Vt;
	EigenJacobi( AtA, AtA_D, AtA_V );
	AtA_V.TransposeTo( AtA_Vt );

	// Square root of AtA_D
	Matrix AtA_SqrtD( AtA_V.Rows(), AtA_V.Cols() );
	AtA_SqrtD.LoadValue( 0 );
    for( size_t k=0 ; k<AtA_SqrtD.Rows() ; k++ )
		AtA_SqrtD.Elt(k,k) = sqrt( AtA_D.Elt(k,0) );

	// S = square root of AtA
	{
		Matrix Sq_Vt;
		AtA_SqrtD.Mult( AtA_Vt, Sq_Vt );
		AtA_V.Mult( Sq_Vt, pSym );
	}

	// Inverse S
	Matrix invS;
	pSym.LU_InverseTo( invS );

	// O = A invS
	pA.Mult( invS, pOrth );
}

//================================================================================
void PolarDecompose_SymOrth( const Matrix & pA, Matrix & pSym, Matrix & pOrth )
{
	// At = Q S ==> A = St Qt = S Qt = S O
	Matrix At, Q;
	pA.TransposeTo( At );
	PolarDecompose_OrthSym( At, Q, pSym );

	// O = Qt
	Q.TransposeTo( pOrth );
}

//================================================================================
void PCA( const List<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] )
{
    vector<Point3D> lCloud;
	pCloud.ToVector( lCloud );

    PCA( lCloud, pMean, pV );
}

//================================================================================
void PCA( const Vector<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] )
{
    vector<Point3D> lCloud( pCloud.Size() );
    for( size_t i = 0 ; i < pCloud.Size() ; i++ )
        lCloud[i] = pCloud[i];

    PCA( lCloud, pMean, pV );
}

//================================================================================
void PCA( const vector<Point3D> & pCloud, Point3D & pMean, Vector3D pV[3] )
{
    const int lNbPoints = (int)pCloud.size();

    // Check
    if( !lNbPoints )
        LzLogException("", "Cannot compute PCA on empty points cloud!")

    // Compute mean point
    pMean.Reset();
    vector<Point3D>::const_iterator iPt;
    for( iPt = pCloud.begin() ; iPt != pCloud.end() ; iPt++ )
    {
        pMean.X() += iPt->X();
        pMean.Y() += iPt->Y();
        pMean.Z() += iPt->Z();
    }

    // Normalize mean
    pMean.X() /= lNbPoints;
    pMean.Y() /= lNbPoints;
    pMean.Z() /= lNbPoints;

    // Fill 0-centered matrix
    Matrix lM( 3, lNbPoints );
    {
        size_t i = 0;

        vector<Point3D>::const_iterator iPt;
        for( iPt = pCloud.begin() ; iPt != pCloud.end() ; iPt++ )
        {
            lM.Elt( 0, i ) = iPt->X() - pMean.X();
            lM.Elt( 1, i ) = iPt->Y() - pMean.Y();
            lM.Elt( 2, i ) = iPt->Z() - pMean.Z();

            i++;
        }
    }

    // Compute covariance matrix
    Matrix lC( 3, 3 );
    {
        lC.LoadValue( 0 );

        // Upper triangular + diagonal
        for( size_t i = 0 ; i < 3 ; i++ )
        for( size_t j = i ; j < 3 ; j++ )
        {
            for( size_t k = 0 ; k < lM.Cols() ; k++ )
                lC.Elt( i, j ) += lM.Elt( i, k ) * lM.Elt( j, k );
        }

        // Lower triangular
        for( size_t i = 0 ; i < 3 ; i++ )
        for( size_t j = 0 ; j < i ; j++ )
            lC.Elt( i, j ) = lC.Elt( j, i );

        // Normalize
        lC.Mult( 1.0 / lNbPoints, lC );
    }

    // Compute SVD
    Matrix lDiag, lVecs;
    EigenJacobi( lC, lDiag, lVecs );

    // Vectors produced by EigenJacobi are already normalized
    pV[0] = Vector3D( lVecs.Elt( 0, 0 ), lVecs.Elt( 1, 0 ), lVecs.Elt( 2, 0 ) );
    pV[1] = Vector3D( lVecs.Elt( 0, 1 ), lVecs.Elt( 1, 1 ), lVecs.Elt( 2, 1 ) );
    pV[2] = Vector3D( lVecs.Elt( 0, 2 ), lVecs.Elt( 1, 2 ), lVecs.Elt( 2, 2 ) );
}

//================================================================================
void WeightedPCA__KindOfBuggy( const Vector<Point3D> & pCloud, const Vector<double> & pWeights, Point3D & pMean, Vector3D pV[3] )
{
//**** DO NOT USE
    LzLogErr("",  "Not a good idea dude!" );
    LzLogErr("",  "This code should not be used!" );
    LzLogErr("",  "C'est vraiment n'imp!" );
    LzLogErr("",  "Seul pMean est fiable, le reste c'est du bullshit!" );
//**** DO NOT USE

    const size_t lNbPoints = pCloud.Size();

    // Check
    if( !lNbPoints )
        LzLogException("",  "Cannot compute PCA on empty points cloud!" );

    // Check
    if( lNbPoints != pWeights.Size() )
        LzLogException("",  "Size mismatch between points and weights! (" << lNbPoints << " != " << pWeights.Size() << ")" );

    // Compute mean point
    pMean.Reset();
    double lTotalWeight = 0;
    for( size_t i = 0 ; i < lNbPoints ; i++ )
    {
        const double lW = pWeights[i];
        const Point3D & lPt = pCloud[i];

        lTotalWeight += lW;

        pMean.X() += lW * lPt.X();
        pMean.Y() += lW * lPt.Y();
        pMean.Z() += lW * lPt.Z();
    }

    // Check total weight ok
    if( LzMath::ToolBox::IsZero( lTotalWeight ) )
        LzLogException("",  "Total point cloud weight is nearly 0! (tot. weight= " << lTotalWeight << ")" );

    // Normalize mean
    pMean /= lTotalWeight;

    // Fill 0-centered matrix
    Matrix lM( 3, lNbPoints );
    for( size_t i = 0 ; i < lNbPoints ; i++ )
    {
        const double lW = pWeights[i];
        const Point3D & lPt = pCloud[i];

        lM.Elt( 0, i ) = lW * ( lPt.X() - pMean.X() );
        lM.Elt( 1, i ) = lW * ( lPt.Y() - pMean.Y() );
        lM.Elt( 2, i ) = lW * ( lPt.Z() - pMean.Z() );
    }

    // Compute covariance matrix
    Matrix lC( 3, 3 );
    {
        lC.LoadValue( 0 );

        // Upper triangular + diagonal
        for( size_t i = 0 ; i < 3 ; i++ )
            for( size_t j = i ; j < 3 ; j++ )
            {
                for( size_t k = 0 ; k < lM.Cols() ; k++ )
                    lC.Elt( i, j ) += lM.Elt( i, k ) * lM.Elt( j, k );
            }

        // Lower triangular
        for(  size_t i = 0 ; i < 3 ; i++ )
            for(  size_t j = 0 ; j < i ; j++ )
                lC.Elt( i, j ) = lC.Elt( j, i );

        // Normalize
        lC.Mult( 1.0 / lNbPoints, lC );
    }

    // Compute SVD
    Matrix lDiag, lVecs;
    EigenJacobi( lC, lDiag, lVecs );

    // Vectors produced by EigenJacobi are already normalized
    pV[0] = Vector3D( lVecs.Elt( 0, 0 ), lVecs.Elt( 1, 0 ), lVecs.Elt( 2, 0 ) );
    pV[1] = Vector3D( lVecs.Elt( 0, 1 ), lVecs.Elt( 1, 1 ), lVecs.Elt( 2, 1 ) );
    pV[2] = Vector3D( lVecs.Elt( 0, 2 ), lVecs.Elt( 1, 2 ), lVecs.Elt( 2, 2 ) );
}

//================================================================================
void WeightedPCA( const List<Point3D> & pCloud, const List<double> & pWeights, Point3D & pMean, Vector3D pV[3] )
{
    Vector<Point3D> lCloud;
    pCloud.ToVector( lCloud );

    Vector<double> lWeights;
    pWeights.ToVector( lWeights );

    WeightedPCA( lCloud, lWeights, pMean, pV );
}

//================================================================================
void WeightedPCA( const Vector<Point3D> & pCloud, const Vector<double> & pWeights, Point3D & pMean, Vector3D pV[3] )
{
    // http://stats.stackexchange.com/questions/113485/weighted-principal-components-analysis

    const size_t lNbPoints = pCloud.Size();

    // Check
    if( lNbPoints < 2 )
        LzLogException("", "Cannot compute PCA on less than 2 points!" );

    // Check
    if( lNbPoints != pWeights.Size() )
        LzLogException("", "Size mismatch between points and weights! ("<<lNbPoints<<" != "<<pWeights.Size()<<")")

    // Compute mean point
    pMean.Reset();
    double lTotalWeight = 0;
    for( size_t i = 0 ; i < lNbPoints ; i++ )
    {
        const double lW = pWeights[i];
        const Point3D & lPt = pCloud[i];

        lTotalWeight += lW;

        pMean.X() += lW * lPt.X();
        pMean.Y() += lW * lPt.Y();
        pMean.Z() += lW * lPt.Z();
    }

    // Check total weight ok
    if( LzMath::ToolBox::IsZero( lTotalWeight ) )
        LzLogException("", "Total point cloud weight is nearly 0! (tot. weight= "<<lTotalWeight<<")")

    // Normalize mean
    pMean /= lTotalWeight;

    // Compute covariance matrix
    Matrix lC( 3, 3 );
#if 1
    {
        Matrix lQ( 3, lNbPoints );
        for( size_t n=0 ; n<lNbPoints ; n++ )
        {
            const Point3D & lPt = pCloud[n];

            lQ.Elt( 0, n ) = lPt.X() - pMean.X();
            lQ.Elt( 1, n ) = lPt.Y() - pMean.Y();
            lQ.Elt( 2, n ) = lPt.Z() - pMean.Z();
        }

        lC.LoadValue( 0 );
        for( int i=0 ; i<3 ; i++ )
        for( int j=0 ; j<3 ; j++ )
        {
            for( size_t n = 0 ; n < lNbPoints ; n++ )
                lC.Elt( i, j ) += lQ.Elt( i, n ) * pWeights[n] * lQ.Elt( j, n );

            lC.Elt( i, j ) /= lTotalWeight;
        }

//lC.Log("C fast");
    }
#else
    {
        Matrix lX( lNbPoints, 3 );
        Matrix lW( lNbPoints, lNbPoints );
        lW.LoadValue( 0 );
        for( size_t i = 0 ; i < lNbPoints ; i++ )
        {
            // X
            const Point3D & lPt = pCloud[i];

            lX.Elt( i, 0 ) = lPt.X() - pMean.X();
            lX.Elt( i, 1 ) = lPt.Y() - pMean.Y();
            lX.Elt( i, 2 ) = lPt.Z() - pMean.Z();

            // W
            lW.Elt( i, i ) = pWeights[i];
        }

        Matrix lXt;
        lX.TransposeTo( lXt );

        Matrix lXtW;
        lXt.Mult( lW, lXtW );

        // C = (Xt W X)/(n - 1)
        lXtW.Mult( lX, lC );
        lC.Mult( 1.0 / ( ( double )lNbPoints - 1.0 ), lC );


//lC.Log("C sloooooooow.....");
    }
#endif

    // Compute SVD
    Matrix lDiag, lVecs;
    EigenJacobi( lC, lDiag, lVecs );

    // Vectors produced by EigenJacobi are already normalized
    pV[0] = Vector3D( lVecs.Elt( 0, 0 ), lVecs.Elt( 1, 0 ), lVecs.Elt( 2, 0 ) );
    pV[1] = Vector3D( lVecs.Elt( 0, 1 ), lVecs.Elt( 1, 1 ), lVecs.Elt( 2, 1 ) );
    pV[2] = Vector3D( lVecs.Elt( 0, 2 ), lVecs.Elt( 1, 2 ), lVecs.Elt( 2, 2 ) );
}

//================================================================================
void PCAMinMax( const vector<Point3D> & pCloud, const Point3D & pO, const Vector3D pV[3], Vector<double> & pMin, Vector<double> & pMax )
{
    Vector<Point3D> lTxVec( pCloud.size() );
    for( size_t i=0 ; i<pCloud.size() ; i++ )
        lTxVec[i] = pCloud[i];

    PCAMinMax( lTxVec, pO, pV, pMin, pMax );
}

//================================================================================
void PCAMinMax( const List<Point3D> & pCloud, const Point3D & pO, const Vector3D pV[3], Vector<double> & pMin, Vector<double> & pMax )
{
    Vector<Point3D> lCloud;
    pCloud.ToVector( lCloud );

    PCAMinMax( lCloud, pO, pV, pMin, pMax );
}

//================================================================================
void PCAMinMax( const Vector<Point3D> & pCloud, const Point3D & pO, const Vector3D pV[3], Vector<double> & pMin, Vector<double> & pMax )
{
    // Check
    if( !pCloud.Size() )
        LzLogException("",  "Cannot compute PCA Min/Max on empty points cloud!" )

    // Reset
    pMin.Assign( 3, +std::numeric_limits<double>::max() );
    pMax.Assign( 3, -std::numeric_limits<double>::max() );

    // Update
    for( size_t i=0 ; i<pCloud.Size() ; i++ )
    {
        const Vector3D lW = pCloud[i] - pO;

        double lScal[3] = { pV[0]*lW, pV[1]*lW, pV[2]*lW };
        for( int j = 0 ; j < 3 ; j++ )
        {
            if( pMin[j] > lScal[j] )
                pMin[j] = lScal[j];

            if( pMax[j] < lScal[j] )
                pMax[j] = lScal[j];
        }
    }
}

//================================================================================
void SVD( const Matrix & pA, Matrix & pU, Matrix & pW, Matrix & pV )
{
    // Check
    const size_t m = pA.Rows();
    const size_t n = pA.Cols();
    if( !m || !n )
        LzLogException("",  "SVD error: empty system!" );

    // Alloc a
    double ** a = new double * [ m + 1 ]; // Indices: 1 to m
    for( size_t i = 1 ; i <= m ; i++ )
        a[i] = new double [ n + 1 ]; // Indices: 1 to n

    // Copy A
    for( size_t i = 0 ; i < m ; i++ )
    for( size_t j = 0 ; j < n ; j++ )
        a[i + 1][j + 1] = pA.Elt( i, j );

    // Alloc w
    double * w = new double [ n + 1 ]; // Indices: 1 to n

    // Alloc v
    double ** v = new double * [ n + 1 ]; // Indices: 1 to n
    for( size_t i = 1 ; i <= n ; i++ )
        v[i] = new double [ n + 1 ]; // Indices: 1 to n

    // Compute
    SVD( a, m, n, w, v );

    // Copy a to U
    pU.SetDims( m, n );
    for( size_t i = 0 ; i < m ; i++ )
    for( size_t j = 0 ; j < n ; j++ )
        pU.Elt( i, j ) = a[i + 1][j + 1];

    // Copy w to W
#if 1
    pW.SetDims( n );
    pW.LoadValue( 0 );
    for( size_t i = 0 ; i < n ; i++ )
        pW.Elt( i ) = w[i + 1];
#else
	// Nope... No need for nxn matrix here!
    pW.SetDims( n, n );
    pW.LoadValue( 0 );
    for( size_t i = 0 ; i < n ; i++ )
        pW.Elt( i, i ) = w[i + 1];
#endif

    // Copy v to V
    pV.SetDims( n, n );
    for( size_t i = 0 ; i < n ; i++ )
    for( size_t j = 0 ; j < n ; j++ )
        pV.Elt( i, j ) = v[i + 1][j + 1];

    // Free a
    for( size_t i = 1 ; i <= m ; i++ )
        delete a[i];
    delete [] a;

    // Free w
    delete [] w;

    // Free v
    for( size_t i = 1 ; i <= n ; i++ )
        delete v[i];
    delete [] v;

#pragma region "Test SVD"
/*
Matrix A( 10, 7 );
A.LoadRandom();
A.Log("A");

Matrix U, W, V;
LzMath::ToolBox::SVD( A, U, W, V );
U.Log("U");
W.Log("W");
V.Log("V");

Matrix UW;
U.Mult( W, UW );
Matrix Vt, UWVt;
V.TransposeTo( Vt );
UW.Mult( Vt, UWVt );
UWVt.Log("UWVt");
*/
#pragma endregion
}

//================================================================================
void SVD( double **a, int m, int n, double *w, double **v )
{
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(a,b) ((a>=b) ? a : b)

    double at, bt, ct;

    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;

    if( m < n )
        LzLogException("", "SVD error: underdetermined system!")

    //rv1=dvector(1,n); /* double *dvector(long nl, long nh) allocate a double vector with subscript range v[nl..nh] */
    rv1 = new double[ n + 1 ];

    for( i = 1; i <= n; i++ )
    {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if( i <= m )
        {
            for ( k = i; k <= m; k++ ) scale += fabs( a[k][i] );
            if ( scale )
            {
                for ( k = i; k <= m; k++ )
                {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN( sqrt( s ), f );
                h = f * g - s;
                a[i][i] = f - g;
                if ( i != n )
                {
                    for ( j = l; j <= n; j++ )
                    {
                        for ( s = 0.0, k = i; k <= m; k++ ) s += a[k][i] * a[k][j];
                        f = s / h;
                        for ( k = i; k <= m; k++ ) a[k][j] += f * a[k][i];
                    }
                }
                for ( k = i; k <= m; k++ ) a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if ( i <= m && i != n )
        {
            for( k = l; k <= n; k++ ) scale += fabs( a[i][k] );
            if( scale )
            {
                for( k = l; k <= n; k++ )
                {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN( sqrt( s ), f );
                h = f * g - s;
                a[i][l] = f - g;
                for( k = l; k <= n; k++ ) rv1[k] = a[i][k] / h;
                if( i != m )
                {
                    for ( j = l; j <= m; j++ )
                    {
                        for( s = 0.0, k = l; k <= n; k++ )
                            s += a[j][k] * a[i][k];

                        for( k = l; k <= n; k++ )
                            a[j][k] += s * rv1[k];
                    }
                }
                for ( k = l; k <= n; k++ )
                    a[i][k] *= scale;
            }
        }
        anorm = MAX( anorm, ( fabs( w[i] ) + fabs( rv1[i] ) ) );
    }
    for ( i = n; i >= 1; i-- )
    {
        if ( i < n )
        {
            if ( g )
            {
                for ( j = l; j <= n; j++ )
                    v[j][i] = ( a[i][j] / a[i][l] ) / g;
                for ( j = l; j <= n; j++ )
                {
                    for ( s = 0.0, k = l; k <= n; k++ ) s += a[i][k] * v[k][j];
                    for ( k = l; k <= n; k++ ) v[k][j] += s * v[k][i];
                }
            }
            for ( j = l; j <= n; j++ ) v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for ( i = n; i >= 1; i-- )
    {
        l = i + 1;
        g = w[i];
        if ( i < n )
            for ( j = l; j <= n; j++ ) a[i][j] = 0.0;
        if ( g )
        {
            g = 1.0 / g;
            if ( i != n )
            {
                for ( j = l; j <= n; j++ )
                {
                    for ( s = 0.0, k = l; k <= m; k++ ) s += a[k][i] * a[k][j];
                    f = ( s / a[i][i] ) * g;
                    for ( k = i; k <= m; k++ ) a[k][j] += f * a[k][i];
                }
            }
            for ( j = i; j <= m; j++ ) a[j][i] *= g;
        }
        else
        {
            for ( j = i; j <= m; j++ ) a[j][i] = 0.0;
        }
        ++a[i][i];
    }
    for ( k = n; k >= 1; k-- )
    {
        for ( its = 1; its <= 30; its++ )
        {
            flag = 1;
            for ( l = k; l >= 1; l-- )
            {
                nm = l - 1;
                if ( ( double )( fabs( rv1[l] ) + anorm ) == anorm )
                {
                    flag = 0;
                    break;
                }
                if ( ( double )( fabs( w[nm] ) + anorm ) == anorm ) break;
            }
            if ( flag )
            {
                c = 0.0;
                s = 1.0;
                for ( i = l; i <= k; i++ )
                {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ( ( double )( fabs( f ) + anorm ) == anorm ) break;
                    g = w[i];
                    h = PYTHAG( f, g );
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = ( -f * h );
                    for ( j = 1; j <= m; j++ )
                    {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];

            if ( l == k )
            {
                if ( z < 0.0 )
                {
                    w[k] = -z;
                    for ( j = 1; j <= n; j++ ) v[j][k] = ( -v[j][k] );
                }
                break;
            }

            if ( its == 50 )
            {
                delete rv1;
                LzLogException("",  "SVD error: too many iterations!" );
            }

            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y );
            g = PYTHAG( f, 1.0 );
            f = ( ( x - z ) * ( x + z ) + h * ( ( y / ( f + SIGN( g, f ) ) ) - h ) ) / x;
            c = s = 1.0;
            for ( j = l; j <= nm; j++ )
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG( f, h );
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for ( jj = 1; jj <= n; jj++ )
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = PYTHAG( f, h );
                w[j] = z;
                if ( z )
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = ( c * g ) + ( s * y );
                x = ( c * y ) - ( s * g );
                for ( jj = 1; jj <= m; jj++ )
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    delete rv1;
    //free_dvector(rv1,1,n);

#undef SIGN
#undef MAX
#undef PYTHAG
}
#pragma endregion


#pragma region "Pressure conversion"
//================================================================================
static const double sPaToUnit[] =
{
    -1,                 // Raw: unused value
    1,                  // 1 Pa = 1 Pa
    0.001,              // 1 Pa = 0.001 kPa
    0.0001,             // 1 Pa = 0.0001 N/cm^2
    0.00001019716213,   // 1 Pa = 0.00001019716213 kg/cm2
    0.01019716213,      // 1 Pa = 0.01019716213 g/cm2
    0.007500615613,     // 1 Pa = 0.007 mmHg
    0.000145037738      // 1 Pa = 0.00145 psi
};

//================================================================================
const char * sPaUnitSuffix[] =
{
    "",
    "Pa",
    "kPa",
    "N/cm^2",
    "kg/cm^2",
    "g/cm^2",
    "mmHg",
    "psi"
};

//================================================================================
double ConvertPa( PaUnit pFrom, PaUnit pTo, double pVal )
{
    if( pFrom == PaUnit::Raw || pTo == PaUnit::Raw )
        LzLogException("",  "Cannot convert pressure values to or from RAW data!" );

    return sPaToUnit[( int )pTo] * pVal / sPaToUnit[( int )pFrom];
}

//================================================================================
std::string PaUnitSuffix( PaUnit pUnit )
{
    return sPaUnitSuffix[( int )pUnit];
}
#pragma endregion


#pragma region "Gradient descent"
/* TEST
{
    std::function<double(const Vector<double> & )> F = []( const Vector<double> & pX )
    {
        // A simple quadratic function
        double lDX = pX[0] - 12;
        double lDY = pX[1] - 34;
        return 5 + 8*lDX*lDX + 5*lDY*lDY;
    };

    Vector<double> lX = {12, 34};
//    Vector<double> lX = {-1000, -1000};
    LzMath::ToolBox::GradResult lRes = LzMath::ToolBox::GradDescent( lX, F, pMaxSteps, pMinGradNorm, 1e-3 );

    if( lRes.mStatus == LzMath::ToolBox::GradResult::Status::StillRunning )
        LzLogWarningMessage("Optimisation failed to converge within the expected accuracy. Res= "<<lRes.toString(), "Warning")
    else
    if( lRes.mStatus == LzMath::ToolBox::GradResult::Status::EnergyRising )
        LzLogWarningMessage("Optimisation failed! Energy rising. Res= "<<lRes.toString(), "Warning")
    else
        LzLogInfoMessage("Optimisation succeeded. Res= "<<lRes.toString(), "Ok")

    // Result
    LzLogM("", "X= "<<lX[0]<<", Y= "<<lX[1])


    return;
}
*/

/* TEST DOMAIN
#if 0
    // Testing domain constrained gradient descent

    // Define domain
    LzMath::ToolBox::GradDomain lDom( 2 );
    lDom.mMins = { 1.0, 2.0 };
    lDom.mMaxs = { 5.0, 6.0 };

    // Define initial guess: must be in domain
    Vector<double> lX;
//    lX = { 0.0, 0.0 }; // Error --- initial guess out of domain
    lX = { 1.5, 2.5 };

    // Define energy
    std::function<double(const Vector<double>&)> Energy = []( const Vector<double> & pX )
        {
#if 0
            // Optimum = 1.23 @ (3.0, 4.0)
            double dX = pX[0] - 3.0;
            double dY = pX[1] - 4.0;
            return 3.0*dX*dX + 5.0*dY*dY + 1.23;
#elif 0
            // Optimum = 4.56 @ (7.0, 3.0) --- out of bounds: should be clamped to (5.0=Max, 3.0)
            double dX = pX[0] - 7.0;
            double dY = pX[1] - 3.0;
            return 5.0*dX*dX + 2.0*dY*dY + 4.56;
#elif 0
            // Optimum = 7.89 @ (10.0, 12.0) --- out of bounds: should be clamped to (5.0=Max, 6.0=Max)
            double dX = pX[0] - 10.0;
            double dY = pX[1] - 12.0;
            return 12.0*dX*dX + 8.0*dY*dY + 7.89;
#elif 0
            // Optimum = 1.11 @ (-1.0, 6.0) --- out of bounds: should be clamped to (1.0=Min, 6.0=Max)
            double dX = pX[0] + 1.0;
            double dY = pX[1] - 6.0;
            return 0.1*dX*dX + 18.0*dY*dY + 1.11;
#elif 0
            // Optimum = 2.22 @ (-1.0, -6.0) --- out of bounds: should be clamped to (1.0=Min, 2.0=Min)
            double dX = pX[0] + 1.0;
            double dY = pX[1] + 6.0;
            return 0.7*dX*dX + 8.0*dY*dY + 2.22;
#elif 1
            // Optimum = 3.33 @ (-1.0, 5.0) --- out of bounds: should be clamped to (1.0=Min, 5.0)
            double dX = pX[0] + 1.0;
            double dY = pX[1] - 5.0;
            return 0.1*dX*dX + 18.0*dY*dY + 3.33;
#endif
        };

    LzMath::ToolBox::GradResult lRes = LzMath::ToolBox::GradDescent( LzMath::ToolBox::GradLog::AllPara, lX, Energy, 2000, 1e-5, lDom );

    // Log res
    LzLogM("", "X= "<<lX.ToString0())

    // Info
    if( lRes.mStatus == LzMath::ToolBox::GradResult::Status::StillRunning )
        LzLogWarningMessage("Optimisation failed to converge within the expected accuracy. Res= "<<lRes.toString(), "Warning")
    else
    if( lRes.mStatus == LzMath::ToolBox::GradResult::Status::EnergyRising )
        LzLogWarningMessage("Optimisation failed! Energy rising. Res= "<<lRes.toString(), "Warning")
    else
        LzLogInfoMessage("Optimisation succeeded. Res= "<<lRes.toString(), "Ok")

    // *** Boum!
    LzLogException("", "***");
#endif
*/

//================================================================================
GradDomain::GradDomain( size_t pDim )
 : mMinMaxMargin(1e-6)
{
    // Check
    if( pDim == 0 )
        LzLogException("", "Cannot define a domain in dimension 0!")

    // Init min and max to NON BOUND
    mMins.Assign( pDim,  0.0 );
    mMaxs.Assign( pDim, -1.0 );
}

//================================================================================
bool GradDomain::IsBound( size_t pDim ) const
{
    // Check
    if( pDim>=mMins.Size() || pDim>=mMaxs.Size() )
        LzLogException("", "Invalid dimension "<<pDim<<"!")

    // Is bound?
    return mMins[pDim] <= mMaxs[pDim];
}

//================================================================================
bool GradDomain::NearMin( const Vector<double> & pX, size_t pDim ) const
{
    // Check
    if( pDim>=pX.Size() || pDim>=mMins.Size() )
        LzLogException("", "Invalid dimension "<<pDim<<"!")

    return pX[pDim] < mMins[pDim] + mMinMaxMargin;
}

//================================================================================
bool GradDomain::NearMax( const Vector<double> & pX, size_t pDim ) const
{
    // Check
    if( pDim>=pX.Size() || pDim>=mMaxs.Size() )
        LzLogException("", "Invalid dimension "<<pDim<<"!")

    return pX[pDim] > mMaxs[pDim] - mMinMaxMargin;
}

//================================================================================
bool GradDomain::IsOutOfBounds( const Vector<double> & pX ) const
{
    // Check
    if( pX.Size() == 0 )
        LzLogException("", "Cannot check bounds in dimension 0!")

    // Check dimensions
    if( mMins.Size()!=pX.Size() || mMaxs.Size()!=pX.Size() )
        LzLogException("", "Size mismatch between X, Mins and/or Maxs! X is "<<pX.Size()<<", Mins is "<<mMins.Size()<<", Maxs is "<<mMaxs.Size()<<".")

    // Check that initial guess is within domain, if it is bound
    for( size_t d=0 ; d<pX.Size() ; d++ )
    {
        // Bound, and out of bounds?
        if( IsBound(d) && (pX[d]<mMins[d] || pX[d]>mMaxs[d]) )
            return true;
    }

    // Within bounds
    return false;
}

//================================================================================
double GradDomain::Clamp( const Vector<double> & pA, const Vector<double> & pB ) const
{
    // Check
    if( pA.Size() != pB.Size() )
        LzLogException("", "Size mismatch between A and B!")

    // Check
    if( pA.Size()!=mMins.Size() || pA.Size()!=mMaxs.Size() )
        LzLogException("", "Size mismatch between A and B, and/or Mins or Maxs!")

//LzLogN("", "Clamping")

    // Assess all dimensions
    double lSqueeze = 1.0;
    for( size_t d=0 ; d<pA.Size() ; d++ )
    {
        // Check if bound
        if( !IsBound(d) )
            continue;

        // Squeezed B
        double lSzB = pA[d] + lSqueeze*(pB[d] - pA[d]);
//LzLogM("", "lSqueeze= "<<lSqueeze<<", lSzB= "<<lSzB)

        // Check if squeezed B is below Min
        if( lSzB < mMins[d] )
        {
            // In case the input assumption is violated
            if( IsZero(pA[d] - lSzB, mMinMaxMargin) )
            {
                // SzB is out of bounds and so is A! (at least for dimension d)
                LzLogErr("", "Input assumption violated! A is nearly outside the optimization domain.")
                continue;
            }

            // Adjust squeeze factor
//LzLogM("", "(pA["<<d<<"] - lSzB)= "<<(pA[d] - lSzB))
            lSqueeze *= (pA[d] - mMins[d]) / (pA[d] - lSzB);
//LzLogM("", "Sq= "<<lSqueeze)
        }
        else
        // Check if squeezed B is above Max
        if( lSzB > mMaxs[d] )
        {
            // In case the input assumption is violated
            if( IsZero(pA[d] - lSzB, mMinMaxMargin) )
            {
                // SzB is out of bounds and so is A! (at least for dimension d)
                LzLogErr("", "Input assumption violated! A is nearly outside the optimization domain.")
                continue;
            }

            // Adjust squeeze factor
//LzLogM("", "(lSzB - pA["<<d<<"])= "<<(lSzB - pA[d]))
            lSqueeze *= (mMaxs[d] - pA[d]) / (lSzB - pA[d]);
//LzLogM("", "Sq= "<<lSqueeze)
        }
    }

    return lSqueeze;
}

//================================================================================
//
// pX is modified in the function but its initial value is restored at the end
//
static void FinDiffGrad( std::function<double(const Vector<double> &)> pEnergy, Vector<double> & pX,  const GradDomain & pDom,
                         double pFiniteDiffStep,
                         Vector<double> & pGrad, double & pGradNorm, bool pCentered=true, double pE0=-1/*needed only if forward diff and not centered*/ )
{
    // Check
    if( !pCentered && pE0<0 )
        LzLogException("", "Local energy must be provided if using forward differences!")

    // Read dimension
    const size_t lDim = pX.Size();

    // Reset
    pGrad.Assign( lDim, 0.0 );
    pGradNorm = 0.0;

    // Compute
    for( size_t d=0 ; d<lDim ; d++ )
    {
        if( pCentered )
        {
            // Centered differences
            //=====================

            // Step forward
            pX[d] += pFiniteDiffStep;
            const double lE_Fw = pEnergy(pX);

            // Step backward
            pX[d] -= 2 * pFiniteDiffStep;
            const double lE_Bk = pEnergy(pX);

            // Step forward, back to starting point
            pX[d] += pFiniteDiffStep;

            //==> Compute gradient
            pGrad[d] = (lE_Fw - lE_Bk) / (2.0 * pFiniteDiffStep);
        }
        else
        {
            // Forward differences
            //====================

            // Shift forward
            pX[d] += pFiniteDiffStep;
            const double lE_Fw = pEnergy(pX);

            // Step backward, back to starting point
            pX[d] -= pFiniteDiffStep;

            //==> Compute gradient
            pGrad[d] = (lE_Fw - pE0) / pFiniteDiffStep;
        }

        // Clamp
        if( pDom.IsBound(d) )
        {
            // Check if near min or max and cancel gradient if needed
            if( (pDom.NearMin(pX, d) && pGrad[d]>0.0) /*we are going apposite to the gradient*/
             || (pDom.NearMax(pX, d) && pGrad[d]<0.0) /*we are going apposite to the gradient*/ )
            {
                // Cancel gradient, variable not moving
                pGrad[d] = 0.0;
            }
        }

        // Update norm
        pGradNorm += pGrad[d]*pGrad[d];
    }

    // Take root of square norm
    pGradNorm = sqrt( pGradNorm );
}

//================================================================================
GradResult GradDescent( GradLog pLogMode, Vector<double> & pX, std::function<double(const Vector<double> &)> pEnergy,
                        size_t pMaxSteps, double pMinGradNorm, const GradDomain & pDom, double pFiniteDiffStep/*=1e-5*/ )
{
    // Log
    LzLogN("", "Computing gradient descent. Max steps= "<<pMaxSteps<<", Min grad norm= "<<pMinGradNorm<<".")

    // Problem dimension
    const size_t lDim = pX.Size();

    // Check
    if( lDim == 0 )
        LzLogException("", "Cannot perform gradient descent in dimension 0!")

    // Check
    if( pDom.IsOutOfBounds(pX) )
        LzLogException("", "Initial guess is out of prescribed domain!")

    // Check
    if( pFiniteDiffStep < 1e-12 )
        LzLogException("", "Cannot perform gradient descent using a finite difference step= "<<pFiniteDiffStep<<"!")

    // Gradient and energy at the end of a setp
    Vector<double> lGrad( lDim );
    double lGradNorm;
    double E1 = 666.666;
    for( size_t iStep=0 ; iStep<pMaxSteps ; iStep++ )
    {
        // Light log
        if( iStep%100 == 0 )
            LzLogM("", "--- Step: "<<iStep)

        // Log
        if( pLogMode == GradLog::All )
            LzLogM("", "--- Step: "<<iStep)

        //---------------------------------
        // Compute E0 from scratch
        // or from E1
        //---------------------------------

        const double E0 = iStep==0 ? pEnergy(pX) : E1 ;

        // Log
        if( pLogMode==GradLog::All || pLogMode==GradLog::AllPara )
            LzLogM("", "E0= "<<E0)

        //---------------------------------
        // Compute gradient
        //---------------------------------

/*
#if 1
*/
//        FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, lGrad, lGradNorm );
        FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, lGrad, lGradNorm, false, E0 );
/*
#else
        // Reset gradient and norm
        lGrad.SetAll( 0.0 );
        lGradNorm = 0.0;
        {
            // Finite diff sample
//            Vector<double> lX_dX = pX;
            //
            for( size_t d=0 ; d<lDim ; d++ )
            {
//gradient mode : centered, forward
                if( 1 )
                {
                    // Centered differences
                    //=====================

                    // Step forward
//lX_dX[d] += pFiniteDiffStep;
pX[d] += pFiniteDiffStep;
                    //
//const double lE_Fw = pEnergy(lX_dX);
const double lE_Fw = pEnergy(pX);

                    // Step backward
//lX_dX[d] -= 2 * pFiniteDiffStep;
pX[d] -= 2 * pFiniteDiffStep;
                    //
//const double lE_Bk = pEnergy(lX_dX);
const double lE_Bk = pEnergy(pX);

                    // Compute gradient
                    lGrad[d] = (lE_Fw - lE_Bk) / (2 * pFiniteDiffStep);

                    // Step forward
//lX_dX[d] += pFiniteDiffStep;
pX[d] += pFiniteDiffStep;
                }
                else
                {
                    // Forward differences
                    //====================

                    // Shift forward
//lX_dX[d] += pFiniteDiffStep;
pX[d] += pFiniteDiffStep;

                    // Compute gradient at d
//lGrad[d] = (pEnergy(lX_dX) - E0) / pFiniteDiffStep;
lGrad[d] = (pEnergy(pX) - E0) / pFiniteDiffStep;

                    // Step backward
//lX_dX[d] -= pFiniteDiffStep;
pX[d] -= pFiniteDiffStep;
                }

                // Clamp
                if( pDom.IsBound(d) )
                {
                    // Check if near min or max and cancel gradient if needed
                    if( (pDom.NearMin(pX, d) && lGrad[d]>0.0) /*we are going apposite to the gradient* /
                     || (pDom.NearMax(pX, d) && lGrad[d]<0.0) /*we are going apposite to the gradient* / )
                    {
                        // Cancel gradient, variable not moving
                        lGrad[d] = 0.0;
                    }
                }

                // Update norm
                lGradNorm += lGrad[d]*lGrad[d];
            }

            // Take root of square norm
            lGradNorm = sqrt( lGradNorm );
        }
#endif
*/
        // Check
        if( lGradNorm <= pMinGradNorm )
        {
            LzLogM("", "Grad norm= "<<lGradNorm<<" is below min= "<<pMinGradNorm<<". Descent finished after "<<iStep<<" step(s).")
            return GradResult( GradResult::Status::Converged, iStep, lGradNorm, E0 );
        }

        //---------------------------------
        // Compute energy at max step
        // Emax
        //---------------------------------

//******
//
const double FULL_STEP_LENGTH = 1.0;
//
//******

// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!
// **** TAKE multiple steps !!!

//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!
//     ***********   ou algo Secant ici !!

        // Compute position at full step
        Vector<double> lXmax = pX;
        for( size_t d=0 ; d<lDim ; d++ )
            lXmax[d] -= FULL_STEP_LENGTH * lGrad[d] / lGradNorm;

        // Adjust Xmax if it is out of bounds
        // Alpha in ]0, 1[ = adjustment factor
        double lAlpha = 1.0;
        if( pDom.IsOutOfBounds(lXmax) )
        {
            // Find alpha and move Xmax
            lAlpha = pDom.Clamp( pX, lXmax );

            // Adjust Xmax
            for( size_t d=0 ; d<lDim ; d++ )
                lXmax[d] = pX[d] + lAlpha*(lXmax[d] - pX[d]);
        }

        // Compute energy
//auto lStart = Chrono_Now;
        const double Emax = pEnergy( lXmax );
//LzLogM("", "Emax computed in "<<Chrono_MSecs(lStart)<<" msec.")

        //---------------------------------
        // Compute parabola
        //---------------------------------

        const double C = E0; //********************** inutile sauf pour debug
        const double B = -(lAlpha * FULL_STEP_LENGTH) * lGradNorm;
        const double A = Emax - B - E0;

        // Log
        if( pLogMode == GradLog::All )
        {
            LzLogM("", "Alpha= "<<lAlpha)
            LzLogM("", "A= "<<A<<", B= "<<B<<", C= "<<C)
        }

        // Log
        if( pLogMode == GradLog::AllPara )
        {
static int count = 0;
            LzLogN("", "===== PARABOLA CHECK ===== "<<count)
            LzLogM("", "t\tEnergy\tParabola")
            LzLogM("", "A= "<<A<<", B= "<<B<<", C= "<<C)
            for( double t=0.0 ; t<1.05 ; t+=0.05 )
            {
                double P = A*t*t + B*t + C;
                double E;
                {
                    Vector<double> lS = pX;
                    for( size_t d=0 ; d<lDim ; d++ )
                        lS[d] -= t * (lAlpha * FULL_STEP_LENGTH) * lGrad[d] / lGradNorm;
                    E = pEnergy( lS );
                }

                // Log
                LzLogM("", t<<"\t"<<E<<"\t"<<P)
            }
            LzLogM("", "==========================")

if( count == 400 )
    LzLogInfoMessage("", "Count reached")

count++;
        }

        //--------------------------------
        // Find minimal energy
        // E1
        //--------------------------------

        // Find minimum in [0,Tmax]
        double t1; //******************************************************* USEFUL???????
        if( A < 1e-6 )
        {
            // Linear decrease or downwards parabola
            t1 = 1.0;
            pX = lXmax;
            E1 = Emax;
        }
        else
        {
            // A > 0 and B < 0 => t1 > 0

            // Upwards parabola
            t1 = -B/(2*A);

            // Clamp
            if( t1 > 1.0 )
            {
                t1 = 1.0; //*************************************** when does that happen ????
                pX = lXmax;
                E1 = Emax;
            }
            else
            {
                for( size_t d=0 ; d<lDim ; d++ )
                    pX[d] -= t1 * (lAlpha * FULL_STEP_LENGTH) * lGrad[d] / lGradNorm;

                // Compute energy at min point
                E1 = pEnergy( pX );

                // Check if no better energy at Emax
                if( E1 > Emax )
                {
//                    //*** see what has happened when we hit this BP
//                    TxLogMsg("Found better energy at Tmax!");
                    t1 = 1.0;
                    pX = lXmax;
                    E1 = Emax;
                }
            }
        }


//********************************** check that we are remaining in the domain
if( pDom.IsOutOfBounds(pX) )
    LzLogException("", "Somehow the solution escaped the prescribed domain!")
//********************************** check that we are remaining in the domain


        // Log
        if( pLogMode==GradLog::All || pLogMode==GradLog::AllPara )
        {
//            LzLogM("", "grad norm= "<<lGradNorm)
            LzLogM("", "t1= "<<t1)
//            LzLogM("", "E1= "<<E1)
        }

        //--------------------------------
        // Check for abnormal energy rise
        //--------------------------------

        if( E0 < E1 )
        {
//            //***** configuration easy to hit with a 2x2x2cube !
//            //******** Log energy and parabola curves here to see what happened

//            //------------------------------------------
//            //*** pour sortir du trou :
//            //**************** diminuer MAX_STEP si echec, puis re-augmenter qd succes !
//            //------------------------------------------

//            // Perf stop
////            mMsecsPerStep = TxServices::Log::PerfStopMsec(lPerfStart);
            LzLogE("", "Energy rising! Current energy E0= "<<E0<<" < next energy E1= "<<E1<<".");
//            return GradResult( GradResult::Status::EnergyRising, iStep, lGradNorm, E0 );
        }

        // Log
        if( pLogMode==GradLog::All || pLogMode==GradLog::AllPara )
            LzLogM("", "---")
    }

    return GradResult( GradResult::Status::StillRunning, pMaxSteps+1, lGradNorm, E1 );
}
#pragma endregion


#pragma region "Conjugate gradient"
//================================================================================
GradResult ConjugateGradient( GradLog /*pLogMode*/, Vector<double> & pX, std::function<double(const Vector<double> &)> pEnergy,
                              size_t pMaxSteps, double pMinGradNorm, const GradDomain & pDom, double pFiniteDiffStep/*=1e-5*/ )
{
    // Log
    LzLogN("", "Computing conjugate gradient descent. Max steps= "<<pMaxSteps<<", Min grad norm= "<<pMinGradNorm<<".")

    // Problem dimension
    const size_t lDim = pX.Size();

    // Check
    if( lDim == 0 )
        LzLogException("", "Cannot perform gradient descent in dimension 0!")

    // Check
    if( pDom.IsOutOfBounds(pX) )
        LzLogException("", "Initial guess is out of prescribed domain!")

    // Check
    if( pFiniteDiffStep < 1e-12 )
        LzLogException("", "Cannot perform gradient descent using a finite difference step= "<<pFiniteDiffStep<<"!")

#if 1
    size_t i = 0; // Step
    size_t k = 0; // Counter for restarting CG (forgetting previous directions)

    Vector<double> r; // Residual
    double unused;
    FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, r, unused );

    // Opposite direction
    for( size_t i=0 ; i<lDim ; i++ )
        r[i] = -r[i];

    //
    // *** No preconditioning here
    //
    Vector<double> s = r;

    Vector<double> d = r; // Conj direction

    // Stop criterion based on norm of residual (****** CHECK
    double delta_new = r * d;
//    double delta_0 = delta_new;

    while( i<pMaxSteps && sqrt(delta_new)>pMinGradNorm )
    {
        // Log
        LzLogN("", "Step "<<i)

        size_t j = 0; // Secant iterations counter
        double delta_d = d * d;

        double sigma_0 = 1.0 / sqrt(delta_d); //*********** CHECK
        double alpha = -sigma_0;

        double eta_prev;
        {
            Vector<double> x_sigma_d = pX;
            for( size_t i=0 ; i<lDim ; i++ )
                x_sigma_d[i] += sigma_0 * d[i];

            Vector<double> df;
            double unused;
            FinDiffGrad( pEnergy, x_sigma_d, pDom, pFiniteDiffStep, df, unused );

            eta_prev = df * d;
        }

        const size_t j_max = 10; //******************* CHECK
        do
        {
            double eta;
            {
                Vector<double> df;
                double unused;
                FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, df, unused );

                eta = df * d;
            }

            alpha = alpha * eta / (eta_prev - eta);

            for( size_t i=0 ; i<lDim ; i++ )
                pX[i] += alpha * d[i];

// **** dans boucle Secant : verifier qu'on ne sort pas du domaine

            eta_prev = eta;

            j++;
        }
        while( j<j_max && alpha*alpha*delta_d>1e-6 ); //******************* CHECK

//*** New energy at x
LzLogM("", "Energy= "<<pEnergy(pX))

        // r <= -f'( x )
        {
            double unused;
            FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, r, unused );

            // Opposite direction
            for( size_t i=0 ; i<lDim ; i++ )
                r[i] = -r[i];
        }

        double delta_old = delta_new;
        double delta_mid = r * s;

        //
        // *** No preconditioning here
        //
        s = r;

        delta_new = r * s;

LzLogM("", "New grad norm= "<<sqrt(delta_new))

        double beta = (delta_new - delta_mid) / delta_old;

        k++;

        if( k==lDim || beta<=0.0 )
        {
            d = s;
            k = 0;
        }
        else
        {
            for( size_t i=0 ; i<lDim ; i++ )
                d[i] = s[i] + beta * d[i];
        }

        i++;
    }


    LzLogM("", "*** clean up finish log")


//    return GradResult( GradResult::Status::StillRunning, pMaxSteps+1, norm_r_i0, pEnergy(pX) );
    return GradResult( GradResult::Status::StillRunning, pMaxSteps+1, 666, pEnergy(pX) );

#else
    //------------
    // Init
    //------------

    Vector<double> d_i0( lDim );
    Vector<double> r_i0( lDim );
    double norm_r_i0;
    {
        // Log
        LzLogN("", "Init")

        FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, r_i0, norm_r_i0 );

        // Point to pposite direction
        for( size_t d=0 ; d<lDim ; d++ )
        {
            r_i0[d] *= -1.0;
            d_i0[d] = r_i0[d];
        }

        // Check
        if( norm_r_i0 <= pMinGradNorm )
        {
            LzLogM("", "Grad norm= "<<norm_r_i0<<" is below min= "<<pMinGradNorm<<". Descent finished after init.")
            return GradResult( GradResult::Status::Converged, 0, norm_r_i0, pEnergy(pX) );
        }
    }

    //------------
    // Loop
    //------------

    Vector<double> r_i1( lDim );
    double sigma_i0 = 1.0 / norm_r_i0; //************ EVITER DE FAIRE UN PAS TROP GRAND !!!

    for( size_t iStep=0 ; iStep<pMaxSteps ; iStep++ )
    {
        // Log
        LzLogN("", "Step "<<iStep)

        /* Here:
         *
         * r_i0: residual at previous iteration
         * d_i0: conjugate direction at previous iteration
         *
         *
         */

        // Line search for alpha --> min pEnergy( pX + alpha*d_i0 )
//double sigma_i1;
        double alpha_i0;
        {
//***** ne pas recalculer ce truc !! voir si pas dispo par ailleurs
            Vector<double> grad_at_X;
            {
                double unused;
                FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, grad_at_X, unused );


//adapter a norme du grad a chaque pas
sigma_i0 = 1.0 / unused; //************ EVITER DE FAIRE UN PAS TROP GRAND !!!

LzLogM("", "grad norm= "<<unused)

            }
//***** ne pas recalculer ce truc !! voir si pas dispo par ailleurs

            // f'( pX + sigma * d_i0 )
            Vector<double> grad_at_X_plus_sigma_d_i0;
            {
                // Shift pX
                Vector<double> X_plus_sigma_d_i0 = pX;
                for( size_t d=0 ; d<lDim ; d++ )
                    X_plus_sigma_d_i0[d] += sigma_i0 * d_i0[d];

                // Compute gradient
                double unused;
                FinDiffGrad( pEnergy, X_plus_sigma_d_i0, pDom, pFiniteDiffStep, grad_at_X_plus_sigma_d_i0, unused );
            }

            // Alpha = -sigma * ( [f'(pX)]^T d_i0 ) / ( [f'(pX+sigma d_i0)]^T d_i0 - [f'(pX)]^T d_i0 )



            alpha_i0 = -sigma_i0 * (grad_at_X * d_i0) / (grad_at_X_plus_sigma_d_i0 * d_i0 - grad_at_X * d_i0);

            //            si les deux grads sont les memes ==> que fait-on ?
            //            si les deux grads sont les memes ==> que fait-on ?
            //            si les deux grads sont les memes ==> que fait-on ?
            //            si les deux grads sont les memes ==> que fait-on ?

//// Store sigma_i1 for next iteration
//sigma_i1 = -alpha_i0;


            // Apply displacement along conjugate direction
            for( size_t d=0 ; d<lDim ; d++ )
                pX[d] = pX[d] + alpha_i0*d_i0[d];
        }

        // Update residual
        double norm_r_i1;
        FinDiffGrad( pEnergy, pX, pDom, pFiniteDiffStep, r_i1, norm_r_i1 );

        // Check
        if( norm_r_i1 <= pMinGradNorm )
        {
            LzLogM("", "Grad norm= "<<norm_r_i1<<" is below min= "<<pMinGradNorm<<". Descent finished after "<<iStep<<" step(s).")
            return GradResult( GradResult::Status::Converged, 0, norm_r_i1, pEnergy(pX) );
        }

//***************************************************** utiliser minus_r_i1 et ne pas changer la direction
//***************************************************** IDEM pour d : mais ne faire la modif qu'en V2 : apres verif du code !!!
        // Point to pposite direction
        for( size_t d=0 ; d<lDim ; d++ )
            r_i1[d] *= -1.0;
//***************************************************** utiliser minus_r_i1 et ne pas changer la direction
//***************************************************** IDEM pour d : mais ne faire la modif qu'en V2 : apres verif du code !!!

        // Compute beta
        double beta_i1;
        if( 1 )
        {
            // Fletcher-Reeves
            // ---------------
            //
            // beta_i1 = r_i1^T r_i1 / r_i0^T r_i0
            //
            // r_i0^T r_i0 = norm_r_i0^2
            // r_i1^T r_i1 = norm_r_i1^2
            //
            beta_i1 = (norm_r_i1 * norm_r_i1) / (norm_r_i0 * norm_r_i0);
        }
        else
        {
            // Polak-Ribiere
            // -------------
            //
            //
            beta_i1 = -666;

            LzLogException("", "TODO: Polak-Ribiere")
        }

        // Update conjugate direction
        for( size_t d=0 ; d<lDim ; d++ )
            //***************************************************** utiliser minus_r_i1 et ne pas changer la direction
            d_i0[d] = r_i1[d] + beta_i1*d_i0[d];
            //***************************************************** utiliser minus_r_i1 et ne pas changer la direction

        // Commit r_i1 ==> r_i0
        r_i0 = r_i1;
        norm_r_i0 = norm_r_i1;


        LzLogM("", "E= "<<pEnergy(pX))
//// Commit sigma
//sigma_i0 = sigma_i1;
    }

    return GradResult( GradResult::Status::StillRunning, pMaxSteps+1, norm_r_i0, pEnergy(pX) );
#endif
}
#pragma endregion


#pragma region "Registration"
//================================================================================
void GetCartesianPermutations( vector<RigidTr3D> & pTo,
                               const List<int> & pForbiddenXs/*={}*/,
                               const List<int> & pForbiddenYs/*={}*//*,
                               const List<int> & pForbiddenZs={}*/ )
{
    // Clear previous if any
    pTo.clear();

    //
    // Produce all alias matrices R'(X', Y', Z') to R(X, Y, Z)
    //
    // X' can be mapped onto +/-X, +/-Y, or +/-Z
    for( int ix=0 ; ix<3 ; ix++ ) // Index of X': 0=X, 1=Y, 2=Z
    for( int sx=0 ; sx<2 ; sx++ ) // Sign of X': 0=-1, 1=+1
    {
        // Build X'
        Vector3D X_p( 0, 0, 0 );
        X_p[ix] = sx==0 ? -1 : +1 ;

        // Skip forbidden Xs
        if( pForbiddenXs.FindElt(X_p[ix] * ix) )
        {
            LzLogM("", "GetCartesianPermutations: skipping X' --> "<<X_p.ToString()<<".")
            continue;
        }

        // Y' can be mapped onto +/-X, +/-Y, or +/-Z
        for( int iy=0 ; iy<3 ; iy++ ) // Index of Y': 0=X, 1=Y, 2=Z
        {
            // ... but it must not be the same axis as X'
            if( iy == ix )
                continue;

            // Explore sign for admissible Y' axis
            for( int sy=0 ; sy<2 ; sy++ ) // Sign of Y': 0=-1, 1=+1
            {
                // Build Y'
                Vector3D Y_p( 0, 0, 0 );
                Y_p[iy] = sy==0 ? -1 : +1 ;

                // Skip forbidden Ys
                if( pForbiddenYs.FindElt(Y_p[iy] * iy) )
                {
                    LzLogM("", "GetCartesianPermutations: skipping Y' --> "<<Y_p.ToString()<<".")
                    continue;
                }

/*
                // Compute Z_p to check if forbidden or not
                const Vector3D Z_p = X_p^Y_p;

                // Skip forbidden Zs
                if( pForbiddenZs.FindElt(666) ) // ******************************************************************
                    continue;
*/

                // Stash transform
                RigidTr3D lAlias_Rp_to_R( X_p, Y_p, Point3D(0,0,0) );
                pTo.push_back( lAlias_Rp_to_R );
            }
        }
    }

    // Final count
    LzLogM("", "Found "<<pTo.size()<<" permutation(s).")
}

//================================================================================
RigidTr3D Arun( const Vector<Point3D> & pP1, const Vector<Point3D> & pP2 )
{
    vector<Point3D> lP1( pP1.Size() );
    for( size_t i = 0 ; i < pP1.Size() ; i++ )
        lP1[i] = pP1[i];

    vector<Point3D> lP2( pP2.Size() );
    for( size_t i = 0 ; i < pP2.Size() ; i++ )
        lP2[i] = pP2[i];

    return Arun( lP1, lP2 );
}

//================================================================================
void LogArunError( const RigidTr3D & pAlibi_P1_to_P2, const Vector<Point3D> & pP1, const Vector<Point3D> & pP2 )
{
    // Check
    const size_t N = pP1.Size();
    if( N < 3 || N != pP2.Size() )
        LzLogException("", "Arun registration check impossible: " << pP1.Size() << " points vs " << pP2.Size() << "points!" )

    // Compute errors
    Vector<double> lErrs( pP1.Size() );
    for( size_t i = 0 ; i < N ; i++ )
        lErrs[i] = pP2[i].DistanceTo( pAlibi_P1_to_P2 * pP1[i] );

    // Compute errors
    double lMean, lMax, lStdDev;
    MeanMaxStdDev( lErrs, lMean, lMax, lStdDev );

    // Log
    LzLogMsg("", "Arun error. Mean= " << lMean << ", Max.= " << lMax << ", std. dev.= " << lStdDev << "." )
}

//================================================================================
RigidTr3D ArunLoop( const vector<Point3D> & pP1, const vector<Point3D> & pP2, double pMaxEulerRotDeg/*=0.1*/, double pMaxTransNorm/*=0.1*/, size_t pMaxSteps/*=300*/ )
{
    Vector<Point3D> lP1( (size_t)pP1.size() );
    for( size_t i = 0 ; i < pP1.size() ; i++ )
        lP1[i] = pP1[i];

    Vector<Point3D> lP2( (size_t)pP2.size() );
    for( size_t i = 0 ; i < pP2.size() ; i++ )
        lP2[i] = pP2[i];

    return ArunLoop( lP1, lP2, pMaxEulerRotDeg, pMaxTransNorm, pMaxSteps );
}

//================================================================================
RigidTr3D ArunLoop( const Vector<Point3D> & pP1, const Vector<Point3D> & pP2, double pMaxEulerRotDeg/*=0.1*/, double pMaxTransNorm/*=0.1*/, size_t pMaxSteps/*=300*/ )
{
    // Copy source points
    Vector<Point3D> lP1 = pP1;

    // Compute ICP
    RigidTr3D lICP;
    Vector<Point3D> lNearest( lP1.Size() );
    size_t iStep = 0;
    while( iStep < pMaxSteps )
    {
        // Compute registration at this step
        RigidTr3D lArun;
        {
            // Find nearest points in destination set
            for( size_t p1 = 0 ; p1 < lP1.Size() ; p1++ )
            {
                const Point3D & A = lP1[p1];

                double lMinDist = +std::numeric_limits<double>::max();
                for( size_t p2 = 0 ; p2 < pP2.Size() ; p2++ )
                {
                    const Point3D & B = pP2[p2];

                    double lDist = A.DistanceTo( B );
                    if( lMinDist > lDist )
                    {
                        lMinDist = lDist;
                        lNearest[p1] = B;
                    }
                }
            }

            // Compute transform
            lArun = Arun( lP1, lNearest );

            // Update source points
            for( size_t p1 = 0 ; p1 < lP1.Size() ; p1++ )
                lP1[p1] *= lArun;
        }

        // Accumulate
        lICP = lArun * lICP;

        // Stop?
        {
            double lDegX, lDegY, lDegZ;
            lArun.GetEulerAngles( lDegX, lDegY, lDegZ );
            const double lAbsDegX = fabs( lDegX );
            const double lAbsDegY = fabs( lDegY );
            const double lAbsDegZ = fabs( lDegZ );
            const Vector3D & lTr = lArun.GetTranslation();

// DEBUG
//LzLogMsg("", "X deg.= "+lAbsDegX+", Y deg.="+lAbsDegY+", Z deg.= "+lAbsDegZ+", Tr= "+lTr.ToString()+".");

            if( lAbsDegX < pMaxEulerRotDeg
			 && lAbsDegY < pMaxEulerRotDeg
			 && lAbsDegZ < pMaxEulerRotDeg
			 && lTr.Norm() < pMaxTransNorm )
            {
                // Log
                LzLogMsg("",  "ICP converged in " << iStep << " iteration(s)." );

                return lICP;
            }
        }

        // Next step
        iStep++;
    }

    // Error
    LzLogException("",  "ICP failed to converge after " << iStep << " steps!" );
}

//================================================================================
RigidTr3D ArunLoop( const vector<Point3D> & pP1, const Tree3D & pDstTree, double pMaxEulerRotDeg/*=0.1*/, double pMaxTransNorm/*=0.1*/, size_t pMaxSteps/*=300*/ )
{
    Vector<Point3D> lP1( (size_t)pP1.size() );
    for( size_t i = 0 ; i < pP1.size() ; i++ )
        lP1[i] = pP1[i];

    return ArunLoop( lP1, pDstTree, pMaxEulerRotDeg, pMaxTransNorm, pMaxSteps );
}

//================================================================================
RigidTr3D ArunLoop( const Vector<Point3D> & pP1, const Tree3D & pDstTree, double pMaxEulerRotDeg/*=0.1*/, double pMaxTransNorm/*=0.1*/, size_t pMaxSteps/*=300*/ )
{
    // Copy source points
    Vector<Point3D> lP1 = pP1;

    // Compute ICP
    RigidTr3D lICP;
    Vector<Point3D> lNearest( lP1.Size() );
    size_t iStep = 0;
    while( iStep < pMaxSteps )
    {
        // Compute registration at this step
        RigidTr3D lArun;
        {
            // Find nearest points in destination set
            for( size_t p1 = 0 ; p1 < lP1.Size() ; p1++ )
            {
                const Point3D & A = lP1[p1];
                lNearest[p1] = pDstTree.OldPoint3D( pDstTree.GetNearestVertex(A, Tree3D::NearestMode::Accurate) );

//double lMinDist = +std::numeric_limits<double>::max();
//for( size_t p2 = 0 ; p2 < pP2.Size() ; p2++ )
//{
//    const Point3D & B = pP2[p2];
//
//    double lDist = A.DistanceTo( B );
//    if( lMinDist > lDist )
//    {
//        lMinDist = lDist;
//        lNearest[p1] = B;
//    }
//}
            }

            // Compute transform
            lArun = Arun( lP1, lNearest );

            // Update source points
            for( size_t p1 = 0 ; p1 < lP1.Size() ; p1++ )
                lP1[p1] *= lArun;
        }

        // Accumulate
        lICP = lArun * lICP;

        // Stop?
        {
            double lDegX, lDegY, lDegZ;
            lArun.GetEulerAngles( lDegX, lDegY, lDegZ );
            const double lAbsDegX = fabs( lDegX );
            const double lAbsDegY = fabs( lDegY );
            const double lAbsDegZ = fabs( lDegZ );
            const Vector3D & lTr = lArun.GetTranslation();

// DEBUG
//LzLogMsg("", "X deg.= "+lAbsDegX+", Y deg.="+lAbsDegY+", Z deg.= "+lAbsDegZ+", Tr= "+lTr.ToString()+".");

            if( lAbsDegX < pMaxEulerRotDeg
             && lAbsDegY < pMaxEulerRotDeg
             && lAbsDegZ < pMaxEulerRotDeg
             && lTr.Norm() < pMaxTransNorm )
            {
                // Log
                LzLogMsg("",  "ICP converged in " << iStep << " iteration(s)." );

                return lICP;
            }
        }

        // Next step
        iStep++;
    }

    // Error
    LzLogException("",  "ICP failed to converge after " << iStep << " steps!" );
}

//================================================================================
Matrix AffinePCAFitPairedPoints( vector<Point3D> & pSrc, const vector<Point3D> & pDst, double pTolToIdentity )
{
	// Create Vectors
    Vector<Point3D> lSrc( (size_t)pSrc.size() );
    for( size_t i=0 ; i<pSrc.size() ; i++ )
		lSrc[i] = pSrc[i];
	//
    Vector<Point3D> lDst( (size_t)pDst.size() );
    for( size_t i=0 ; i<pDst.size() ; i++ )
		lDst[i] = pDst[i];

	// Just do it already!
	return AffinePCAFitPairedPoints( lSrc, lDst, pTolToIdentity );
}

//================================================================================
Matrix AffinePCAFitPairedPoints( Vector<Point3D> & pSrc, const Vector<Point3D> & pDst, double pTolToIdentity )
{
static const size_t sMAX_ITER_COUNT = 100;

    // Check
    if( pSrc.Size() != pDst.Size() )
        LzLogException("",  "Points count mismatch! " << pSrc.Size() << " != " << pDst.Size() << "." );

    // Check
    if( pSrc.Size() < 3 )
        LzLogException("",  "Cannot register sets containing only " << pSrc.Size() << " point(s)! (minimum= 3 points)" );

	// Precompute PCA axis as they are only influenced by rigid rotation and not PCA scaling
    Point3D lG;
    Vector3D lV[3];
	//
	PCA( pSrc, lG, lV );

    // Initialize matrix
    Matrix lRes( 4, 4 );
    lRes.LoadIdentity();

    // Loop iterations
    bool lRigidAndAffineIsIdentity;
    size_t lIterCount = 0;
    do
    {
        // Affine iteration
        RigidTr3D lArun;
        Matrix lAffine;
		lRigidAndAffineIsIdentity = AffinePCAFitPairedPoints_Iteration( pSrc, pDst, lArun, lAffine, pTolToIdentity, lG, lV );

        // Assemble matrix pipeline
        Matrix lAffArun;
        {
            // Arun
            Matrix lM;
            lArun.ToMatrix4x4( lM );

            // Affine
            lAffine.Mult( lM, lAffArun );
        }

        // Accumulate and commit
        Matrix lNewRes;
        lAffArun.Mult( lRes, lNewRes );
        lRes = lNewRes;

        // Stop criterion
        lIterCount++;
    }
    while( lIterCount<sMAX_ITER_COUNT && !lRigidAndAffineIsIdentity );

    // Log
    LzLogMsg("",  "Finished after " << lIterCount << " iteration(s)." );

    // Finish
    return lRes;
}

//================================================================================
//bool AffinePCAFitPairedPoints_Iteration( Vector<Point3D> & pSrc, const Vector<Point3D> & pDst, RigidTr3D & pArun, Matrix & pAffine, double pTolToIdentity )
//{
//    //-----------------------
//    // Step 1: Arun
//    //-----------------------
//
//    // Compute transform
//    pArun = Arun( pSrc, pDst );
//
//    // Apply transform
//    for( size_t i = 0 ; i < pSrc.Size() ; i++ )
//        pSrc[i] *= pArun;
//
//    //-----------------------
//    // Step 2: Affine PCA
//    //-----------------------
//
//    // PCA on new source points
//    Vector3D lV[3];
//    Point3D lG;
//    PCA( pSrc, lG, lV ); // Achtung: scaling along principal axes does not affect the principal vectors. Principle vectors are only affected by the rigid rotation (Arun).
//
//    // Compute optimale scale factors
//    Matrix lX( 3 );
//    Matrix lMi[3];
//    {
//        // Tensor products
//        for( int i=0 ; i<3 ; i++ )
//            lV[i].TensorProd( lV[i], lMi[i] );
//
//        // Build linear system
//        Matrix A( 3, 3 ), B( 3 );
//        A.LoadValue( 0 );
//        B.LoadValue( 0 );
//
//        for( size_t p=0 ; p<pSrc.Size() ; p++ )
//        {
//            Vector3D GPp = pSrc[p] - lG;
//            Vector3D QpG = lG - pDst[p];
//
//            for( int i = 0 ; i < 3 ; i++ )
//            {
//                // Update B
//                B.Elt( i ) -= ( lMi[i] * GPp ) * QpG;
//
//                // Update A
//                for( int j = 0 ; j < 3 ; j++ )
//                    A.Elt( i, j ) += ( lMi[i] * GPp ) * ( lMi[j] * GPp );
//            }
//        }
//
//
//#pragma region "Solve linear system"
//        try
//        {
//            Matrix A_inv;
//            A.M3x3_InverseTo( A_inv );
//            A_inv.Mult( B, lX );
//        }
//        catch( ... )
//        {
//            // Log
//            LzLogMsg("",  "A is not inversible!" );
//            A.Log( "A" );
//
//            // Rank = 2?
//            Matrix A2x2( 2, 2 ), B2( 2 );
//            A.SubMatrixTo( 0, 0, A2x2 );
//            B.SubMatrixTo( 0, 0, B2 );
//
//            try
//            {
//                Matrix A2x2_inv, X2( 2 );
//                A2x2.M2x2_InverseTo( A2x2_inv );
//                A2x2_inv.Mult( B2, X2 );
//
//                // Copy scale factors
//                lX.Elt( 0 ) = X2.Elt( 0 );
//                lX.Elt( 1 ) = X2.Elt( 1 );
//                lX.Elt( 2 ) = 1;
//            }
//            catch( ... )
//            {
//                // Log
//                LzLogMsg("",  "A2x2 is not inversible!" );
//                A2x2.Log( "A2x2" );
//
//                // Rank = 1?
//                double A1x1 = A.Elt( 0, 0 );
//                double B1x1 = B.Elt( 0 );
//
//                try
//                {
//                    double X3;
//                    X3 = B1x1 / A1x1;
//
//                    // Copy scale factors
//                    lX.Elt( 0 ) = X3;
//                    lX.Elt( 1 ) = 1;
//                    lX.Elt( 2 ) = 1;
//                }
//                catch( ... )
//                {
//                    LzLogException("",  "Nope...." );
//                }
//            }
//        }
//#pragma endregion
//
//
////// Log solution
////lX.Log( "X: PCA affine scale factors" );
//    }
//
//    // Compute transform
//    pAffine.SetDims( 4, 4 );
//    {
//        // Assemble sum{ Xi Mi }
//        Matrix lXiMi( 3, 3 );
//        lXiMi.LoadValue( 0 );
//        for( int i = 0 ; i < 3 ; i++ )
//            lXiMi.AddMult( lX.Elt( i ), lMi[i], lXiMi );
//
//        // Square mat
//        pAffine.SubMatrixFrom( 0, 0, lXiMi );
//
//        // Last col
//        Point3D lXiMiG = lXiMi * lG;
//        for( int i = 0 ; i < 3 ; i++ )
//            pAffine.Elt( i, 3 ) = lG.mV[i] - lXiMiG.mV[i];
//
//        // Last row
//        pAffine.Elt( 3, 0 ) = pAffine.Elt( 3, 1 ) = pAffine.Elt( 3, 2 ) = 0;
//        pAffine.Elt( 3, 3 ) = 1;
//    }
//
//    // Apply transform
//    for( size_t i = 0 ; i < pSrc.Size() ; i++ )
//        pSrc[i] *= pAffine;
//
//    // Returns true if the affine fit is identity
//    return pArun.IsIdentity( pTolToIdentity )
//		&& fabs(lX.Elt(0) - 1)<=pTolToIdentity
//		&& fabs(lX.Elt(1) - 1)<=pTolToIdentity
//		&& fabs(lX.Elt(2) - 1)<=pTolToIdentity;
//}

//================================================================================
bool AffinePCAFitPairedPoints_Iteration( Vector<Point3D> & pSrc, const Vector<Point3D> & pDst, RigidTr3D & pArun, Matrix & pAffine, double pTolToIdentity,
                                         Point3D & pSrcG, Vector3D pSrcV[3], List<size_t> * ppLockedScales/*=nullptr*/ )
{
//// Log
//{
//    LzLogN("", "*** OLD INPUTS: "<<pSrc.Size()<<" --> "<<pDst.Size()<<"")
//    for( int i=0 ; i<10 ; i++ )
//        LzLogM("", pSrc[i].ToString()<<" --> "<<pDst[i].ToString())
//}

//// Log
//{
//    LzLogN("", "*** OLD PCA")
//    LzLogM("", "G= "<<pSrcG.ToString())
//    for( int i=0 ; i<3 ; i++ )
//        LzLogM("", "V= "<<pSrcV[i].ToString())
//}


    //-----------------------
    // Step 1: Arun
    //-----------------------

    // Compute transform
    pArun = Arun( pSrc, pDst );

    // Apply transform
    for( size_t i = 0 ; i < pSrc.Size() ; i++ )
        pSrc[i] *= pArun;

    // Update mass center
    pSrcG *= pArun;
	
	// Update eigenvectors
	for( int i=0 ; i<3 ; i++ )
		pSrcV[i] *= pArun;

    //-----------------------
    // Step 2: PCA scale
    //-----------------------

	Matrix lX( 3 );
    Matrix lMi[3];
    {
        // Tensor products
        for( int i=0 ; i<3 ; i++ )
            pSrcV[i].TensorProd( pSrcV[i], lMi[i] );

        // Build linear system
        Matrix A( 3, 3 ), B( 3 );
        A.LoadValue( 0 );
        B.LoadValue( 0 );

        for( size_t p=0 ; p<pSrc.Size() ; p++ )
        {
            Vector3D GPp = pSrc[p] - pSrcG;
            Vector3D QpG = pSrcG - pDst[p];

            for( int i=0 ; i<3 ; i++ )
            {
                // Update B
                B.Elt( i ) -= ( lMi[i] * GPp ) * QpG;

                // Update A
                for( int j = 0 ; j < 3 ; j++ )
                    A.Elt( i, j ) += ( lMi[i] * GPp ) * ( lMi[j] * GPp );
            }
        }

//// Log
//{
//    LzLogN("", "*** OLD SYSTEM")
//    A.Log("A");
//    B.Log("B");
//}


#pragma region "Solve linear system"
        try
        {
            Matrix A_inv;
            A.M3x3_InverseTo( A_inv );
            A_inv.Mult( B, lX );
        }
        catch( ... )
        {
            // Log
            LzLogMsg("",  "A is not inversible!" );
            A.Log( "A" );

            // Rank = 2?
            Matrix A2x2( 2, 2 ), B2( 2 );
            A.SubMatrixTo( 0, 0, A2x2 );
            B.SubMatrixTo( 0, 0, B2 );

            try
            {
                Matrix A2x2_inv, X2( 2 );
                A2x2.M2x2_InverseTo( A2x2_inv );
                A2x2_inv.Mult( B2, X2 );

                // Copy scale factors
                lX.Elt( 0 ) = X2.Elt( 0 );
                lX.Elt( 1 ) = X2.Elt( 1 );
                lX.Elt( 2 ) = 1;
            }
            catch( ... )
            {
                // Log
                LzLogMsg("",  "A2x2 is not inversible!" );
                A2x2.Log( "A2x2" );

                // Rank = 1?
                double A1x1 = A.Elt( 0, 0 );
                double B1x1 = B.Elt( 0 );

                try
                {
                    double X3;
                    X3 = B1x1 / A1x1;

                    // Copy scale factors
                    lX.Elt( 0 ) = X3;
                    lX.Elt( 1 ) = 1;
                    lX.Elt( 2 ) = 1;
                }
                catch( ... )
                {
                    LzLogException("",  "Nope...." );
                }
            }
        }
#pragma endregion


//// Log solution
//lX.Log( "X: PCA affine scale factors" );
    }

	// Lock scale?
	if( ppLockedScales )
	{
		BrowseList( iL, *ppLockedScales )
		{
			// Get dim
            size_t lDim = ppLockedScales->GetAt(iL);

			// Check and lock
			if( lDim < 3 )
				lX.Elt( lDim ) = 1.0; 
//*************** plutot que 1 : ==> valeur input imposee par caller
//*************** plutot que 1 : ==> valeur input imposee par caller
//*************** plutot que 1 : ==> valeur input imposee par caller
//*************** plutot que 1 : ==> valeur input imposee par caller
//*************** plutot que 1 : ==> valeur input imposee par caller
//*************** plutot que 1 : ==> valeur input imposee par caller

			else
                LzLogException("", "Unexpected dimension in locking list! Dim= "<<lDim<<".");
		}
	}

    //-----------------------
    // Compute transform
    //-----------------------

    pAffine.SetDims( 4, 4 );
    {
        // Assemble sum{ Xi Mi }
        Matrix lXiMi( 3, 3 );
        lXiMi.LoadValue( 0 );
        for( int i = 0 ; i < 3 ; i++ )
            lXiMi.AddMult( lX.Elt( i ), lMi[i], lXiMi );

        // Square mat
        pAffine.SubMatrixFrom( 0, 0, lXiMi );

        // Last col
        Point3D lXiMiG = lXiMi * pSrcG;
        for( int i = 0 ; i < 3 ; i++ )
            pAffine.Elt( i, 3 ) = pSrcG.mV[i] - lXiMiG.mV[i];

        // Last row
        pAffine.Elt( 3, 0 ) = pAffine.Elt( 3, 1 ) = pAffine.Elt( 3, 2 ) = 0;
        pAffine.Elt( 3, 3 ) = 1;
    }

    // Apply transform
    for( size_t i = 0 ; i < pSrc.Size() ; i++ )
        pSrc[i] *= pAffine;

    // Returns true if the affine fit is identity
    return pArun.IsIdentity( pTolToIdentity )
		&& fabs(lX.Elt(0) - 1)<=pTolToIdentity
		&& fabs(lX.Elt(1) - 1)<=pTolToIdentity
		&& fabs(lX.Elt(2) - 1)<=pTolToIdentity;
}
#pragma endregion


#pragma region "Pivot calibration"
//================================================================================
static double T_A2[3], T_B2[3], T_C2[3];
static double T_AB[3], T_AC[3], T_BC[3];
static double T_A [3], T_B [3], T_C [3], T_D [3];
static double S_A [3], S_B [3], S_C [3], S_D [3];
// Arrays ( 0=X , 1=Y , 2=Z )
//static void FillArrays( const CContainer1<CTransform3D> & p_TrRefA2RefB )
static void FillArrays( const Vector<RigidTr3D> & pTrRefA2RefB )
{
    // Reset
    for( int j = 0 ; j < 3 ; j++ )
    {
        T_A2[j] = T_B2[j] = T_C2[j] = 0;
        T_AB[j] = T_AC[j] = T_BC[j] = 0;
        T_A [j] = T_B [j] = T_C [j] = T_D [j] = 0;
        S_A [j] = S_B [j] = S_C [j] = S_D [j] = 0;
    }

    // Fill-in
    for( size_t i = 0 ; i < pTrRefA2RefB.Size() ; i++ )
    {
        // Access : M[row][col]
        Matrix M;
        /*const CMatrix4 & M = */pTrRefA2RefB[i].ToMatrix4x4( M ); //*************.GetMatrix4();

        for( int j = 0 ; j < 3 ; j++ )
        {
            T_A2[j] += M.Elt( j, 0 ) * M.Elt( j, 0 ); //M[j][0]*M[j][0];
            T_B2[j] += M.Elt( j, 1 ) * M.Elt( j, 1 ); //M[j][1]*M[j][1];
            T_C2[j] += M.Elt( j, 2 ) * M.Elt( j, 2 ); //M[j][2]*M[j][2];

            T_AB[j] += 2 * M.Elt( j, 0 ) * M.Elt( j, 1 ); //2*M[j][0]*M[j][1];
            T_AC[j] += 2 * M.Elt( j, 0 ) * M.Elt( j, 2 ); //2*M[j][0]*M[j][2];
            T_BC[j] += 2 * M.Elt( j, 1 ) * M.Elt( j, 2 ); //2*M[j][1]*M[j][2];

            T_A [j] += 2 * M.Elt( j, 0 ) * M.Elt( j, 3 ); //2*M[j][0]*M[j][3];
            T_B [j] += 2 * M.Elt( j, 1 ) * M.Elt( j, 3 ); //2*M[j][1]*M[j][3];
            T_C [j] += 2 * M.Elt( j, 2 ) * M.Elt( j, 3 ); //2*M[j][2]*M[j][3];
            T_D [j] +=   M.Elt( j, 3 ) * M.Elt( j, 3 ); //  M[j][3]*M[j][3];

            S_A [j] += M.Elt( j, 0 ); //M[j][0];
            S_B [j] += M.Elt( j, 1 ); //M[j][1];
            S_C [j] += M.Elt( j, 2 ); //M[j][2];
            S_D [j] += M.Elt( j, 3 ); //M[j][3];
        }
    }
}

//================================================================================
static void FillVariances( double n, double alf, double bet, double gam, double pVar[3] )
{
    // Compute variances
    for( int i = 0 ; i < 3 ; i++ )
    {
        double E_X2_, _EX_2;
        E_X2_ = ( alf * alf * T_A2[i] + alf * T_A[i] + alf * bet * T_AB[i] + alf * gam * T_AC[i]
                  + bet * bet * T_B2[i] + bet * T_B[i] + bet * gam * T_BC[i]
                  + gam * gam * T_C2[i] + gam * T_C[i]
                  + T_D[i]
                ) / n;

        _EX_2 = ( alf * alf * S_A[i] * S_A[i] + 2 * alf * S_A[i] * S_D[i] + 2 * alf * bet * S_A[i] * S_B[i] + 2 * alf * gam * S_A[i] * S_C[i]
                  + bet * bet * S_B[i] * S_B[i] + 2 * bet * S_B[i] * S_D[i] + 2 * bet * gam * S_B[i] * S_C[i]
                  + gam * gam * S_C[i] * S_C[i] + 2 * gam * S_C[i] * S_D[i]
                  + S_D[i] * S_D[i]
                ) / ( n * n ); // = p_PivotInRefB[i]^2

        pVar[i] = E_X2_ - _EX_2;
    }
}

//================================================================================
void BuildPivot6D( const Vector<RigidTr3D> & pTrRefA2RefB, Point3D & pPivotInRefA, Point3D & pPivotInRefB, double pVar[3] )
{

//**************************************************************************
//**************************************************************************
//
// Implementer un estimateur sans biais !!!!! /(n-1) .... a checker !
//
//**************************************************************************
//**************************************************************************



    // At least 2 elements ?
    //double n = pTrRefA2RefB.Size();
    const size_t n = pTrRefA2RefB.Size();
    if( n < 3 )
        LzLogException("", "Not enough transforms to compute pivot point!")

    // Fill-in
    FillArrays( pTrRefA2RefB );

    // Build matrix and vector
//CMatrix3 Mat( 0,0,0 , 0,0,0 , 0,0,0 );
//CCoord3 Vec( 0,0,0 );
    Matrix Mat( 3, 3 );
    Mat.LoadValue( 0 );

    Matrix Vec( 3 );
    Vec.LoadValue( 0 );

    for( int j = 0 ; j < 3 ; j++ )
    {
        Mat.Elt( 0, 0 )/*[0][0]*/ += 2 * ( T_A2[j] - S_A[j] * S_A[j] / n );
        Mat.Elt( 0, 1 )/*[0][1]*/ += T_AB[j] - 2 * S_A[j] * S_B[j] / n;
        Mat.Elt( 0, 2 )/*[0][2]*/ += T_AC[j] - 2 * S_A[j] * S_C[j] / n;

        Mat.Elt( 1, 0 )/*[1][0]*/ += T_AB[j] - 2 * S_A[j] * S_B[j] / n;
        Mat.Elt( 1, 1 )/*[1][1]*/ += 2 * ( T_B2[j] - S_B[j] * S_B[j] / n );
        Mat.Elt( 1, 2 )/*[1][2]*/ += T_BC[j] - 2 * S_B[j] * S_C[j] / n;

        Mat.Elt( 2, 0 )/*[2][0]*/ += T_AC[j] - 2 * S_A[j] * S_C[j] / n;
        Mat.Elt( 2, 1 )/*[2][1]*/ += T_BC[j] - 2 * S_B[j] * S_C[j] / n;
        Mat.Elt( 2, 2 )/*[2][2]*/ += 2 * ( T_C2[j] - S_C[j] * S_C[j] / n );

        Vec.Elt( 0 )/*[0]*/ -= T_A[j] - 2 * S_A[j] * S_D[j] / n;
        Vec.Elt( 1 )/*[1]*/ -= T_B[j] - 2 * S_B[j] * S_D[j] / n;
        Vec.Elt( 2 )/*[2]*/ -= T_C[j] - 2 * S_C[j] * S_D[j] / n;
    }

    // Log
    {
        LzLogNode("",  "Solving:")
//CString S[] = { "X", "Y", "Z" };
        Mat.Log( "M" );
        Vec.Log( "V" );
//char S[] = { 'X', 'Y', 'Z' };
//          for( int i=0 ; i<3 ; i++ )
//          {
//              LzLogM("", "| %.3f\t%.3f\t%.3f\t| | %s | %s | %.3f\t|",
//                  Mat[i][0],Mat[i][1],Mat[i][2],S[i],i==1?"=":" ",Vec[i]);
//          }
    }

    // Solve
//CMatrix3 InvMat;
//if( !Mat.GetInverse(InvMat) )
//{
//  LOGERR("[CGraphTools::BuildPivot6D] : Could not solve system !");
//  return FALSE;
//}
    Matrix InvMat;
    try
    {
        Mat.M3x3_InverseTo( InvMat );
    }
    catch( ... )
    {
        LzLogException("",  "Could not solve system! Cannot invert matrix." );
    }


//CCoord3 Sol = InvMat * Vec;
    Matrix Sol;
    InvMat.Mult( Vec, Sol );

    double alf, bet, gam;
    pPivotInRefA.X() = alf = Sol.Elt( 0 ); //[0];
    pPivotInRefA.Y() = bet = Sol.Elt( 1 ); //[1];
    pPivotInRefA.Z() = gam = Sol.Elt( 2 ); //[2];


//********************************************************************************* not necessarily correct !!!!
//********************************************************************************* not necessarily correct !!!!
    // Transfer solution to ref B
    pPivotInRefB = Point3D( 0, 0, 0 );
    for( size_t i = 0 ; i < n ; i++ )
    {
        Point3D lPivInB = pTrRefA2RefB[i] * pPivotInRefA;

        pPivotInRefB.X() += lPivInB.X();
        pPivotInRefB.Y() += lPivInB.Y();
        pPivotInRefB.Z() += lPivInB.Z();
//      p_PivotInRefB.X() += (p_TrRefA2RefB[i] * p_PivotInRefA).X();
//      p_PivotInRefB.Y() += (p_TrRefA2RefB[i] * p_PivotInRefA).Y();
//      p_PivotInRefB.Z() += (p_TrRefA2RefB[i] * p_PivotInRefA).Z();
    }
    pPivotInRefB.X() /= n;
    pPivotInRefB.Y() /= n;
    pPivotInRefB.Z() /= n;
//********************************************************************************* not necessarily correct !!!!
//********************************************************************************* not necessarily correct !!!!

    // Compute variances
    FillVariances( n, alf, bet, gam, pVar );

    // Log
    LzLogMsg("",  "Solution in A = " << pPivotInRefA.ToString() );
    LzLogMsg("",  "Solution in B = " << pPivotInRefB.ToString() );
    LzLogMsg("",  "Variances = ( " << pVar[0] << ", " << pVar[1] << ", " << pVar[2] << " )" );

    // Log PCA bounding box
    {
        Vector<Point3D> lPivotsInB( n );
        for( size_t i = 0 ; i < n ; i++ )
            lPivotsInB[i] = pTrRefA2RefB[i] * pPivotInRefA;

        // Compute PCA
        Point3D lO;
        Vector3D lV[3];
        PCA( lPivotsInB, lO, lV );

        // Compute min/max
        Vector<double> lMin, lMax;
        PCAMinMax( lPivotsInB, lO, lV, lMin, lMax );

        // Log
        LzLogN("",  "Min/Max PCA in Ref B:" );
        LzLogM("",  "Span X: " << ( lMax[0] - lMin[0] ) );
        LzLogM("",  "Span Y: " << ( lMax[1] - lMin[1] ) );
        LzLogM("",  "Span Z: " << ( lMax[2] - lMin[2] ) );
    }

    ////            // Proba loi normale uniformisee: P(-u<X<u)=90% => u=1.65 / 95% => u=1.96 / 99% => u=3.0
    ////            // http://www.statsoft.com/textbook/sttable.html
    ////            double u = 1.65;
    ////
    ////            // Trust interval sizes
    ////            double TX = 2*u*sqrt(Var[0]);
    ////            double TY = 2*u*sqrt(Var[1]);
    ////            double TZ = 2*u*sqrt(Var[2]);
}

//================================================================================
//BOOL CGraphTools::BuildPivot5D( const CContainer1<CTransform3D> & p_TrRefA2RefB,
//                              CPoint3D & p_PivotInRefA,
//                              CPoint3D & p_PivotInRefB,
//                              double p_Var[3],
//                              EPivotAxis p_Axis/*=PIVOT_Z*/ )
//{
//  // At least 2 elements ?
//  double n = p_TrRefA2RefB.Width();
//  if( n < 2 )
//  {
//      LOGERR("[CGraphTools::BuildPivot5D] : Not enough transforms to compute pivot point !");
//      return FALSE;
//  }

//  // Fill-in
//  FillArrays( p_TrRefA2RefB );

//  // Compute solution = numerator/denominator
//  double Num = 0;
//  double Den = 0;
//  int j;
//  switch( p_Axis )
//  {
//  // Solution in X
//  case PIVOT_X:
//      for( j=0 ; j<3 ; j++ )
//      {
//          Num -= T_A[j] - 2*S_A[j]*S_D[j]/n;
//          Den += 2*( T_A2[j] - S_A[j]*S_A[j]/n );
//      }
//      break;
//  // Solution in Y
//  case PIVOT_Y:
//      for( j=0 ; j<3 ; j++ )
//      {
//          Num -= T_B[j] - 2*S_B[j]*S_D[j]/n;
//          Den += 2*( T_B2[j] - S_B[j]*S_B[j]/n );
//      }
//      break;
//  // Solution in Z
//  case PIVOT_Z:
//      for( j=0 ; j<3 ; j++ )
//      {
//          Num -= T_C[j] - 2*S_C[j]*S_D[j]/n;
//          Den += 2*( T_C2[j] - S_C[j]*S_C[j]/n );
//      }
//      break;
//  // Unknown
//  default:
//      LOGERR("[CGraphTools::BuildPivot5D] : Unknown axis !");
//      return FALSE;
//      break;
//  }

//  // Solvable ?
//  if( Den == 0 )
//  {
//      LOGERR("[CGraphTools::BuildPivot5D] : Could not solve system !");
//      return FALSE;
//  }

//  // Compute solution
//  double alf = 0;
//  double bet = 0;
//  double gam = 0;
//  p_PivotInRefA = CPoint3D(0,0,0);
//  switch( p_Axis )
//  {
//  // Solution in X
//  case PIVOT_X:
//      p_PivotInRefA.X() = alf = Num/Den;
//      break;
//  // Solution in Y
//  case PIVOT_Y:
//      p_PivotInRefA.Y() = bet = Num/Den;
//      break;
//  // Solution in Z
//  case PIVOT_Z:
//      p_PivotInRefA.Z() = gam = Num/Den;
//      break;
//  }

//  // Transfer solution to ref B
//  p_PivotInRefB = CPoint3D(0,0,0);
//  for( int i=0 ; i<n ; i++ )
//  {
//      p_PivotInRefB.X() += (p_TrRefA2RefB[i] * p_PivotInRefA).X();
//      p_PivotInRefB.Y() += (p_TrRefA2RefB[i] * p_PivotInRefA).Y();
//      p_PivotInRefB.Z() += (p_TrRefA2RefB[i] * p_PivotInRefA).Z();
//  }
//  p_PivotInRefB.X() /= n;
//  p_PivotInRefB.Y() /= n;
//  p_PivotInRefB.Z() /= n;

//  // Compute variances
//  FillVariances( n, alf, bet, gam, p_Var );

////    // Compute variances
////    for( i=0 ; i<3 ; i++ )
////    {
////        double E_X2_, _EX_2;
////        E_X2_ = ( alf*alf*T_A2[i] + alf*T_A[i] + alf*bet*T_AB[i] + alf*gam*T_AC[i]
////                + bet*bet*T_B2[i] + bet*T_B[i] + bet*gam*T_BC[i]
////                + gam*gam*T_C2[i] + gam*T_C[i]
////                + T_D[i]
////                ) / n;
////
////        _EX_2 = ( alf*alf*S_A[i]*S_A[i] + 2*alf*S_A[i]*S_D[i] + 2*alf*bet*S_A[i]*S_B[i] + 2*alf*gam*S_A[i]*S_C[i]
////                + bet*bet*S_B[i]*S_B[i] + 2*bet*S_B[i]*S_D[i] + 2*bet*gam*S_B[i]*S_C[i]
////                + gam*gam*S_C[i]*S_C[i] + 2*gam*S_C[i]*S_D[i]
////                + S_D[i]*S_D[i]
////                ) / (n*n); // = p_PivotInRefB[i]^2
////
////        p_Var[i] = E_X2_ - _EX_2;
////    }

//  // Log
//  LOGMSG("[CGraphTools::BuildPivot5D] : Solution in A = %s",p_PivotInRefA.AsString());
//  LOGMSG("[CGraphTools::BuildPivot5D] : Solution in B = %s",p_PivotInRefB.AsString());
//  LOGMSG("[CGraphTools::BuildPivot5D] : Variances = ( %f, %f, %f )",p_Var[0],p_Var[1],p_Var[2]);


////            // Proba loi normale uniformisee: P(-u<X<u)=90% => u=1.65 / 95% => u=1.96 / 99% => u=3.0
////            // http://www.statsoft.com/textbook/sttable.html
////            double u = 1.65;
////
////            // Trust interval sizes
////            double TX = 2*u*sqrt(Var[0]);
////            double TY = 2*u*sqrt(Var[1]);
////            double TZ = 2*u*sqrt(Var[2]);



//  return TRUE;
//}
#pragma endregion

}
}
