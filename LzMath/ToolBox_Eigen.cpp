#include "ToolBox.h"
#include "Matrix.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <LzServices/LzLog.h>
#include <random>


namespace LzMath {
namespace ToolBox
{
using LzGeom::BBox;


//================================================================================
void FitCartesianEllipsoid( const vector<Point3D> & pPts, Point3D & pCenter, double & pRadX, double & pRadY, double & pRadZ )
{
static const int sMaxSteps = 100;

// Log
LzLogTimeN("", "FitCartesianEllipsoid")

    // Compute BBox
    BBox lBB;
    lBB.Update( pPts );

    // BBox center
    const Point3D & lCenterBB = lBB.Center();

    // Initial guess
    Eigen::VectorXd X( 6 );
    X(0) = lCenterBB.X();
    X(1) = lCenterBB.Y();
    X(2) = lCenterBB.Z();
    X(3) = 0.5*lBB.Size(0);
    X(4) = 0.5*lBB.Size(1);
    X(5) = 0.5*lBB.Size(2);

    //----------------------------
    // Opti loop
    //----------------------------

    for( int iStep=0 ; iStep<sMaxSteps ; iStep++ )
    {
LzLogM("", "Step "<<iStep<<", X=\n"<<X)

            // System at current step: JJt * dX = b
            Eigen::MatrixXd JJt( 6, 6 );
            Eigen::VectorXd b( 6 );

            // Reset
            JJt.setZero();
            b.setZero();

            // Assemble system
            for( const Point3D & iPt : pPts )
            {
                const double DX = X(0) - iPt.X();
                const double DY = X(1) - iPt.Y();
                const double DZ = X(2) - iPt.Z();
                //
                const double DX_2 = DX * DX;
                const double DY_2 = DY * DY;
                const double DZ_2 = DZ * DZ;
                //
                const double RX = X(3);
                const double RY = X(4);
                const double RZ = X(5);
                //
                const double RX_2 = RX * RX;
                const double RY_2 = RY * RY;
                const double RZ_2 = RZ * RZ;
                //
                const double d_i = RY_2 * RZ_2 * DX_2
                                 + RX_2 * RZ_2 * DY_2
                                 + RX_2 * RY_2 * DZ_2
                                 - RX_2 * RY_2 * RZ_2;

                // Compute local system
                Eigen::VectorXd J_i( 6 );
                J_i(0) = RY_2 * RZ_2 * DX;
                J_i(1) = RX_2 * RZ_2 * DY;
                J_i(2) = RX_2 * RY_2 * DZ;
                J_i(3) = RX * (RZ_2 * DY_2 + RY_2 * DZ_2 - RY_2 * RZ_2);
                J_i(4) = RY * (RZ_2 * DX_2 + RX_2 * DZ_2 - RX_2 * RZ_2);
                J_i(5) = RZ * (RY_2 * DX_2 + RX_2 * DY_2 - RX_2 * RY_2);
                //
                const Eigen::MatrixXd JJt_i = 4/*=2x2 factors missing from J_i above*/ * J_i * J_i.transpose();

                // Accumulate local system to global system
                JJt += JJt_i;
                b -= d_i * J_i;
            }

            // Solve for increment
            Eigen::VectorXd dX = JJt.inverse() * b;
//            LzLogM("", "Step "<<iStep<<", dX=\n"<<dX)

            // Apply increment
            X += dX;

            // Stop criterion: compute infinite norm
            const double inf_norm_dX = dX.lpNorm<Eigen::Infinity>();
            if( inf_norm_dX < 1e-6 )
            {
LzLogM("", "Increment too small. Breaking loop...")
                break;
            }
    }

    // Unpack solution
    pCenter = { X(0), X(1), X(2) };
    pRadX = X(3);
    pRadY = X(4);
    pRadZ = X(5);
}

//================================================================================
void UnitTest_FitCartesianEllipsoid()
{
//    LZTriModel::Mesh lBB;
//    LzTriModel::GenEllipsoid( lBB, 1.0, 1.0, 1.0, 40, 80 );
//    //
//    MeshTopology lBBTopo;
//    lBBTopo.Set( lBB );

    std::uniform_real_distribution<double> lPosDist( -100.0, +100.0 );
    std::uniform_real_distribution<double> lRadDist(   +1.0, +100.0 );
    std::uniform_real_distribution<double> lNoiseDist( -5.0,   +5.0 );

//see D:\LinkZ\Tech\_Code\Twinsight\kPCA\kPCA\mainwindow.cpp:2555
}

//================================================================================
RigidTr3D Arun( const vector<Point3D> & pP1, const vector<Point3D> & pP2 )
{
    // Check
    const size_t N = pP1.size();
    if( N<3 || N!=pP2.size() )
        LzLogException("", "Arun registration impossible: "<<pP1.size()<<" points vs "<<pP2.size()<<" points!")

    // Compute centroids
    const Point3D lG1 = Centroid( pP1 );
    const Point3D lG2 = Centroid( pP2 );

    // Build matrix H
    Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
    for( size_t i=0 ; i<N ; i++ )
    {
        const Point3D Q1(pP1[i].X() - lG1.X(),
                         pP1[i].Y() - lG1.Y(),
                         pP1[i].Z() - lG1.Z());

        const Point3D Q2(pP2[i].X() - lG2.X(),
                         pP2[i].Y() - lG2.Y(),
                         pP2[i].Z() - lG2.Z());

        for( int i1=0 ; i1<3 ; i1++ )
        for( int i2=0 ; i2<3 ; i2++ )
            H( i1, i2 ) += Q1.mV[i1] * Q2.mV[i2];
    }

    // SVD: H = USVt
    //
    // https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
    // --> Singular values are always sorted in decreasing order
    //
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Compute X
    Eigen::Matrix3d X = svd.matrixV() * svd.matrixU().transpose();

    //===> Check 1
    const double lDetX = X.determinant();
    Eigen::Matrix3d V = svd.matrixV();
    if( lDetX < 0 )
    {
        // Look for the smallest singular value
        unsigned int lMinS;
#if 1
        // --> Singular values are always sorted in decreasing order
        lMinS = 2;
#else
        {
            const double lS0 = svd.singularValues()( 0 );
            const double lS1 = svd.singularValues()( 1 );
            const double lS2 = svd.singularValues()( 2 );

            if( lS1 < lS2 )
                lMinS = lS1 < lS0 ? 1 : 0 ;
            else
                lMinS = lS2 < lS0 ? 2 : 0 ;
        }
#endif
        // Invert the corresponding vector
        V( 0, lMinS ) *= -1.0;
        V( 1, lMinS ) *= -1.0;
        V( 2, lMinS ) *= -1.0;

        // Update X
        X = V * svd.matrixU().transpose();

        //===> Check 2
        const double lDetX = X.determinant();
        if( lDetX < 0 )
        {
            // Log singular values
            LzLogErr("", "Singular values\n"<<svd.singularValues());
            LzLogException("", "Arun registration failed! (|X|= "<<lDetX<<")")
        }
    }

    // Set rotation
    RigidTr3D lP1toP2;
    for( int i=0 ; i<3 ; i++ )
    for( int j=0 ; j<3 ; j++ )
        lP1toP2.mR[i * 3 + j] = X( i, j );

    // Compute and set translation
    lP1toP2.SetTranslation( lG2 - lP1toP2*lG1 );

    // Done
    return lP1toP2;

/*
    OLD IMPLEMENTATION

    // Check
    const size_t N = pP1.size();
    if( N<3 || N!=pP2.size() )
        LzLogException("", "Arun registration impossible: "<<pP1.size()<<" points vs "<<pP2.size()<<"points!")

    // Compute centroids
    Point3D lG1 = Centroid( pP1 );
    Point3D lG2 = Centroid( pP2 );

    // Shift points to centroid
    vector<Point3D> lQ1( N );
    vector<Point3D> lQ2( N );
    for( size_t i=0 ; i<N ; i++ )
    {
        // 1
        lQ1[i].X() = pP1[i].X() - lG1.X();
        lQ1[i].Y() = pP1[i].Y() - lG1.Y();
        lQ1[i].Z() = pP1[i].Z() - lG1.Z();

        // 2
        lQ2[i].X() = pP2[i].X() - lG2.X();
        lQ2[i].Y() = pP2[i].Y() - lG2.Y();
        lQ2[i].Z() = pP2[i].Z() - lG2.Z();
    }

    // Build matrix H
    Matrix H( 3, 3 );
    H.LoadValue( 0 );
    for( size_t i=0 ; i<N ; i++ )
    {
        const Point3D & Q1 = lQ1[i];
        const Point3D & Q2 = lQ2[i];

        for( int i1=0 ; i1<3 ; i1++ )
        for( int i2=0 ; i2<3 ; i2++ )
            H.Elt( i1, i2 ) += Q1.mV[i1] * Q2.mV[i2];
    }

    // Singular value decomposition
    Matrix U, W, V;
    SVD( H, U, W, V );



    // Compute X
    Matrix Ut, X;
    U.TransposeTo( Ut );
    V.Mult( Ut, X );

    //===> Check 1
    double lDetX = X.M3x3_Det();
    if( lDetX < 0 )
    {
        // Look for the smallest singular value
        unsigned int lMinS;
        {
            const double lS0 = W.Elt( 0 );
            const double lS1 = W.Elt( 1 );
            const double lS2 = W.Elt( 2 );
            if( lS1 < lS2 )
                lMinS = lS1 < lS0 ? 1 : 0 ;
            else
                lMinS = lS2 < lS0 ? 2 : 0 ;
        }

        // Invert the corresponding vector
        V.Elt( 0, lMinS ) *= -1.0;
        V.Elt( 1, lMinS ) *= -1.0;
        V.Elt( 2, lMinS ) *= -1.0;

        // Update X
        V.Mult( Ut, X );

        //===> Check 2
        double lDetX = X.M3x3_Det();
        if( lDetX < 0 )
        {
            // Log singular values
            W.Log( "W: singular values" );
            LzLogException("", "Arun registration failed! (|X|= "<<lDetX<<")")
        }
    }

    // Set rotation
    RigidTr3D lP1toP2;
    lP1toP2.FromMatrix3x3( X );

    // Compute and set translation
    lP1toP2.SetTranslation( lG2 - lP1toP2 * lG1 );

    // Done
    return lP1toP2;
*/
}

//================================================================================
void InertialPointCloudPCA( const vector<Point3D> & pCloud, const vector<double> & pWeights, Point3D & pMean, Vector3D pV[3] )
{
    // https://en.wikipedia.org/wiki/Moment_of_inertia

    const unsigned int lNbPoints = pCloud.size();

    // Check
    if( lNbPoints < 2 )
        LzLogException("", "Cannot compute inertial PCA on less than 2 points!" );

    // Check
    if( lNbPoints != pWeights.size() )
        LzLogException("", "Size mismatch between points and weights! ("<<lNbPoints<<" != "<<pWeights.size()<<")")


    //----------------------------
    // Compute centroid
    //----------------------------

    pMean.Reset();
    double lTotWeight = 0;
    //
    for( size_t i=0 ; i<pCloud.size() ; i++ )
    {
        const Point3D & lP = pCloud[i];
        const double lW = pWeights[i];

        lTotWeight += lW;

        pMean.X() += lW * lP.X();
        pMean.Y() += lW * lP.Y();
        pMean.Z() += lW * lP.Z();
    }

    // Check
    if( IsZero( lTotWeight ) )
        LzLogException("",  "Could not normalize weighted centroid by weight= " << lTotWeight << "!" );

    // Normalize
    pMean /= lTotWeight;

// Log
LzLogM("", "Mass center= "<<pMean.ToString())
LzLogM("", "Tot weight= "<<lTotWeight)


    //----------------------------
    // Compute inertial matrix
    //----------------------------

    // Set empty I
    Eigen::Matrix3d I;
    I.setZero();

    // Compute I
    //
    // I_ij = sum{k=1 .. N} [ M_k * (|r_k|^2 * dirac_ij - x_i^k * x_j^k) ]
    //
//******************************** precompute all points first (substract Mean.XYZ and compute Rk
//******************************** precompute all points first (substract Mean.XYZ and compute Rk
//******************************** precompute all points first (substract Mean.XYZ and compute Rk
//******************************** precompute all points first (substract Mean.XYZ and compute Rk
//******************************** precompute all points first (substract Mean.XYZ and compute Rk
    for( int i=0 ; i<3 ; i++ )
    for( int j=i ; j<3 ; j++ )
    {
        // Computing upper triangular: i<=j

        // Accumulate terms
        for( size_t k=0 ; k<pCloud.size() ; k++ )
        {
            // Position vector from Mean
            const Vector3D lRk = pCloud[k] - pMean; //******************** valider qu'il faut bien se mettre en Mean
//            const Vector3D lRk = pCloud[k] - Point3D(0,0,0);

            // Common term
            double lTerm = -lRk.mV[i] * lRk.mV[j];

            // Diagonal?
            if( i == j )
            {
                // Update term
                lTerm += lRk.SquaredNorm();
            }

            // Update tensor
            I(i, j) += pWeights[k] * lTerm;
        }
    }
    //
    // Add lower triangular
    I( 1, 0 ) = I( 0, 1 );
    I( 2, 0 ) = I( 0, 2 );
    I( 2, 1 ) = I( 1, 2 );

LzLogM("", "I=\n"<<I);

    // SVD: H = USVt
    //
    // https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
    // --> Singular values are always sorted in decreasing order
    //
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(I, Eigen::ComputeFullU);

    // Read eigenvectors
    const Eigen::Matrix3d & U = svd.matrixU();

    // Vectors produced by EigenJacobi are already normalized
    pV[0] = Vector3D( U( 0, 0 ), U( 1, 0 ), U( 2, 0 ) );
    pV[1] = Vector3D( U( 0, 1 ), U( 1, 1 ), U( 2, 1 ) );
    pV[2] = Vector3D( U( 0, 2 ), U( 1, 2 ), U( 2, 2 ) );
}

}
}
