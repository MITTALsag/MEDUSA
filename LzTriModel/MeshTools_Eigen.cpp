#include "MeshTools.h"
#include <LzMath/ToolBox.h>
#include <LzGeom/Line3D.h>
//
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <LzServices/LzLog.h>
//
#if 0
#include <LzMath/Polynom.h> // Used in TetraTensor_V2 and TEST_InertiaTensor
#endif
//
#include <iomanip>


namespace LzTriModel
{



////================================================================================
//RigidTr3D SurfacicPCA( const Mesh & pObject )
//{
//    // Get face centroids and areas
//    vector<Point3D> lCents;
//    vector<double> lAreas;
//    pObject.GetFacesCentroids( lCents, lAreas );

// Eigen weighted PCA ??? if not ==> LzMath::ToolBox::WeightedPCA( ... )

//    return RigidTr3D();
//}

//================================================================================
/* Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates
 * F. Tonon
 * Department of Geology and Geophysics, University of Utah
 * Journal of Mathematics and Statistics 1 (1): 8-11, 2004
 * ISSN 1549-3644
 *
 *     |  a  -b' -c' |
 * I = | -b'  b  -a' |
 *     | -c' -a'  c  |
 */
static void TetraTensor( Eigen::Matrix3d & pTo, const vector<Point3D> & pVers, double pMu )
{
    // Here: assuming pVers holds 4 vertices
    // =====

    // Here: assuming the vertices have been offset by the mass center's coordinates
    // =====

    // Adopt the paper's notations
    const Point3D & A1 = pVers[0];
    const Point3D & A2 = pVers[1];
    const Point3D & A3 = pVers[2];
    const Point3D & A4 = pVers[3];
    //
    const double x1 = A1.X();
    const double y1 = A1.Y();
    const double z1 = A1.Z();
    //
    const double x2 = A2.X();
    const double y2 = A2.Y();
    const double z2 = A2.Z();
    //
    const double x3 = A3.X();
    const double y3 = A3.Y();
    const double z3 = A3.Z();
    //
    const double x4 = A4.X();
    const double y4 = A4.Y();
    const double z4 = A4.Z();

    // Jacobian = 6 * Vol, possibly negative value
    const double DetJ = /*fabs(*/ ((A2 - A1)^(A3 - A1)) * (A4 - A1) /*)*/; // No fabs to account for inverted tets

//// Log
//LzLogM("", "DetJ= "<<DetJ)

    // Compute coefficients
    const double a =
        pMu * DetJ * ( y1*y1 + y1*y2 + y2*y2 + y1*y3 + y2*y3
                     + y3*y3 + y1*y4 + y2*y4 + y3*y4 + y4*y4 + z1*z1 + z1*z2
                     + z2*z2 + z1*z3 + z2*z3 + z3*z3 + z1*z4 + z2*z4 + z3*z4 + z4*z4 ) / 60;
    //
    const double b =
        pMu * DetJ * ( x1*x1 + x1*x2 + x2*x2 + x1*x3 + x2*x3 + x3*x3
                     + x1*x4 + x2*x4 + x3*x4 + x4*x4 + z1*z1 + z1*z2 + z2*z2 + z1*z3
                     + z2*z3 + z3*z3 + z1*z4 + z2*z4 + z3*z4 + z4*z4 ) / 60;
    //
    const double c =
        pMu * DetJ * ( x1*x1 + x1*x2 + x2*x2 + x1*x3 + x2*x3 + x3*x3 + x1*x4
                     + x2*x4 + x3*x4 + x4*x4 + y1*y1 + y1*y2 + y2*y2 + y1*y3
                     + y2*y3 + y3*y3 + y1*y4 + y2*y4 + y3*y4 + y4*y4 ) / 60;
    //
    const double a_p =
        pMu * DetJ * ( 2*y1*z1 + y2*z1 + y3*z1 + y4*z1 + y1*z2
                     + 2*y2*z2 + y3*z2 + y4*z2 + y1*z3 + y2*z3 + 2*y3*z3
                     + y4*z3 + y1*z4 + y2*z4 + y3*z4 + 2*y4*z4 ) / 120;
    //
    const double b_p =
        pMu * DetJ * ( 2*x1*z1 + x2*z1 + x3*z1 + x4*z1 + x1*z2
                     + 2*x2*z2 + x3*z2 + x4*z2 + x1*z3 + x2*z3 + 2*x3*z3
                     + x4*z3 + x1*z4 + x2*z4 + x3*z4 + 2*x4*z4 ) / 120;
    //
    const double c_p =
        pMu * DetJ * ( 2*x1*y1 + x2*y1 + x3*y1 + x4*y1 + x1*y2
                     + 2*x2*y2 + x3*y2 + x4*y2 + x1*y3 + x2*y3 + 2*x3*y3
                     + x4*y3 + x1*y4 + x2*y4 + x3*y4 + 2*x4*y4 ) / 120;

    // Set tensor values
    //
    // Diagonal
    pTo( 0, 0 ) = a;
    pTo( 1, 1 ) = b;
    pTo( 2, 2 ) = c;
    //
    // Off-diagonal
//    pTo( 1, 0 ) = pTo( 0, 1 ) = -b_p; //nope, this term relates x and z ==> 0 and 2
//    pTo( 2, 0 ) = pTo( 0, 2 ) = -c_p; //nope, this term relates x and y ==> 0 and 1
    pTo( 1, 0 ) = pTo( 0, 1 ) = -c_p;
    pTo( 2, 0 ) = pTo( 0, 2 ) = -b_p;
    //
    pTo( 2, 1 ) = pTo( 1, 2 ) = -a_p;
}

//================================================================================
// Solid Body Rotation and the Inertia Tensor
// https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Classical_Mechanics_(Tatum)/02%3A_Moments_of_Inertia/2.17%3A_Solid_Body_Rotation_and_the_Inertia_Tensor
//
// Why the inertia tensor is the inertia tensor.
// http://number-none.com/blow/inertia/deriving_i.html
//
//
static void TetraTensor_V2( Eigen::Matrix3d & pTo, const vector<Point3D> & pVers, double pMu )
{
    // Here: assuming pVers holds 4 vertices
    // =====

    // Adopt the paper's notations
    const Point3D & A1 = pVers[0];
    const Point3D & A2 = pVers[1];
    const Point3D & A3 = pVers[2];
    const Point3D & A4 = pVers[3];

    // Jacobian = 6 * Vol, possibly negative value
    const double DetJ = /*fabs(*/ ((A2 - A1)^(A3 - A1)) * (A4 - A1) /*)*/; // No fabs to account for inverted tets

#if 1
    // Symmetric terms
    //
    // Integral_over_the_Tet{ (A + B*_X_ + C*_Y_ + D*_Z_)^2 } dX dY dZ
    //
    #define Sym( A, B, C, D )\
        ( 10*A*A + 5*A*B + 5*A*C + 5*A*D +\
             B*B +   B*C +   B*D +   C*C +\
             C*D +   D*D ) / 60.0

    // Skew terms
    //
    // Integral_over_the_Tet{ (A + B*_X_ + C*_Y_ + D*_Z_)*
    //                        (E + F*_X_ + G*_Y_ + H*_Z_) } dX dY dZ
    //
    #define Skew( A, B, C, D, E, F, G, H )\
        ( 20*A*E + 5*A*F + 5*A*G + 5*A*H +\
           5*B*E + 2*B*F + B*G   + B*H +\
           5*C*E +   C*F + 2*C*G + C*H +\
           5*D*E +   D*F + D*G   + 2*D*H ) / 120.0

    // Precompute stuff
    const double Ax = A1.X();
    const double Bx = A2.X() - Ax;
    const double Cx = A3.X() - Ax;
    const double Dx = A4.X() - Ax;
    //
    const double Ay = A1.Y();
    const double By = A2.Y() - Ay;
    const double Cy = A3.Y() - Ay;
    const double Dy = A4.Y() - Ay;
    //
    const double Az = A1.Z();
    const double Bz = A2.Z() - Az;
    const double Cz = A3.Z() - Az;
    const double Dz = A4.Z() - Az;

    // Sym
    const double lInt_x_2 = Sym( Ax, Bx, Cx, Dx );
    const double lInt_y_2 = Sym( Ay, By, Cy, Dy );
    const double lInt_z_2 = Sym( Az, Bz, Cz, Dz );

    // Skew
    const double lInt_xy = Skew( Ax, Bx, Cx, Dx, Ay, By, Cy, Dy );
    const double lInt_xz = Skew( Ax, Bx, Cx, Dx, Az, Bz, Cz, Dz );
    const double lInt_yz = Skew( Ay, By, Cy, Dy, Az, Bz, Cz, Dz );
#else
    //
    const double x1 = A1.X();
    const double y1 = A1.Y();
    const double z1 = A1.Z();
    //
    const double x2 = A2.X();
    const double y2 = A2.Y();
    const double z2 = A2.Z();
    //
    const double x3 = A3.X();
    const double y3 = A3.Y();
    const double z3 = A3.Z();
    //
    const double x4 = A4.X();
    const double y4 = A4.Y();
    const double z4 = A4.Z();

    // Mapping (X,Y,Z) --> (x,y,z)
    const LzMath::Polynom x = x1 + _X_*(x2 - x1) + _Y_*(x3 - x1) + _Z_*(x4 - x1);
    const LzMath::Polynom y = y1 + _X_*(y2 - y1) + _Y_*(y3 - y1) + _Z_*(y4 - y1);
    const LzMath::Polynom z = z1 + _X_*(z2 - z1) + _Y_*(z3 - z1) + _Z_*(z4 - z1);

    // Integrands (elements of the covariance matrix)
    const LzMath::Polynom x_2 = x*x;
    const LzMath::Polynom y_2 = y*y;
    const LzMath::Polynom z_2 = z*z;
    //
    const LzMath::Polynom xy = x*y;
    const LzMath::Polynom xz = x*z;
    const LzMath::Polynom yz = y*z;

    // Integrate x_2
    double lInt_x_2;
    {
        LzMath::Polynom x_2_dZ = x_2.IntZ();
        x_2_dZ = x_2_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ ); // No need to subtract Compose( _X_, _Y_, 0.0 ) since it is 0

        LzMath::Polynom x_2_dZ_dY = x_2_dZ.IntY();
        x_2_dZ_dY = x_2_dZ_dY.Compose( _X_, 1.0-_X_, 0.0 ); // No need to subtract Compose( _X_, 0.0, 0.0 ) since it is 0

        LzMath::Polynom x_2_dZ_dY_dX = x_2_dZ_dY.IntX();
        lInt_x_2 = x_2_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 ); // No need to subtract EvalAt( 0.0, 0.0, 0.0 ) since it is 0
    }

    // Integrate y_2
    double lInt_y_2;
    {
        LzMath::Polynom y_2_dZ = y_2.IntZ();
        y_2_dZ = y_2_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ );

        LzMath::Polynom y_2_dZ_dY = y_2_dZ.IntY();
        y_2_dZ_dY = y_2_dZ_dY.Compose( _X_, 1-_X_, 0.0 );

        LzMath::Polynom y_2_dZ_dY_dX = y_2_dZ_dY.IntX();

        lInt_y_2 = y_2_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 );
    }

    // Integrate z_2
    double lInt_z_2;
    {
        LzMath::Polynom z_2_dZ = z_2.IntZ();
        z_2_dZ = z_2_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ );

        LzMath::Polynom z_2_dZ_dY = z_2_dZ.IntY();
        z_2_dZ_dY = z_2_dZ_dY.Compose( _X_, 1.0-_X_, 0.0 );

        LzMath::Polynom z_2_dZ_dY_dX = z_2_dZ_dY.IntX();

        lInt_z_2 = z_2_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 );
    }

    // Integrate xy
    double lInt_xy;
    {
        LzMath::Polynom xy_dZ = xy.IntZ();
        xy_dZ = xy_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ );

        LzMath::Polynom xy_dZ_dY = xy_dZ.IntY();
        xy_dZ_dY = xy_dZ_dY.Compose( _X_, 1.0-_X_, 0.0 );

        LzMath::Polynom xy_dZ_dY_dX = xy_dZ_dY.IntX();

        lInt_xy = xy_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 );
    }

    // Integrate xz
    double lInt_xz;
    {
        LzMath::Polynom xz_dZ = xz.IntZ();
        xz_dZ = xz_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ );

        LzMath::Polynom xz_dZ_dY = xz_dZ.IntY();
        xz_dZ_dY = xz_dZ_dY.Compose( _X_, 1.0-_X_, 0.0 );

        LzMath::Polynom xz_dZ_dY_dX = xz_dZ_dY.IntX();

        lInt_xz = xz_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 );
    }

    // Integrate yz
    double lInt_yz;
    {
        LzMath::Polynom yz_dZ = yz.IntZ();
        yz_dZ = yz_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ );

        LzMath::Polynom yz_dZ_dY = yz_dZ.IntY();
        yz_dZ_dY = yz_dZ_dY.Compose( _X_, 1.0-_X_, 0.0 );

        LzMath::Polynom yz_dZ_dY_dX = yz_dZ_dY.IntX();

        lInt_yz = yz_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 );
    }

#endif

    // Set tensor values
    //
    // Using C: covariance matrix
    //
    //     | x2 xy xz |
    // C = | xy y2 yz |
    //     | xz yz z2 |
    //
    // I = trace(C)*Id_3x3 - C
    //
    const double K = pMu * DetJ;
    //
    // Diagonal
    pTo( 0, 0 ) = K*( lInt_y_2 + lInt_z_2 );
    pTo( 1, 1 ) = K*( lInt_x_2 + lInt_z_2 );
    pTo( 2, 2 ) = K*( lInt_x_2 + lInt_y_2 );
    //
    // Off-diagonal
    pTo( 1, 0 ) = pTo( 0, 1 ) = -K*lInt_xy;
    pTo( 2, 0 ) = pTo( 0, 2 ) = -K*lInt_xz;
    pTo( 2, 1 ) = pTo( 1, 2 ) = -K*lInt_yz;
}

//================================================================================
void TEST_InertiaTensor( const Mesh & pObject, Point3D & pMean, Eigen::Matrix3d & pI )
{
#if 1
{
    {
        const int TESTS = 1000000;

        std::uniform_real_distribution<double> lDist( -200, +200 );

        vector<double> lRnd_A( TESTS );
        vector<RigidTr3D> lRnd_T( TESTS );
        for( int test=0 ; test<TESTS ; test++ )
        {
            lRnd_A[test] = lDist( LzServices::RandomEngine() );
            lRnd_T[test] = RigidTr3D( lDist( LzServices::RandomEngine() ),
                                      lDist( LzServices::RandomEngine() ),
                                      lDist( LzServices::RandomEngine() ),
                                      Vector3D(lDist( LzServices::RandomEngine() ),
                                               lDist( LzServices::RandomEngine() ),
                                               lDist( LzServices::RandomEngine() )) );
        }

        LzLogTimeN("", "Tensors")

        {
            LzLogTimeN("", "OLD")
            for( int test=0 ; test<TESTS ; test++ )
            {
                // Compute intertia tensor for tetrahedron
                const double a = lRnd_A[test];
                const Vector3D lCenter( a*0.25, a*0.25, a*0.25 );
                const vector<Point3D> lVers =
                {
                    lRnd_T[test]*(Point3D(0, 0, 0) - lCenter),
                    lRnd_T[test]*(Point3D(a, 0, 0) - lCenter),
                    lRnd_T[test]*(Point3D(0, a, 0) - lCenter),
                    lRnd_T[test]*(Point3D(0, 0, a) - lCenter)
                };

                Eigen::Matrix3d TetI;
                TetraTensor( TetI, lVers, 1.0 );

                if( test+1 == TESTS )
                    LzLogM("", "TetI=\n"<<TetI)
            }
        }

        {
            LzLogTimeN("", "NEW")
            for( int test=0 ; test<TESTS ; test++ )
            {
                // Compute intertia tensor for tetrahedron
                const double a = lRnd_A[test];
                const Vector3D lCenter( a*0.25, a*0.25, a*0.25 );
                const vector<Point3D> lVers =
                {
                    lRnd_T[test]*(Point3D(0, 0, 0) - lCenter),
                    lRnd_T[test]*(Point3D(a, 0, 0) - lCenter),
                    lRnd_T[test]*(Point3D(0, a, 0) - lCenter),
                    lRnd_T[test]*(Point3D(0, 0, a) - lCenter)
                };

                Eigen::Matrix3d TetI;
                TetraTensor_V2( TetI, lVers, 1.0 );

                if( test+1 == TESTS )
                    LzLogM("", "TetI=\n"<<TetI)
            }
        }
    }

exit(0);


#if 0
        std::uniform_real_distribution<double> lDist(-100.0, +100.0);

//#define DIAG

        double lMaxErr = 0.0;
        for( int t=0 ; t<1000 ; t++ )
        {
#ifndef DIAG
            const double A = lDist(LzServices::RandomEngine());
            const double B = lDist(LzServices::RandomEngine());
            const double C = lDist(LzServices::RandomEngine());
            const double D = lDist(LzServices::RandomEngine());
            //
            const double E = A; //lDist(LzServices::RandomEngine());
            const double F = B; //lDist(LzServices::RandomEngine());
            const double G = C; //lDist(LzServices::RandomEngine());
            const double H = D; //lDist(LzServices::RandomEngine());
#else
            const double A = lDist(LzServices::RandomEngine());
            const double B = lDist(LzServices::RandomEngine());
            const double C = lDist(LzServices::RandomEngine());
            const double D = lDist(LzServices::RandomEngine());
            //
            const double E = A;
            const double F = B;
            const double G = C;
            const double H = D;
#endif

            const LzMath::Polynom P1 = A + _X_*B + _Y_*C + _Z_*D;
            const LzMath::Polynom P2 = E + _X_*F + _Y_*G + _Z_*H;
            const LzMath::Polynom P = P1*P2;

            // Symbolic integration
            LzMath::Polynom P_dZ = P.IntZ();
            P_dZ = P_dZ.Compose( _X_, _Y_, 1.0-_X_-_Y_ ); // No need to subtract Compose( _X_, _Y_, 0.0 ) since it is 0

            LzMath::Polynom P_dZ_dY = P_dZ.IntY();
            P_dZ_dY = P_dZ_dY.Compose( _X_, 1.0-_X_, 0.0 ); // No need to subtract Compose( _X_, 0.0, 0.0 ) since it is 0

            LzMath::Polynom P_dZ_dY_dX = P_dZ_dY.IntX();
            const double lSymVal = P_dZ_dY_dX.EvalAt( 1.0, 0.0, 0.0 ); // No need to subtract EvalAt( 0.0, 0.0, 0.0 ) since it is 0

            // Explicit integration
#ifndef DIAG
            const double lExpVal = (1.0/120.0)*(
                                        20*A*E + 5*A*F + 5*A*G + 5*A*H +
                                        5*B*E + 2*B*F + B*G + B*H +
                                        5*C*E + C*F + 2*C*G + C*H +
                                        5*D*E + D*F + D*G + 2*D*H );
#else
            // Diagonal terms: A = E, B = F, C = G, D = H
            const double lExpVal = (1.0/60.0)*(
                                        10*A*A +
                                        5*A*B + 5*A*C + 5*A*D +
                                        B*B + B*C + B*D + C*C + C*D + D*D );
#endif

//            LzLogM("", "Symbolic val= "<<lSymVal)
//            LzLogM("", "Explicit val= "<<lExpVal)
//            LzLogM("", "")

            double lErr = std::abs( lSymVal - lExpVal );
            if( lMaxErr < lErr )
                lMaxErr = lErr;
        }

        LzLogM("", "*** MAX ERROR= "<<lMaxErr)

    exit(0);
#endif




/// TESTS UNITAIRES : CUBE et CYLINDRE !!
/// TESTS UNITAIRES : CUBE et CYLINDRE !!
/// TESTS UNITAIRES : CUBE et CYLINDRE !!
//    https://hepweb.ucsd.edu/ph110b/110b_notes/node26.html
//    https://scienceworld.wolfram.com/physics/MomentofInertiaCylinder.html
/// TESTS UNITAIRES : CUBE et CYLINDRE !!
/// TESTS UNITAIRES : CUBE et CYLINDRE !!
/// TESTS UNITAIRES : CUBE et CYLINDRE !!


    // Compute inertial tensor of a cube
    Mesh lMesh;
    lMesh.Load(R""(D:\LinkZ\Data\TwInsight\TwinFit\TEST Inertia tensor\CENTERED_UNIT_CUBE.obj)"");
//    lMesh.Load(R""(D:\LinkZ\Data\TwInsight\TwinFit\TEST Inertia tensor\CENTERED_UNIT_CUBE__4.obj)"");
//    lMesh.Load(R""(D:\LinkZ\Data\TwInsight\TwinFit\TEST Inertia tensor\Calibrate_Cyl_HIGH_RES.obj)"");
//    lMesh.Load(R""(D:\LinkZ\Data\TwInsight\TwinFit\TEST Inertia tensor\Calibrate_Cyl_HIGH_RES_ROT_Z.obj)"");
//
//_tmp_DST_T2_in_Anat.obj
//_tmp_SRC_T4_in_Anat.obj
//    lMesh.Load(R""(D:\LinkZ\Data\TwInsight\TwinFit\TEST Inertia tensor\_tmp_SRC_T4_in_Anat.obj)"");
//    lMesh.Load(R""(D:\LinkZ\Data\TwInsight\TwinFit\TEST Inertia tensor\_tmp_DST_T2_in_Anat.obj)"");







    const vector<Point3D> & lVers = lMesh.mVertices;


//    const Point3D lCenter(0, 0, 0);
    const Point3D lCenter(0.5, 0.5, 0.5);
    const Vector3D lCenterVec = lCenter - Point3D(0, 0, 0);


    Eigen::Matrix3d MeshI;
    MeshI.setZero();

    for( size_t t=0 ; t<lMesh.mTriangles.size() ; t++ )
    {
//        // Log
//        LzLogN("", "Tri "<<t)

        const LzTriModel::Triangle & lTri = lMesh.mTriangles[t];

        const vector<Point3D> lTet =
        {
            lCenter,
            lVers[ lTri.mIdxV[0] ] - lCenterVec,
            lVers[ lTri.mIdxV[1] ] - lCenterVec,
            lVers[ lTri.mIdxV[2] ] - lCenterVec,
        };

        Eigen::Matrix3d TetI;
//        TetraTensor( TetI, lTet, 1.0 );
        TetraTensor_V2( TetI, lTet, 1.0 );

//        // Log
//        LzLogM("", "TetI=\n"<<TetI)

        MeshI += TetI;

//        // Log
//        LzLogM("", "MeshI=\n"<<MeshI)
    }

    // Log
    LzLogM("", "MeshI=\n"<<MeshI)

    // Find eigenvectors of inertia tensor
#if 1
// ????? LzMath::ToolBox::EigenJacobi( );
// ????? LzMath::ToolBox::EigenJacobi( );
// ????? LzMath::ToolBox::EigenJacobi( );


    // SVD: H = USVt
    //
    // https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
    // --> Singular values are always sorted in decreasing order
    //
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(MeshI, Eigen::ComputeFullU /*| Eigen::ComputeFullV*/);

    const Eigen::Matrix3cd & U = svd.matrixU();

    // Check that EVs are in decreasing order
    //
    // For discussion regarding sing values vs eigenvalues
    // https://rampure.org/resources/data100/notes/eigen-singular.html
    //
    LzLogM("", "Sing values=\n"<<svd.singularValues())
#else

//dans le cas du cube I calcule au coin, les vecteurs ne sont pas orthogonaux !

    Eigen::EigenSolver<Eigen::Matrix3d> eig( MeshI );
    const Eigen::Matrix3cd U = eig.eigenvectors();

    LzLogM("", "Values=\n"<<eig.eigenvalues())

#endif
    const Eigen::Matrix3cd Ut = U.transpose();
    LzLogM("", "UtU=\n"<<Ut*U)



    LzLogM("", "Eig=\n"<<U)

    Vector3D lEigV[3];
    {
        for( int j=0 ; j<3 ; j++ )
        {
            lEigV[j].X() = U(0,j).real();
            lEigV[j].Y() = U(1,j).real();
            lEigV[j].Z() = U(2,j).real();

            LzLogM("", "Eig["<<j<<"]= "<<lEigV[j].ToString())
            LzLogM("", "Eig["<<j<<"].Norm= "<<lEigV[j].Norm())
        }
    }

    LzLogM("", "Scal= "<<lEigV[0]*lEigV[1])
    LzLogM("", "Scal= "<<lEigV[0]*lEigV[2])
    LzLogM("", "Scal= "<<lEigV[1]*lEigV[2])


    // Convert
    // Use weakest eigen value --> X
    // Use 2nd weakest ev --> y
    // Strongest ev --> z
    RigidTr3D lAlias_I_to_W( lEigV[2], lEigV[1], Point3D(0,0,0) );
    lMesh.RigidTransform( lAlias_I_to_W.Inverse() );
    lMesh.Save(LzServices::StartUpPath()+"/_tmp_NEW.obj");



    {
        LzLogN("", "Checking eigen vector")

//        const Eigen::Vector3d V( 1, 0, 0 );
        const Eigen::Vector3d V( 1, 1, 0 );
        const Eigen::Vector3d W = MeshI * V;

        LzLogM("", "V=\n"<<V)
        LzLogM("", "W=\n"<<W)
    }
    LzLogException("", "*** TESTING")
}
#endif

}

//================================================================================
void InertialPCA( const Mesh & pObject, Point3D & pMean, Vector3D pV[3] )
{
//******************************************************* unit test on unit tet
#if 0
    // Compute intertia tensor for tetrahedron
    const double a = -1.0;
    const Vector3D lCenter( a*0.25, a*0.25, a*0.25 );
    const vector<Point3D> lVers =
    {
        Vector3D(0, 0, 0) - lCenter,
        Vector3D(a, 0, 0) - lCenter,
        Vector3D(0, a, 0) - lCenter,
        Vector3D(0, 0, a) - lCenter,
    };

    Eigen::Matrix3d TetI;
TetraTensor( TetI, lVers, 1.0 );
TetraTensor_V2( TetI, lVers, 1.0 );

    // Log
    LzLogM("", "TetI=\n"<<TetI)

LzLogException("", "*** TESTING")
#endif
//******************************************************* unit test on unit tet


    //---------------------
    // Mass center
    //---------------------

    {
        // Reset mean and total volume
        pMean.Reset();
        double lTotVol = 0;

        // Pick a reasonable center for tetrahedrons' 4th vertex
        const Point3D lD = LzMath::ToolBox::Centroid( pObject.mVertices );

        // Compute volumetric mass center
        for( size_t t=0 ; t<pObject.mTriangles.size() ; t++ )
        {
            // Triangle and vertices
            const Triangle & lTri = pObject.mTriangles[t];
            //
            const Point3D & lA = pObject.mVertices[ lTri.mIdxV[0] ];
            const Point3D & lB = pObject.mVertices[ lTri.mIdxV[1] ];
            const Point3D & lC = pObject.mVertices[ lTri.mIdxV[2] ];

            // Build plane and comp signed tet volume
            const Vector3D lTriNormal = (lC - lA)^(lB - lA); // AC x AB
            const Vector3D lAD = lD - lA;

            // Store signed volume
            const double lSgnVol = (lAD * lTriNormal) / 6.0;

            // Add
            lTotVol += lSgnVol;

            // Centroid
            pMean.X() += lSgnVol * ( lA.X() + lB.X() + lC.X() + lD.X() ) / 4.0;
            pMean.Y() += lSgnVol * ( lA.Y() + lB.Y() + lC.Y() + lD.Y() ) / 4.0;
            pMean.Z() += lSgnVol * ( lA.Z() + lB.Z() + lC.Z() + lD.Z() ) / 4.0;
        }

        // Check
        if( LzMath::ToolBox::IsZero(lTotVol) )
            LzLogException("", "Cannot compute intertia tensor! Found a null volume= "<<lTotVol)

        // Normalize
        pMean /= lTotVol;

//// Log
//LzLogM("", "Mass center= "<<pMean.ToString())
//LzLogM("", "Tot vol= "<<lTotVol)
    }


    //---------------------
    // Inertia tensor
    //---------------------

    // Set tensor
    Eigen::Matrix3d I;
    I.setZero();

    // Scale coordinates to mm
    const double SCALE = pObject.IsMeters() ? 1000.0 : 1.0 ;

    // Mean offset vector
    const Vector3D lMean = pMean - Point3D(0,0,0);

    // Accumulate inertia tensors of all tetrahedrons
    for( size_t t=0 ; t<pObject.mTriangles.size() ; t++ )
    {
        // Tetrahedron vertices corrected for mass center
        const Triangle & lTri = pObject.mTriangles[t];
        //
        const vector<Point3D> lVers
        {
            Point3D(0, 0, 0), // pMean - pMean
            (pObject.mVertices[ lTri.mIdxV[0] ] - lMean)*SCALE,
            (pObject.mVertices[ lTri.mIdxV[1] ] - lMean)*SCALE,
            (pObject.mVertices[ lTri.mIdxV[2] ] - lMean)*SCALE,
        };

        // Compute intertia tensor for tetrahedron
        Eigen::Matrix3d TetI;
        TetraTensor_V2( TetI, lVers, 0.001/* A density of 0.001 will produce TetI matrices with reasonable values */ );

//// Log
//LzLogM("", "TetI=\n"<<TetI)

        // Update global tensor
        I += TetI;
    }

//// Log
//LzLogM("", "Total inertia tensor=\n"<<I)


#if 0
    NON !! les vecteurs ne sont pas orthogonaux

    // Find eigenvectors of inertia tensor
    Eigen::EigenSolver<Eigen::Matrix3d> eig( I );
LzLogM("", "Eig=\n"<<eig.eigenvectors())

    const Eigen::Matrix3cd & U = eig.eigenvectors();

    pV[0] = Vector3D( U( 0, 0 ).real(), U( 1, 0 ).real(), U( 2, 0 ).real() );
    pV[1] = Vector3D( U( 0, 1 ).real(), U( 1, 1 ).real(), U( 2, 1 ).real() );
    pV[2] = Vector3D( U( 0, 2 ).real(), U( 1, 2 ).real(), U( 2, 2 ).real() );

//    pV[0] = Vector3D( U( 0, 0 ).real(), U( 0, 1 ).real(), U( 0, 2 ).real() );
//    pV[1] = Vector3D( U( 1, 0 ).real(), U( 1, 1 ).real(), U( 1, 2 ).real() );
//    pV[2] = Vector3D( U( 2, 0 ).real(), U( 2, 1 ).real(), U( 2, 2 ).real() );

//    Eigen::Matrix3cd D = eig.eigenvalues();
//LzLogM("", "D=\n"<<eig.eigenvalues()) //********************************* eigenvalues are not sorted!!!

//LzLogM("", "First prod=\n"<<I*eig.eigenvectors()(0))
//LzLogM("", "Check prod=\n"<<eig.eigenvalues()(0)*eig.eigenvectors()(0))

//        eig.eigenvalues()(0);
#else

    // SVD: H = USVt
    //
    // https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
    // --> Singular values are always sorted in decreasing order
    //
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(I, Eigen::ComputeFullU /*| Eigen::ComputeFullV*/);

    // Read eigenvectors
    const Eigen::Matrix3d & U = svd.matrixU();
#endif

    // Vectors produced by EigenJacobi are already normalized
    pV[0] = Vector3D( U( 0, 0 ), U( 1, 0 ), U( 2, 0 ) );
    pV[1] = Vector3D( U( 0, 1 ), U( 1, 1 ), U( 2, 1 ) );
    pV[2] = Vector3D( U( 0, 2 ), U( 1, 2 ), U( 2, 2 ) );
}

//================================================================================
void TEST_InertialPCA()
{
    // Using data from the article to check implementation

    vector<Point3D> lTet =
    {
        Point3D(8.33220, -11.86875, 0.93355),
        Point3D(0.75523, 5.00000, 16.37072),
        Point3D(52.61236, 5.00000, -5.38580),
        Point3D(2.00000, 5.00000, 3.00000),
    };

    Point3D lG( 0, 0, 0 );
    for( const Point3D & iP : lTet )
    {
        lG.X() += iP.X();
        lG.Y() += iP.Y();
        lG.Z() += iP.Z();
    }
    lG /= 4;
    LzLogM("", "G= "<<lG.ToString())

//    const Point3D lG( 15.92492, 0.78281, 3.72962 );

    // Offset by G
    for( Point3D & iP : lTet )
    {
        iP.X() -= lG.X();
        iP.Y() -= lG.Y();
        iP.Z() -= lG.Z();
    }

    // Inertia tensor
    Eigen::Matrix3d I;
    TetraTensor( I, lTet, 1.0 );

    // Log
    LzLogM("", "I=\n"<<std::setprecision(12)<<I)
}

}
