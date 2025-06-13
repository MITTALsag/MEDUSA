#include "OutlinerTools_Meshing.h"
#include "Outliner.h"
#include <LzMath/ToolBox.h>
#include <LzGeom/Line3D.h>
#include <Eigen/Dense>
#include <iostream>
#include <csignal>
#include <unistd.h>
#include <map>
#define USE_ASYNC_RECURSION
#ifdef USE_ASYNC_RECURSION
#include <LzAsync/DispatchQ.h>
#endif


namespace LzTriModel
{
namespace OutlinerTools
{

//================================================================================
static double MinSegToSeg( const Point3D & pA, const Point3D & pB,
                           const Point3D & pC, const Point3D & pD,
                           Point3D & pI, Point3D & pJ )
{
//// Log
//LzLogN("", "MinSegToSeg")
//LzLogM("", "A= "<<pA.ToString())
//LzLogM("", "B= "<<pB.ToString())
//LzLogM("", "C= "<<pC.ToString())
//LzLogM("", "D= "<<pD.ToString())

    // Vectors
    const Vector3D AB = pB - pA;
    const Vector3D CD = pD - pC;

    // Try to compute the angle between the vectors
    double lAlpha_Rad;
    const bool lHasAlpha = AB.ShortPositiveAngleTo_NoExc( CD, lAlpha_Rad );

//LzLogM("", "AB= "<<AB.Norm())
//LzLogM("", "CD= "<<CD.Norm())
//LzLogM("", "Angle= "<<LzMath::RAD_2_DEG*AB.ShortPositiveAngleTo(CD))

    // Parallel segments?
    if( !lHasAlpha
     || std::abs(lAlpha_Rad -          0)<1e-6
     || std::abs(lAlpha_Rad - LzMath::PI)<1e-6 )
    {
//LzLogM("", "PARALLEL SEGMENTS")

        // Project C and D onto the AB line
        const Line3D AB_line( pA, pB );
        const Point3D C_p = AB_line.Projection( pC );
        const Point3D D_p = AB_line.Projection( pD );

        // Get length of AB from the matrix M
        const double AB_norm = AB.Norm(); //std::sqrt( M(0, 0) );
        const double CD_norm = CD.Norm(); //std::sqrt( M(1, 1) );

        // Mid-points on line (A, B)
        const Point3D AB_mid = pA.MidPointTo( pB );
        const Point3D CD_p_mid = C_p.MidPointTo( D_p );

        // Test for overlapping of projections
        if( AB_mid.DistanceTo(CD_p_mid) <= 0.5*(AB_norm + CD_norm) )
        {
            //==> Yes, parallel and projections are overlapping

            const Point3D * lpOuters[6][2] =
            {
                { &pA, &C_p }, { &pA, &D_p },
                { &pB, &C_p }, { &pB, &D_p },
                { &pA, &pB }, { &C_p, &D_p },
            };

            // Find max dist between non-overlapping ends
            double lMaxDist = 0.0;
            int lMax_e = 0;
            for( int e=0 ; e<6 ; e++ )
            {
                const double lDist = lpOuters[e][0]->DistanceTo( *lpOuters[e][1] );
                if( lMaxDist < lDist )
                {
                    lMaxDist = lDist;
                    lMax_e = e;
                }
            }

            // Inner ends on line (AB) corresponding to the outer end points
            // We must use projected points here since we will be computing the
            // midpoint between inner points (on (AB))
            const Point3D * lpInners[6][2] =
            {
                { &pB, &D_p }, { &pB, &C_p },
                { &pA, &D_p }, { &pA, &C_p },
                { &C_p, &D_p }, { &pA, &pB },
            };

            // Point I is the midpoint between the inner points on (AB)
            pI = lpInners[lMax_e][0]->MidPointTo( *lpInners[lMax_e][1] );

            // Point J is point I projected back on (CD)
            pJ = Line3D{pC, pD}.Projection( pI );

            // Distance can be 0 if A, B, C, D are aligned, or
            // it will be the distance between the 2 lines (A, B) and (C, D)
            return pI.DistanceTo( pJ );
        }
        else
        {
            //==> No, parallel but not overlapping

            // Possible pairs of segments' ends
            const Point3D * lpEnds[4][2] =
            {
                { &pA, &pC }, { &pA, &pD },
                { &pB, &pC }, { &pB, &pD },
            };

            // Find min dist between non-overlapping ends
            double lMinDist = +std::numeric_limits<double>::max();
            for( int e=0 ; e<4 ; e++ )
            {
                const double lDist = lpEnds[e][0]->DistanceTo( *lpEnds[e][1] );
                if( lMinDist > lDist )
                {
                    lMinDist = lDist;
                    pI = *lpEnds[e][0];
                    pJ = *lpEnds[e][1];
                }
            }

            return lMinDist;
        }
    }
    else
    {
//LzLogM("", "INTERSECTING SEGMENTS")

        // Build matrix
        Eigen::Matrix2d M;
        M(0, 0) = AB.SquaredNorm();
        M(0, 1) = M(1, 0) = -AB*CD;
        M(1, 1) = CD.SquaredNorm();

//// Compute det
//const double lDet = M.determinant();
//LzLogM("", "M=\n"<<M)
//LzLogM("", "Det= "<<lDet)

        // Assemble right hand side
        Eigen::Vector2d Y;
        Y(0) = -( (pA.X() - pC.X())*AB.X() + (pA.Y() - pC.Y())*AB.Y() + (pA.Z() - pC.Z())*AB.Z() );
        Y(1) = +( (pA.X() - pC.X())*CD.X() + (pA.Y() - pC.Y())*CD.Y() + (pA.Z() - pC.Z())*CD.Z() );

        // Solve
        Eigen::Vector2d X = M.inverse() * Y;

//LzLogM("", "Y(0)= "<<Y(0))
//LzLogM("", "Y(1)= "<<Y(1))
////
//LzLogM("", "X(0)= "<<X(0))
//LzLogM("", "X(1)= "<<X(1))
////
//LzLogM("", "M*X=\n"<<M*X)

        // Clamp
        if( X(0) < 0.0 ) X(0) = 0.0;
        if( X(0) > 1.0 ) X(0) = 1.0;
        //
        if( X(1) < 0.0 ) X(1) = 0.0;
        if( X(1) > 1.0 ) X(1) = 1.0;

        // Compute I et J
        pI = pA + X(0)*AB;
        pJ = pC + X(1)*CD;

//LzLogM("", "X(0)= "<<X(0)<<" --> pI= "<<pI.ToString())
//LzLogM("", "X(1)= "<<X(1)<<" --> pJ= "<<pJ.ToString())
//LzLogM("", "Dist= "<<pI.DistanceTo( pJ ))

        return pI.DistanceTo( pJ );
    }
}

//================================================================================
static double MinPntToSeg( const Point3D & pPnt,
                           const Point3D & pA, const Point3D & pB,
                           Point3D & pJ )
{
    // Check
    // ... not null dist

    // Buid line and project
    const Vector3D u_AB = pB - pA;
    //
    const Line3D AB( pA, u_AB );
    pJ = AB.Projection( pPnt );

    // Vector A to J
    const Vector3D u_AJ = pJ - pA;

    // Check position of J w.r.t. AB
    if( u_AB*u_AJ >= 0 )
    {
        //==> J is pointed to by AB, meaning J is in [A, B)

        // Outside of AB?
        if( u_AJ.SquaredNorm() > u_AB.SquaredNorm() )
        {
            // Clamp
            pJ = pB;
        }
        else
        {
            // Leave pJ where it is: inside the segment
        }
    }
    else
    {
        //==> J is NOT in [A, B), meaning J is strictly "below" A

        // Meaning J is outside of AB ==> clamp
        pJ = pA;
    }

    return pPnt.DistanceTo( pJ );
}

//================================================================================
static bool HasSelfIntersections( const List<Point3D> & pOut, double pEpsilon )
{
    // Check closedness
    const bool lIsClosed = Outliner::IsClosed( pOut );

    // Here points count = N >= 2 (open or closed outline)
    const size_t N = pOut.Count();

    // If N == 2 and closed ==> same point: no self intersection
    // If N == 2 and open ==> single segment: no self intersection
    if( pOut.Count() == 2 )
        return false;

    // Here count = N >= 3

    // Check all segments, using Epsilon to define what is Zero and what is not
    for( size_t i=0 ; i<=N-3 ; i++ )
    {
        // Read first segment
        const Point3D & A = pOut[  i  ];
        const Point3D & B = pOut[i + 1];

        for( size_t j=i+1 ; j<=N-2 ; j++ )
        {
            // Read second segment
            const Point3D & C = pOut[  j  ];
            const Point3D & D = pOut[j + 1];

            // Compute min distance seg to seg, and projection points
            Point3D I, J;
            const double lMin_AB_to_CD = MinSegToSeg( A, B, C, D, I, J );

            // Below threshold?
            if( LzMath::ToolBox::IsZero(lMin_AB_to_CD, pEpsilon) )
            {
                // Possible intersection! Check exclusions

                // Exclusion 1: B = C: AB -- CD = A(BC)D
                if( j == i+1 )
                {
                    // Check position of I and J: must be in C (within pEpsilon)
                    if( LzMath::ToolBox::IsZero( C.DistanceTo(I), pEpsilon )
                     && LzMath::ToolBox::IsZero( C.DistanceTo(J), pEpsilon ) )
                    {
                        // NO collision
                        continue;
                    }
                }
                else
                // Exclusion 2: CLOSED outline and A = D: CD -- AB = C(DA)B
                if( lIsClosed && i==0 && j==N-2 )
                {
                    // Check position of I and J: must be in A (within pEpsilon)
                    if( LzMath::ToolBox::IsZero( A.DistanceTo(I), pEpsilon )
                     && LzMath::ToolBox::IsZero( A.DistanceTo(J), pEpsilon ) )
                    {
                        // NO collision
                        continue;
                    }
                }

                // Intersection confirmed
//LzLogM("", "HasSelfIntersections (closed= "<<std::boolalpha<<lIsClosed<<"): YES")
                return true;
            }
        }
    }

    // No self intersections
//LzLogM("", "HasSelfIntersections (closed= "<<std::boolalpha<<lIsClosed<<"): NO")
    return false;
}

//================================================================================
#ifdef USE_ASYNC_RECURSION
static std::shared_ptr<LzAsync::DispatchQ> spDQ;
static std::mutex sMeshTriMex;
#endif
//
//#define ADD_MIDPOINT_ON_SPLIT

//================================================================================
template<class T> class FlipFlopper
{
public:
    FlipFlopper( T pMin, T pMax )
    {
        mMin = pMin;
        mMax = pMax;
        mI0 = (pMin + pMax) / 2; // Integer division
        mDelta = 0;
        mSgn = Sgn::Neg;
    }

    T Init() const { return mI0; }

    bool Next( T & pNext )
    {
        switch( mSgn )
        {
        case Sgn::Pos:
            pNext += mDelta++;
            mSgn = Sgn::Neg;
            return pNext <= mMax;

        case Sgn::Neg:
            pNext -= mDelta++;
            mSgn = Sgn::Pos;
            return pNext >= mMin;
        }
    }

protected:
    enum class Sgn { Pos, Neg };
    Sgn mSgn;

    T mMin;
    T mMax;

    T mI0;
    T mDelta;
};


//================================================================================
static bool IsInner( const Point3D & pA, const Point3D & pB,
                     const Point3D & pA_prev, const Point3D & pA_next,
                     const Vector3D & pNor )
{
    // Fast computation of INNER test -- when there are NO self-intersections

    // A_prev --> A
    const Vector3D A_prev_A = pA - pA_prev;
    try
    {
        // Compute angles
        const double alpha = A_prev_A.SignedAngleTo((pA_next - pA), pNor);
        const double beta  = A_prev_A.SignedAngleTo((pB - pA), pNor);

        // Test
        return beta>alpha && beta<LzMath::PI;
    }
    catch(...)
    {
        // Could not compute angles: NOT INNER
        return false;
    }
}

//================================================================================
static bool IsInner( const Point3D & pA, const Point3D & pB,
                     const vector<size_t> & pIdx,
                     const Mesh & pTo )
{
    // Slow but safer computation of INNER test -- when there ARE self-intersections

    //
    // Here: pIdx represents a closed outline but the first point IS NOT REPEATED at the end.
    // =====
    //

    // Get the normal vector
    const Vector3D & lNor = pTo.mNormals[0];

    // Get the number of indices
    const size_t N = pIdx.size();

    // Sample 3 points along the segment [A, B]
    // All 3 points must be IN
    const Vector3D AB = pB - pA;
    const Point3D lSamplePts[3] =
    {
        pA + 0.25*AB,  // 1 / 4
        pA + 0.50*AB,  // 2 / 4
        pA + 0.75*AB,  // 3 / 4
    };

    // All sample points must be Inner
    for( int s=0 ; s<3 ; s++ )
    {
        // Read sample point
        const Point3D & lSampPt = lSamplePts[s];

        // Compute sum of angles to determine in/out
        double lAngleSum = 0;

        // Check all segments in closed outline: NEED to repeat last (w to go around the list
        for( size_t i=0 ; i<N ; i++ )
        {
            // Two points
            const Point3D & lV = pTo.mVertices[ pIdx[i] ];
            const Point3D & lW = pTo.mVertices[ pIdx[(i + 1) % N] ];

            // Two vectors
            const Vector3D lPtV = lV - lSampPt;
            const Vector3D lPtW = lW - lSampPt;

            // Check if point lies near one of the two vertices v or w
            if( LzMath::ToolBox::IsZero(lPtV.Norm() * lPtW.Norm()) )
            {
                // Point on boundary: NOT inner
                return false;
            }

            // Compute angle: it is safe since the product of both norms above is not null (within the SAME default tolerance)
            const double lAngle = lPtV.SignedAngleTo( lPtW, lNor );

            // Check if point lies on edge [v,w] iff. angle close to +Pi or -Pi
            if( LzMath::ToolBox::IsZero( lAngle - LzMath::PI, 1e-6 )
             || LzMath::ToolBox::IsZero( lAngle + LzMath::PI, 1e-6 ) )
            {
                // Point on boundary: NOT inner
                return false;
            }

            // Point inside or outside
            lAngleSum += lAngle;
        }

        // Absolute value
        double lAbsAngleSum = std::abs( lAngleSum );

        // For the sample point to be "inside", the sum of all angles must be close to +/-2 Pi

        // Reject point if it is Outer i.e., does not comply with this constraint
        if( !LzMath::ToolBox::IsZero( lAbsAngleSum - LzMath::DPI, 1e-6 ) )
        {
            // Outer
            return false;
        }
    }

    // All tests succeeded, split [A, B] is Inner
    return true;
}

//================================================================================
static void MeshOutline_Recursive( map<std::thread::id, Mesh> & pTo,
                                   const std::shared_ptr<const vector<size_t>> ppIdx,
                                   MeshOutlineOptions & pOptions )
{
    //--------------------------------------
    // Reached a leaf?
    //--------------------------------------

    if( ppIdx->size() == 3 )
    {
#ifdef USE_ASYNC_RECURSION
        // Acquire mutex
        // LzAsync::MutexNeed(*spDQ, "sMeshTriMex");
        // std::unique_lock<std::mutex> lLock( sMeshTriMex );
        // LzAsync::MutexHave(*spDQ, "sMeshTriMex");


#endif
        // Create final triangle and leave (no pun intended)
        pTo.at(std::this_thread::get_id()).mTriangles.emplace_back((*ppIdx)[0], (*ppIdx)[1], (*ppIdx)[2], 0, 0, 0);
        // LzAsync::MutexFree(*spDQ, "sMeshTriMex");
        return;
    }
    // std::this_thread::sleep_for( std::chrono::milliseconds(2) );

    //
    // Nope... need to subdivide
    //

    //--------------------------------------
    // Read some useful data
    //--------------------------------------
    TagStart(*spDQ, "Balise Test FlipFlopper");

    // Number of (unique) indices (at this stage the first point is not repeated anymore at the end)
    const size_t N = ppIdx->size();

    // Read normal (used in internal check)
    const Vector3D & lNor = pTo.at(std::this_thread::get_id()).mNormals[0];

    // Log

//    // Optional Node: simple
//    // Optional Node: simple + timer
//     LzOptLogTimeN(ENABLED,THREAD,MSG)\
//        if( ENABLED ) LINKZ_LOG_TXT("/-- "<<THREAD_INFO(THREAD)<<" "<<MSG)\
//        LzServices::LzLogNodeObject lLzLogNodObj(ENABLED, true);


//// Log
//LzLogN("", "MeshOutline_Recursive. N= "<<N<<" --- CallId= "<<pOptions.mCallId++)

    //--------------------------------------
    // Future best split: to be found!
    //--------------------------------------

    std::pair<size_t, size_t> lBestSplit;

    //--------------------------------------
    // Evaluate all possible splits
    //--------------------------------------
    // Score for a split = MINIMAL distance between split and all outline edges
    using Split = LzServices::Sortable<std::pair<size_t, size_t>, double>;
    List<Split> lSplits;
    {
//// Log
//LzLogN("", "Evaluating all possible splits")

        // Unused projection points
        Point3D I, J;
        for( size_t i=0 ; i<=N-3 ; i++ )
        {
            // Get point
            const Point3D & A = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[i] ];
            const size_t last_j = i==0 ? N-2 : N-1 ;
#if 1
            FlipFlopper<size_t> lFF( i+2, last_j );
            for( size_t j=lFF.Init() ; lFF.Next(j) ;)
#else
            for( size_t j=i+2 ; j<=last_j ; j++ )
#endif
            {
                // Get point
                const Point3D & B = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[j] ];

                // Compute max of min dists
                double lMinDist = +std::numeric_limits<double>::max();
                for( size_t k0=0 ; k0<N ; k0++ )
                {
                    // k+1 modulo N
                    const size_t k1 = (k0 + 1) % N;

                    // Compute dist between edge [k0, k1] and split [i, j]
                    double lDist;

//*** NB: some of these computations are done multiple times! Reduce redundancy
//*** NB: some of these computations are done multiple times! Reduce redundancy
//*** NB: some of these computations are done const Vector3D & lNor = pTo[std::this_thread::get_id()].mNormals[0];multiple times! Reduce redundancy
//*** NB: some of these computations are done multiple times! Reduce redundancy

                    // Split-edge sharing vertex?
                    if( k0==i || k0==j )
                    {
                        // Get point, compute dist
                        const Point3D & C = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[k1] ];
                        lDist = MinPntToSeg( C, A, B, I );
                    }
                    else
                    if( k1==i || k1==j )
                    {
                        // Get point, compute dist
                        const Point3D & C = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[k0] ];
                        lDist = MinPntToSeg( C, A, B, I );
                    }
                    else
                    {
                        // Get points, compute dist
                        const Point3D & C = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[k0] ];
                        const Point3D & D = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[k1] ];
                        lDist = MinSegToSeg(A, B, C, D, I, J);
                    }

                    // Update min dist
                    if( lMinDist > lDist )
                        lMinDist = lDist;
                }

                // Create new split
                const Split lNewSplit{{i, j}, lMinDist};

                //....................
                // Early stop?
                //....................

                if( pOptions.mEarlyStopEvaluation
                 && lMinDist >= pOptions.mEarlyStopMinDist )
                {
//                    // Log
//                    LzLogN("", "Using early stop. Min dist= "<<lMinDist)

                    // Fast or Slow inner test?
                    bool lIsInner;
                    if( pOptions.mHasSelfInters )
                    {
                        //===> Outline has self-intersections: SLOW MODE

                        // Test
                        lIsInner = IsInner( A, B, *ppIdx, pTo.at(std::this_thread::get_id()) );
                    }
                    else
                    {
                        //===> Outline has NO self-intersections: FAST MODE

                        // Read prev and next points
                        const Point3D & A_prev = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[(N + i - 1) % N] ];
                        const Point3D & A_next = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[(    i + 1) % N] ];

                        // Test
                        lIsInner = IsInner( A, B, A_prev, A_next, lNor );
                    }

                    // Check
                    if( lIsInner )
                    {
                        // Set split as best split and skip:
                        // - the rest of the exhaustive evaluation
                        // - the search for best internal split below
                        lBestSplit = lNewSplit.mT;
                        goto FoundBestInternalSplit;
                    }
                }
                //.................... Early stop

                // Stash split: BIG min dists FIRST, SMALL min dists LAST
                lSplits.AddIntoDecrList( lNewSplit, /*pUnique= whatever... who cares... entries will always be unique*/true );
            }
        }
    }


    //--------------------------------------
    // Find the best INTERNAL split
    //--------------------------------------

    {
//        // Log
//        LzLogN("", "Looking for the best possible internal split")

        // Start from the top i.e., the ones with the GREATEST (safest) distance to split
        BrowseList( iS, lSplits )
        {
            // Pick split
            lBestSplit = lSplits.GetAt(iS).mT;

            // Test if split is an INNER split

            // Read indices
            const size_t i = lBestSplit.first;
            const size_t j = lBestSplit.second;

            // Read split points
            const Point3D & A = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[i] ];
            const Point3D & B = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[j] ];

            // Fast or Slow inner test?
            bool lIsInner;
            if( pOptions.mHasSelfInters )
            {
                //===> Outline has self-intersections: SLOW MODE

                // Test
                lIsInner = IsInner( A, B, *ppIdx, pTo.at(std::this_thread::get_id()) );
            }
            else
            {
                //===> Outline has NO self-intersections: FAST MODE

                // Read prev and next points
                const Point3D & A_prev = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[(N + i - 1) % N] ];
                const Point3D & A_next = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[(    i + 1) % N] ];

                // Test
                lIsInner = IsInner( A, B, A_prev, A_next, lNor );
            }

            // Check
            if( lIsInner )
                goto FoundBestInternalSplit;
        }
    }

    //**********************************
    //***                            ***
    //*** Could not find inner split ***
    //***                            ***
    //**********************************

    // Possible self intersections?
    if( pOptions.mHasSelfInters )
    {
        // Log
        LzLogM("", "Could not find an internal split. Self-intersections are possible => picking the best non-internal split.")

        // Pick the first one (since they are all external anyway)
        lBestSplit = lSplits.GetHead().mT;
    }
    else
    {
        // This should not happen
        LzLogErr("", "Could not find an internal split in a mesh without self-intersections!")

        // Convert local indices to outline in order to compute orientation
        List<Point3D> lLocalOut;
        for( size_t v=0 ; v<N ; v++ )
            lLocalOut.AddTail( pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[v] ] );
        lLocalOut.AddTail( pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[0] ] );

        // Check orientation CCW/CW/?
        const LzTriModel::Outliner::Orientation lOr = LzTriModel::Outliner::GetOrientation( lLocalOut, lNor );
        if( lOr != LzTriModel::Outliner::Orientation::CCW )
            LzLogErr("", "Unexpected orientation: "<<LzTriModel::Outliner::OrientationToString(lOr)<<" in local outline.")

        // Log vertex coordinates for debug
        LzLogN("", "*** Faulty outline: N= "<<N<<" unique vertices.")
        for( size_t v=0 ; v<N ; v++ )
            LzLogM("", pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[v] ].ToString())

#if 1
        // FAIL noisily
        LzLogE("", "NOISY fail! Fatal error...")

        // Fatal error
        LzLogException("", "Could not find an internal split in a mesh without self-intersections!")
#else
        // FAIL silently
        LzLogE("", "SILENT fail! Picking the first split in the list.")

        // Pick the first one (since they are all external anyway)
        lBestSplit = lSplits.GetHead().mT;
#endif
    }

FoundBestInternalSplit:

    //--------------------------------------
    // Extract sub-indices and recurse
    //--------------------------------------

    // Direct indices
    std::shared_ptr<vector<size_t>> lpSubIdx_0 = std::make_shared<vector<size_t>>();
    for( size_t k=lBestSplit.first ; k<=lBestSplit.second ; k++ )
        lpSubIdx_0->push_back( (*ppIdx)[k] );
    //
    // Job name
    const string lJobName_0 = "Mesh ["+std::to_string(ppIdx->back())+", "+std::to_string(ppIdx->front())+"] ==> "
                                   "["+std::to_string(lpSubIdx_0->back())+", "+std::to_string(lpSubIdx_0->front())+"]";

    // Circular indices
    std::shared_ptr<vector<size_t>> lpSubIdx_1 = std::make_shared<vector<size_t>>();
    for( size_t k=lBestSplit.second ; k<=(N + lBestSplit.first) ; k++ )
        lpSubIdx_1->push_back( (*ppIdx)[k % N] );
    //
    // Job name
    const string lJobName_1 = "Mesh ["+std::to_string(ppIdx->back())+", "+std::to_string(ppIdx->front())+"] ==> "
                                   "["+std::to_string(lpSubIdx_1->back())+", "+std::to_string(lpSubIdx_1->front())+"]";


#ifdef ADD_MIDPOINT_ON_SPLIT // Add split mid point
    {
        static int sIdx = 0;

        // Read indices
        const size_t i = lBestSplit.first;
        const size_t j = lBestSplit.second;

        // Read split points
        const Point3D & A = pTo.mVertices[ (*ppIdx)[i] ];
        const Point3D & B = pTo.mVertices[ (*ppIdx)[j] ];

        // Criterion = exclusion radius
        const Point3D & C = A.MidPointTo( B );

        bool lAllVersFarEnough = true;
        for( size_t v=0 ; v<ppIdx->size() ; v++ )
        {
            const Point3D & lVer = pTo.at(std::this_thread::get_id()).mVertices[ (*ppIdx)[v] ];

            if( C.DistanceTo( lVer ) < 2.0 )
            {
                lAllVersFarEnough = false;
                break;
            }
        }

        if( lAllVersFarEnough )
        {
//            LzLogN("", "Best split= "<<i<<" -- "<<j<<".")
//            LzLogM("", "Dist A["<<(*ppIdx)[i]<<"] B["<<(*ppIdx)[j]<<"]= "<<A.DistanceTo(B))

#ifdef USE_ASYNC_RECURSION
            // Acquire mutex
            std::unique_lock<std::mutex> lLock( sMeshTriMex );
#endif
            const size_t lMidPtIdx = pTo.mVertices.size();
            pTo.at(std::this_thread::get_id()).mVertices.push_back( C );

            lpSubIdx_0->push_back( lMidPtIdx );
            lpSubIdx_1->push_back( lMidPtIdx );

//pTo.Save(LzServices::StartUpPath()+"/__DEBUG__pTo__"+std::to_string(sIdx)+".obj");
//{
//    LzLogN("", "lpSubIdx_0 "<<sIdx)
//    for( size_t i : *lpSubIdx_0 )
//        LzLogM("", i)
//}
//{
//    LzLogN("", "lpSubIdx_1 "<<sIdx)
//    for( size_t i : *lpSubIdx_1 )
//        LzLogM("", i)
//}

sIdx++;

        }
    }
#endif



//    LzLogM("", "Splitting "<<N<<" into "<<lpSubIdx_0->size()<<" + "<<lpSubIdx_1->size()<<".")




    // Recursion
#ifdef USE_ASYNC_RECURSION
    TagStop(*spDQ, "Balise Test FlipFlopper");
    spDQ->dispatch( { {lJobName_0, [&pTo, lpSubIdx_0, &pOptions]{ MeshOutline_Recursive( pTo, lpSubIdx_0, pOptions ); }},
                      {lJobName_1, [&pTo, lpSubIdx_1, &pOptions]{ MeshOutline_Recursive( pTo, lpSubIdx_1, pOptions ); }} } );

//list of jobs ... ou std::initializer list ???
//voir si on peut tout passer par params (et non options) et lister les indices plutot que copier des tableaux


//    spDQ->dispatch( [&pTo, lpSubIdx_1, &pOptions]{ MeshOutline_Recursive( pTo, lpSubIdx_1, pOptions ); } );

#else
    {
        MeshOutlineOptions lOpt_0 = pOptions;
        lOpt_0.mJobName = lJobName_0;
        MeshOutline_Recursive( pTo, lpSubIdx_0, lOpt_0 );
    }
    {
        MeshOutlineOptions lOpt_1 = pOptions;
        lOpt_1.mJobName = lJobName_1;
        MeshOutline_Recursive( pTo, lpSubIdx_1, lOpt_1 );
    }

//chrono de chaque tache en non //


//faire meme test sur un contour tout bete (cercle)

#endif
}

//================================================================================
void MeshOutline( Mesh & pTo,
                  const List<Point3D> & pOut,
                  const Vector3D & pNormalCCW,
                  const MeshOutlineOptions & pOptions/*=MeshOutlineOptions()*/ )
{

    // Check
    if( !LzTriModel::Outliner::IsClosed(pOut) )
        LzLogException("", "Cannot mesh an open outline!")

    // Check
    if( pOut.Count() < 3 )
        LzLogException("", "Cannot mesh an outline with less than 3 (unique) vertices!")

// find usage ccw cw et undefined
// find usage ccw cw et undefined
// find usage ccw cw et undefined

    // Restrict to actual N UNIQUE points (ignore the last repeated point)
    const size_t N = pOut.Count() - 1;

    // Log
    LzLogTimeN("", "Meshing closed outline containing "<<N<<" unique vertices.")



#ifdef USE_ASYNC_RECURSION
    //--------------------------------------
    // Set-up async mode
    //--------------------------------------

    // Create DQ
    spDQ = std::make_shared<LzAsync::DispatchQ>( "MeshOutline -- DQ", 8 );
//    spDQ = std::make_shared<LzAsync::DispatchQ>( "MeshOutline -- DQ", 16 );

    //--------------------------------------
    // Initialize map of meshes
    //--------------------------------------
    map<std::thread::id, Mesh> lTo;
    vector<std::thread::id> lThreadIds = spDQ->GetThreadIds();

    for ( const std::thread::id & iTh : lThreadIds )
    {   
        // LzLogM("test", "id :" << iTh)
        // Clean previous, if any
        lTo[iTh].Free();
        //
        // Vertices
        lTo[iTh].mVertices.resize( N );
        size_t v = 0;
        BrowseList( iP, pOut )
        {
            lTo[iTh].mVertices[v++] = pOut.GetAt(iP);
            if( v == N )
                break;
        }
        //
        lTo[iTh].mNormals = { pNormalCCW };
        //
        // lTolTo[iTh].mTriangles: left EMPTY
    }
#endif


    //--------------------------------------
    // Set-up options
    //--------------------------------------

    // Copy user-provided options, add private options
    MeshOutlineOptions lOptions = pOptions;
    lOptions.mCallId = 0;
    lOptions.mHasSelfInters = HasSelfIntersections( pOut, lOptions.mSelfInterEpsilon );

    // Log
    LzLogM("", "Found self intersections= "<<std::boolalpha<<lOptions.mHasSelfInters)

    //--------------------------------------
    // Create initial index table
    //--------------------------------------

    // Compute orientation
    const LzTriModel::Outliner::Orientation lOrient = LzTriModel::Outliner::GetOrientation( pOut, pNormalCCW );

    // Log
    LzLogM("", "Found orientation= "<<LzTriModel::Outliner::OrientationToString(lOrient))

    // Init to 'identity' mapping
    std::shared_ptr<vector<size_t>> lpIdx = std::make_shared<vector<size_t>>( N );
    if( lOrient == LzTriModel::Outliner::Orientation::CCW )
    {
        for( size_t v=0 ; v<N ; v++ )
            (*lpIdx)[v] = v;
    }
    else // CW or Undefined (the latter will be patched after the meshing)
    {
        for( size_t v=0 ; v<N ; v++ )
            (*lpIdx)[v] = N-1-v;

        // Check that an Undefined orientation can only be obtained if there are self-intersections
        if (lOrient == LzTriModel::Outliner::Orientation::Undefined 
            && !lOptions.mHasSelfInters)
        {

            // This should not happen
            LzLogErr("", "Found an Undefined orientation but NO self intersections!")

            // Convert local indices to outline in order to compute orientation
            List<Point3D> lLocalOut;
            for( size_t v=0 ; v<N ; v++ )
                lLocalOut.AddTail( pTo.mVertices[ (*lpIdx)[v] ] );
            lLocalOut.AddTail( pTo.mVertices[ (*lpIdx)[0] ] );

            // Check orientation CCW/CW/?
            const LzTriModel::Outliner::Orientation lOr = LzTriModel::Outliner::GetOrientation( lLocalOut, pNormalCCW );
            if( lOr != LzTriModel::Outliner::Orientation::CCW )
                LzLogErr("", "Unexpected orientation: "<<LzTriModel::Outliner::OrientationToString(lOr)<<" in outline.")

            // Log vertex coordinates for debug
            LzLogN("", "*** Faulty outline: N= "<<N<<" unique vertices.")
            for( size_t v=0 ; v<N ; v++ )
                LzLogM("", pTo.mVertices[ (*lpIdx)[v] ].ToString())

            // Fatal error
            LzLogException("", "Found an Undefined orientation but NO self intersections!")
        }
    }

    //--------------------------------------
    // Launch recursive meshing
    //--------------------------------------

#ifdef USE_ASYNC_RECURSION
    // Start profiling
    // spDQ->startProfiling("./json_files/testOutlinerTools",LzAsync::LogsMod::REALTIME);
    spDQ->startProfiling("./json_files/testOutlinerTools");

    {
        LzLogTimeN("", "----------- NARROW TIMER")

        // Mesh outline
        spDQ->dispatch( {"MeshOutline -- Init ["+std::to_string(lpIdx->back())+", "+std::to_string(lpIdx->front())+"]",
                         [&]{ MeshOutline_Recursive( lTo, lpIdx, lOptions ); }} );

        // Exit when no more work to be done
        spDQ->waitForIdle( /*pRestartPool=*/false );

        // spDQ->ExportLogsJson("./json_files/testOutlinerTools.json");
    }

    // // Log activities
    // spDQ->logActivities("MeshOutline -- activities");

    // // Check errors
    // vector<string> lDQErrs = spDQ->getLogEntries();
    // const size_t lNbErrs = LzAsync::DispatchQ::GetNbErrors( lDQErrs );
    // if( lNbErrs != 0 )
    // {
    //     // Log
    //     LzLogN("", "Found "<<lNbErrs<<" error(s) during meshing of outline:")
    //     for( const string & lErr : lDQErrs )
    //         LzLogE("", lErr)

    //     // Fail
    //     LzLogException("", "Found "<<lNbErrs<<" error(s) during meshing of outline!")
    // }
#else
    // Mesh outline
    lOptions.mJobName = "MeshOutline -- Init ["+std::to_string(lpIdx->back())+", "+std::to_string(lpIdx->front())+"]";
    MeshOutline_Recursive( pTo, lpIdx, lOptions );
#endif

#ifdef USE_ASYNC_RECURSION
    //--------------------------------------
    // Merge lTo into one Mesh
    //--------------------------------------
    for(auto& [iThId, iTo] : lTo)
    {
        pTo.Merge(iTo);
    }

#endif


    //--------------------------------------
    // Adjust mesh orientation if undefined
    //--------------------------------------

    if( lOrient == LzTriModel::Outliner::Orientation::Undefined )
    {
        LzLogN("", "Adjusting mesh orientation due to undefined orientation of outline -- probably due to self-intersections.")

        // Add-up all CCW normals weighted by triangle areas (x2) and check against prescribed normal
        Vector3D lTriN{0, 0, 0};
        for( size_t t=0 ; t<pTo.mTriangles.size() ; t++ )
        {
            // Read triangle
            const LzTriModel::Triangle & lTri = pTo.mTriangles[ t ];

            // Read vertices
            const Point3D & A = pTo.mVertices[ lTri.mIdxV[0] ];
            const Point3D & B = pTo.mVertices[ lTri.mIdxV[1] ];
            const Point3D & C = pTo.mVertices[ lTri.mIdxV[2] ];

            // Compute and add tri normal
            lTriN += (B - A)^(C - A);
        }

        // Compare to prescribed, and flip tris orientation if misoriented
        if( lTriN*pNormalCCW < 0.0 )
        {
            // Log
            LzLogM("", "Swapping triangles' orientation.")

            // Swap
            pTo.SwapTrianglesOrientation();
        }
        else
        {
            // Log
            LzLogM("", "Leaving triangles' orientation unchanged.")
        }
    }
}


//================================================================================
void Test_MinSegToSeg()
{
    // Overlapping RIGHT
    {
        const Point3D A{  0,  0, 0 };
        const Point3D B{ 10,  0, 0 };
        const Point3D C{  5, 10, 0 };
        const Point3D D{ 45, 10, 0 };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
        LzLogM("", "")
    }

    // Overlapping LEFT
    {
        const Point3D A{   0,  0, 0 };
        const Point3D B{ -10,  0, 0 };
        const Point3D C{   0, 10, 0 };
        const Point3D D{ -45, 10, 0 };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
        LzLogM("", "")
    }


    // Touching RIGHT
    {
        const Point3D A{  0,  0, 0 };
        const Point3D B{ 20, 20, 0 };
        const Point3D C{ 20, 20, 0 };
        const Point3D D{ 50, 50, 0 };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
        LzLogM("", "")
    }

    // Touching LEFT
    {
        const Point3D A{    0,   0, 0 };
        const Point3D B{  -20,   0, 0 };
        const Point3D C{  -20, -20, 0 };
        const Point3D D{ -500, -20, 0 };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
        LzLogM("", "")
    }

    // Nested 1
    {
        const Point3D A{   0,   0, 0 };
        const Point3D B{ 100,   0, 0 };
        const Point3D C{  10, -20, 0 };
        const Point3D D{  70, -20, 0 };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
        LzLogM("", "")
    }

    // Nested 2
    {
        const Point3D A{  0,  0, 0 };
        const Point3D B{ 10,  0, 0 };
        const Point3D C{  0, 20, 0 };
        const Point3D D{ 70, 20, 0 };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
        LzLogM("", "")
    }

    //--------------------------------- UNIT TEST -- BUG INTER

    {
        const Point3D A{ 76.9918, -139.958899478519, -523.811934472356 };
        const Point3D B{ 76.9918, -104.854110624386, -534.735270672913 };
        const Point3D C{ 76.9918, -109.643644867059, -533.723023410258 };
        const Point3D D{ 76.9918, -109.860178736663, -532.94760937819 };

        Point3D I, J;

        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(A, B, C, D, I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())
    }





}

//================================================================================
void Test_HasSelfIntersections()
{
    // Debug self inter: MUST HAVE SELF INTER
    {
        const List<Point3D> lOut =
        {
            { 0, -106.812438089509, -522.338343394609 },
            { 0, -106.831450289573, -522.359780424974 },
            { 0, -106.830421471258, -522.367470636691 },
            { 0, -106.833724173759, -522.339382027619 },
        };

        Point3D I, J;
        LzLogM("", "MinSegToSeg= "<<MinSegToSeg(lOut[0], lOut[1], lOut[2], lOut[3], I, J))

        LzLogM("", "I= "<<I.ToString())
        LzLogM("", "J= "<<J.ToString())

        LzLogM("", "Has self intersections= "<<std::boolalpha<<HasSelfIntersections(lOut, 1e-6))
        LzLogM("", "")
//return; //******************************
    }

    // Yes self inters 1
    {
        const List<Point3D> lOut
        {
            {  0,  0, 0 },
            { 10,  0, 0 },
            { 10, 10, 0 },
            {  5, -5, 0 },
        };

        LzLogM("", "Has self intersections= "<<std::boolalpha<<HasSelfIntersections(lOut, 1e-6))
        LzLogM("", "")
    }

    // Yes self inters 2
    {
        const List<Point3D> lOut
        {
            {   0,   0, 0 },
            {  10,   0, 0 },
            {  10,  10, 0 },
            { -10, -10, 0 },
            {   0,   0, 0 },
        };

        LzLogM("", "Has self intersections= "<<std::boolalpha<<HasSelfIntersections(lOut, 1e-6))
        LzLogM("", "")
    }

    // Yes self inters 3
    {
        const List<Point3D> lOut
        {
            {   0,   0, 0 },
            {  10,   0, 0 },
            {  10,  10, 0 },
            { -10, -10, 0 },
        };

        LzLogM("", "Has self intersections= "<<std::boolalpha<<HasSelfIntersections(lOut, 1e-6))
        LzLogM("", "")
    }

    // No self inters 1
    {
        const List<Point3D> lOut
        {
            {  0,  0, 0 },
            {  5,  0, 0 },
            { 10,  0, 0 },
            { 10,  5, 0 },
            { 10, 10, 0 },
            {  5,  5, 0 },
            {  0,  0, 0 },
        };

        LzLogM("", "Has self intersections= "<<std::boolalpha<<HasSelfIntersections(lOut, 1e-6))
        LzLogM("", "")
    }
}

//================================================================================
void Test_FlipFlopper()
{
    // Fan iterator
    class FlipFlopper
    {
    public:
        FlipFlopper( int pMin, int pMax )
        {
            mMin = pMin;
            mMax = pMax;

            // Integer division
    //        mI0 = static_cast<int>( (pMin + pMax) / 2.0 );
            mI0 = (pMin + pMax) / 2;
            mDelta = 0;
            mSgn = Sgn::Neg;
        }

        int Init() const { return mI0; }

        bool Next( int & pNext )
        {
            switch( mSgn )
            {
            case Sgn::Pos:
                pNext += mDelta++;
                mSgn = Sgn::Neg;
                return pNext <= mMax;

            case Sgn::Neg:
                pNext -= mDelta++;
                mSgn = Sgn::Pos;
                return pNext >= mMin;
            }
        }

    protected:
        enum class Sgn { Pos, Neg };
        Sgn mSgn;

        int mMin;
        int mMax;

        int mI0;
        int mDelta;
    };

    {
        // Pair to Pair
        LzLogM("", "Even --> Even: -2 --> 4\n")
        FlipFlopper lFF(-2, 4);
        for( int x=lFF.Init() ; lFF.Next(x) ;)
            LzLogM("", x)
    }
    LzLogM("", "-------------------------------")
    {
        // Impair to Pair
        LzLogM("", "Odd --> Even: 5 --> 14\n")
        FlipFlopper lFF(5, 14);
        for( int x=lFF.Init() ; lFF.Next(x) ;)
            LzLogM("", x)
    }
    LzLogM("", "-------------------------------")
    {
        // Pair to Impair
        LzLogM("", "Even --> Odd: 2 --> 7\n")
        FlipFlopper lFF(2, 7);
        for( int x=lFF.Init() ; lFF.Next(x) ;)
            LzLogM("", x)
    }
    LzLogM("", "-------------------------------")
    {
        // Impair to Impair
        LzLogM("", "Odd --> Odd: 5 --> 13\n")
        FlipFlopper lFF(5, 13);
        for( int x=lFF.Init() ; lFF.Next(x) ;)
            LzLogM("", x)
    }
    LzLogM("", "-------------------------------")
    {
        // Self to self: Impair
        LzLogM("", "Same odd --> Same odd: 1 --> 1\n")
        FlipFlopper lFF(1, 1);
        for( int x=lFF.Init() ; lFF.Next(x) ;)
            LzLogM("", x)
    }
    LzLogM("", "-------------------------------")
    {
        // Self to self: Pair
        LzLogM("", "Same even --> Same even: 8 --> 8\n")
        FlipFlopper lFF(8, 8);
        for( int x=lFF.Init() ; lFF.Next(x) ;)
            LzLogM("", x)
    }
}

//================================================================================
void Test_MeshOutline()
{
#if 1
{
    static constexpr size_t N_points = 100;
    static constexpr double R = 10.00;

    // Generate points
    List<Point3D> lPts;
    {
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
        //srerttzete
    }

    // Mesh
    Mesh lOutMesh;
    LzTriModel::OutlinerTools::MeshOutline( lOutMesh, lPts, Vector3D{0,0,1} );
}
return;
#endif

#if 0
{
    class A
    {
    public:
        A() { LzLogM("", "A::A") }
        A( const A & ) { LzLogM("", "A::A(const A &)") }
        A( A && ) { LzLogM("", "A::A(A &&)") }
        void operator=( A && ) { LzLogM("", "A::operator=(A &&)") }
        ~A() { LzLogM("", "A::~A") }
    };

    std::map<string, A> lMap;

//    A a;
    lMap["test"] = A();

    LzLogM("", "-----")
    lMap.erase( lMap.begin() );

    LzLogM("", "-----")
    return;
}
#endif



#if 0
//--------------------------------- UNIT TEST -- BUG FEMUR
{
    const List<Point3D> lPts_BUG_FEM =
    {
        { 131.173805237234, -114.602790605832, -528.995 },
        { 130.870205399918, -115.550445679137, -528.995 },
        { 114.155702585816, -126.286538994374, -528.995 },
        { 115.189287108929, -160.239815850542, -528.995 },
        { 135.84692903162, -109.66567411534, -528.995 },
        { 131.173805237234, -114.602790605832, -528.995 },
    };
    //
    const Vector3D lNormalCCW{ 0, 0, 1 };

    Mesh lTriMesh;
    MeshOutline( lTriMesh, lPts_BUG_FEM, lNormalCCW );
    lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__BUG_FEM.obj");
}
#endif


#if 0  // UNIT TEST 1 TRI CCW
    {
        const List<Point3D> lPts_1_tri__CCW =
        {
            { 0, 0, 0 },
            { 100, 0, 0 },
            { 0, 100, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW_one_tri{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_1_tri__CCW, lNormalCCW_one_tri );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__1_tri__CCW.obj");

        return;
    }
#endif


#if 0  // UNIT TEST 1 TRI CW
    {
        const List<Point3D> lPts_1_tri__CW =
        {
            { 0, 0, 0 },
            { 0, 100, 0 },
            { 100, 0, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW_one_tri{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_1_tri__CW, lNormalCCW_one_tri );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__1_tri__CW.obj");

        return;
    }
#endif


#if 0  // UNIT TEST 2 TRIS CCW
    {
        const List<Point3D> lPts_2_tris__CCW =
        {
            { 0, 0, 0 },
            { 100, 0, 0 },
            { 100, 100, 0 },
            { 0, 100, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW_one_tri{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_2_tris__CCW, lNormalCCW_one_tri );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__2_tris__CCW.obj");

        return;
    }
#endif


#if 0  // UNIT TEST 2 TRIS CW
    {
        const List<Point3D> lPts_2_tris__CW =
        {
            { 0, 0, 0 },
            { 0, 100, 0 },
            { 100, 100, 0 },
            { 100, 0, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW_one_tri{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_2_tris__CW, lNormalCCW_one_tri );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__2_tris__CW.obj");

        return;
    }
#endif


#if 0  // UNIT TEST 3 TRIS CCW
    {
        const List<Point3D> lPts_3_tris__CCW =
        {
            { 0, 0, 0 },
            { 100, 0, 0 },
            { 100, 100, 0 },
            { 50, 150, 0 },
            { 0, 100, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_3_tris__CCW, lNormalCCW );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__3_tris__CCW.obj");

        return;
    }
#endif


#if 0  // UNIT TEST 3 TRIS CW
    {
        const List<Point3D> lPts_3_tris__CW =
        {
            { 0, 0, 0 },
            { 0, 100, 0 },
            { 50, 150, 0 },
            { 100, 100, 0 },
            { 100, 0, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_3_tris__CW, lNormalCCW );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__3_tris__CW.obj");

        return;
    }
#endif




#if 0  // UNIT TEST SELF INTERSECT 3 TRIS CCW
    {
        const List<Point3D> lPts_self_isect_3_tris__CW =
        {
            { 0, 0, 0 },
            { 100, 0, 0 },
            { 100, 80, 0 },
            { -20, 80, 0 },
            { 0, 100, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_self_isect_3_tris__CW, lNormalCCW );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__self_isect_3_tris__CCW.obj");

        // Compare with OLD
        LzLogInfoMessage("Trying old meshing code...", "")
        LzTriModel::Outliner::MeshOutline( lPts_self_isect_3_tris__CW, lNormalCCW, lTriMesh );

        return;
    }
#endif


#if 0  // UNIT TEST SELF INTERSECT 3 TRIS CW
    {
        const List<Point3D> lPts_self_isect_3_tris__CW =
        {
            { 0, 0, 0 },
            { 0, 100, 0 },
            { -20, 80, 0 },
            { 100, 80, 0 },
            { 100, 0, 0 },
            { 0, 0, 0 },
        };
        //
        const Vector3D lNormalCCW{ 0, 0, 1 };

        Mesh lTriMesh;
        MeshOutline( lTriMesh, lPts_self_isect_3_tris__CW, lNormalCCW );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__self_isect_3_tris__CW.obj");

        // Compare with OLD
        LzLogInfoMessage("Trying old meshing code...", "")
        LzTriModel::Outliner::MeshOutline( lPts_self_isect_3_tris__CW, lNormalCCW, lTriMesh );

        return;
    }
#endif


#if 0  // UNIT TEST FEMUR POST COND WITH INVERTED TRIS CW
    {
        // Outline points
        const double lPts__CW[181][3] =
        {
            //-- Out[1]: 181 points.
            { -40.6339026477, 35.5232803296, -230.792706406 },
            { -40.3995726378, 35.4766658168, -230.775675241 },
            { -39.9344826316, 35.3841469626, -230.74187238 },
            { -39.2594002735, 35.2498550073, -230.69280722 },
            { -37.998222744, 34.9989730308, -230.601144513 },
            { -37.5367427081, 34.9071722963, -230.567604025 },
            { -37.4132025366, 34.8825968486, -230.558625094 },
            { -37.308989735, 34.8618661323, -230.55105088 },
            { -37.0268410756, 34.8057392099, -230.530544243 },
            { -35.5323779231, 34.5084504771, -230.421926277 },
            { -34.8305185167, 34.3688318501, -230.370914955 },
            { -34.4308019343, 34.2893175201, -230.341863451 },
            { -33.0141500208, 34.007507526, -230.238900825 },
            { -32.6541720358, 33.935898267, -230.212737532 },
            { -31.2479499958, 33.656163054, -230.110532953 },
            { -31.2002634994, 33.6466769332, -230.107067086 },
            { -31.0653754119, 33.6198440812, -230.097263385 },
            { -29.8775943818, 33.3835626334, -230.010935154 },
            { -29.3686436861, 33.2823187138, -229.973944487 },
            { -28.6710770003, 33.1435540238, -229.923245161 },
            { -26.4443097423, 32.7005903982, -229.761403143 },
            { -26.0055454878, 32.6133084408, -229.729513644 },
            { -25.5131243401, 32.5153526908, -229.693724349 },
            { -25.4043861808, 32.4937217597, -229.685821231 },
            { -25.2190603215, 32.4568554845, -229.6723517 },
            { -24.9497669865, 32.4032858302, -229.652779391 },
            { -24.737414651, 32.3610432653, -229.637345569 },
            { -24.3161472701, 32.2772419046, -229.606727747 },
            { -23.6824373165, 32.1511800283, -229.560669545 },
            { -23.6691619284, 32.1485391981, -229.559704686 },
            { -22.9107987751, 31.9976804631, -229.504586657 },
            { -22.7671859241, 31.969112022, -229.494148838 },
            { -21.988425883, 31.8141957998, -229.437548358 },
            { -21.7458814011, 31.7659472087, -229.419920163 },
            { -20.3196521482, 31.4822320254, -229.316261454 },
            { -20.1615090067, 31.4507731205, -229.304767569 },
            { -19.3377259287, 31.2869006109, -229.244894803 },
            { -19.2065268186, 31.2608015953, -229.235359218 },
            { -18.7974415765, 31.179423588, -229.205626798 },
            { -18.4690843271, 31.1141045399, -229.181761708 },
            { -17.0466078923, 30.8311358926, -229.078375755 },
            { -16.7748367363, 30.7770733335, -229.058623358 },
            { -16.2222379185, 30.6671466337, -229.018460333 },
            { -15.6608345785, 30.5554684788, -228.977657394 },
            { -15.3190085516, 30.4874701302, -228.952813391 },
            { -9.045057888, 29.2394133688, -228.496821042 },
            { -8.73818604423, 29.1783683431, -228.474517518 },
            { -8.2887897982, 29.0889713979, -228.441855283 },
            { -7.87705230548, 29.007065787, -228.411930096 },
            { -7.35347058008, 28.9029113638, -228.373876041 },
            { -6.84652005525, 28.8020653317, -228.337030747 },
            { -5.8964118333, 28.6130633688, -228.267976637 },
            { -5.49440752823, 28.5330939495, -228.238758861 },
            { -5.44884323814, 28.5240299923, -228.235447236 },
            { -5.51577924756, 28.5433835292, -227.194911488 },
            { -5.50452399889, 28.5430910647, -226.857091072 },
            { -5.61632994656, 28.5725536598, -225.614962211 },
            { -5.64718224145, 28.5820799207, -225.030475324 },
            { -5.26699627332, 28.5157681623, -223.389701891 },
            { -5.25703950697, 28.5139582213, -223.359420024 },
            { -5.20927290993, 28.5056017448, -223.157612088 },
            { -4.96814805434, 28.4665336517, -221.599542141 },
            { -4.99781745534, 28.4737434858, -221.375277255 },
            { -4.86155899584, 28.4497204573, -220.831707269 },
            { -4.56532975117, 28.3963509473, -219.847837464 },
            { -4.45185054094, 28.3750559638, -219.618141927 },
            { -4.46566048791, 28.37737176, -219.693829911 },
            { -4.4023312915, 28.3670678049, -219.292076671 },
            { -4.4065427062, 28.3717525119, -218.626353793 },
            { -4.29033273108, 28.3507877067, -218.245246243 },
            { -4.26296810846, 28.3507037772, -217.315334678 },
            { -4.0064776911, 28.3001511655, -217.21528634 },
            { -3.96990563026, 28.2974949293, -216.412944556 },
            { -3.77260824519, 28.2633047216, -215.522985786 },
            { -3.76897976495, 28.2647016085, -215.15590931 },
            { -3.50775587505, 28.2200960311, -213.862871657 },
            { -3.49522161025, 28.2190909012, -213.604293386 },
            { -3.40189141254, 28.2036655233, -213.053789677 },
            { -3.29448888925, 28.18778626, -212.09618615 },
            { -3.25681853956, 28.1828429057, -211.651912092 },
            { -3.17656867351, 28.1730605157, -210.575871721 },
            { -2.9794073071, 28.1400233589, -209.490977423 },
            { -2.95833602461, 28.1379433448, -209.123855631 },
            { -2.87293821157, 28.1274275511, -207.997118426 },
            { -2.86343413125, 28.1267900306, -207.779476781 },
            { -2.86661722122, 28.1277531149, -207.72259492 },
            { -3.10962070998, 28.1844310884, -206.296669406 },
            { -3.1575998664, 28.1959106414, -205.965109099 },
            { -3.35522649086, 28.2420214309, -204.802597249 },
            { -3.41236660437, 28.2548097749, -204.560617924 },
            { -3.48515391156, 28.2704346517, -204.367579806 },
            { -3.93846019026, 28.3666669875, -203.351765406 },
            { -4.09864745528, 28.4007323769, -202.982542871 },
            { -4.85415729891, 28.5585674498, -201.7313568 },
            { -4.96855307719, 28.5825775777, -201.52260779 },
            { -5.81680693798, 28.7582251453, -200.3884082 },
            { -5.94760199191, 28.7851847989, -200.234994545 },
            { -5.98460749451, 28.7927364622, -200.2047409 },
            { -7.27217552794, 29.054896971, -199.254551557 },
            { -7.72561388753, 29.1468622694, -198.98205629 },
            { -8.83618510601, 29.3720063662, -198.331853394 },
            { -9.22670839596, 29.4509404501, -198.144073528 },
            { -10.3650011121, 29.6809097331, -197.615236014 },
            { -10.8461443201, 29.7779898716, -197.413383145 },
            { -11.8330170354, 29.9770841185, -197.004002941 },
            { -12.6206959334, 30.1361120689, -196.656550693 },
            { -12.9319141027, 30.1987317333, -196.556239793 },
            { -14.1803791258, 30.4490831409, -196.301010858 },
            { -14.3349704819, 30.4800016906, -196.283473799 },
            { -14.4897064523, 30.5109640943, -196.26333531 },
            { -15.6816265537, 30.748551123, -196.266478589 },
            { -15.9735815944, 30.8063189302, -196.341345833 },
            { -16.8487092341, 30.9800449746, -196.467325159 },
            { -17.3820070434, 31.0856802549, -196.584329274 },
            { -17.7027616802, 31.149308653, -196.638508947 },
            { -18.6322529295, 31.3333814032, -196.849408338 },
            { -19.3064685588, 31.466517864, -197.068612868 },
            { -19.8647436092, 31.5765157422, -197.29235795 },
            { -20.5722916205, 31.7157939095, -197.598646041 },
            { -21.0843499291, 31.8163202256, -197.867130701 },
            { -21.5825012469, 31.9143422203, -198.089217797 },
            { -22.2765480463, 32.0505243225, -198.465562844 },
            { -23.0057478031, 32.1942949451, -198.741335735 },
            { -23.3288401503, 32.2577758868, -198.901726837 },
            { -24.5248524253, 32.492331517, -199.570936375 },
            { -24.527636336, 32.4928795176, -199.572141808 },
            { -24.5296087476, 32.4932594588, -199.574436204 },
            { -25.7430949447, 32.7292566947, -200.597139784 },
            { -25.9176112573, 32.7629946953, -200.779140917 },
            { -26.8763728969, 32.9490128973, -201.663437099 },
            { -27.0116597109, 32.9750517433, -201.824469973 },
            { -27.7624201256, 33.1210839261, -202.452823157 },
            { -28.3302949394, 33.2313431271, -202.962644074 },
            { -28.4778543694, 33.2602066857, -203.05819693 },
            { -29.0907218732, 33.3798503945, -203.496104114 },
            { -29.6871670674, 33.4959766902, -203.976200464 },
            { -29.8793534239, 33.5332563077, -204.15488537 },
            { -30.5315718842, 33.6593485738, -204.834471087 },
            { -30.9296104609, 33.7364208691, -205.228403411 },
            { -31.6133785535, 33.8687681219, -205.913931762 },
            { -32.091451608, 33.9613216192, -206.389824181 },
            { -32.599905259, 34.0592493286, -206.983802217 },
            { -33.2336171069, 34.1825016246, -207.51635394 },
            { -33.2823909164, 34.1919372578, -207.566087567 },
            { -33.8014845916, 34.2918132549, -208.189976678 },
            { -34.716883231, 34.468123979, -209.258348931 },
            { -34.7271855865, 34.4700997797, -209.271842594 },
            { -34.7345299084, 34.4715043445, -209.282144051 },
            { -35.5382306447, 34.6250342373, -210.439547242 },
            { -36.1227983966, 34.7364248802, -211.329607156 },
            { -36.3337100542, 34.7767179192, -211.632843743 },
            { -36.6163608361, 34.830837099, -212.01828468 },
            { -37.112628561, 34.9244203777, -212.943846461 },
            { -37.4311766318, 34.9843411456, -213.563779648 },
            { -37.8071027173, 35.0547886664, -214.341497965 },
            { -38.2477618982, 35.1363320785, -215.432313176 },
            { -38.3286202074, 35.1509275934, -215.696047144 },
            { -38.7049470423, 35.2202855613, -216.676233043 },
            { -38.755197416, 35.2295201895, -216.811728104 },
            { -38.7716332167, 35.2325825575, -216.848787443 },
            { -39.1480767415, 35.3019714371, -217.827648891 },
            { -39.4726777454, 35.3603324152, -218.926538037 },
            { -39.5725458495, 35.3783807416, -219.248570008 },
            { -39.6963747502, 35.4005577983, -219.682752073 },
            { -39.9097532092, 35.437353031, -220.676708092 },
            { -40.1292254261, 35.4765548365, -221.464324917 },
            { -40.1944642871, 35.487054658, -221.898076184 },
            { -40.3939268178, 35.5194425126, -223.174806148 },
            { -40.3998508369, 35.5202102931, -223.246336029 },
            { -40.5178636977, 35.5369645239, -224.418652089 },
            { -40.5747460185, 35.5454079084, -224.920026848 },
            { -40.5809283699, 35.5423427681, -225.664073745 },
            { -40.6889348396, 35.5621801452, -225.957245406 },
            { -40.9170286379, 35.5998733321, -227.30361124 },
            { -40.9368968277, 35.6027767113, -227.486659554 },
            { -40.9353420113, 35.6021331786, -227.544413851 },
            { -40.9050933638, 35.5876304751, -229.011315876 },
            { -40.7761611113, 35.5567804163, -229.902574388 },
            { -40.6987970441, 35.5385110206, -230.395504342 },
            { -40.6559442875, 35.5283914134, -230.668543339 },
            { -40.6339026477, 35.5232803296, -230.792706406 },
        };
        //
        const Vector3D lNormalCCW{ 0.195496785964395, 0.980687983661217, -0.00566439582728056 };

        // Build outline
        const size_t N = sizeof(lPts__CW) / sizeof(*lPts__CW);
        LzLogM("", "N= "<<N)
        //
        List<Point3D> lOut;
        for( int p=0 ; p<N ; p++ )
            lOut.AddTail( {lPts__CW[p]} );

        // Meshing options
        MeshOutlineOptions lOptions;
//lOptions.mEarlyStopEvaluation = true;
////lOptions.mEarlyStopMinDist = 1.0;
//lOptions.mEarlyStopMinDist = 0.5;

        // Generate mesh
        Mesh lTriMesh;
        MeshOutline( lTriMesh, lOut, lNormalCCW, lOptions );
        lTriMesh.Save(LzServices::StartUpPath()+"/__DEBUG__TriMeshOutline__Femur_post_cond_invert_tris__CW.obj");

        // Compare with OLD
        LzLogInfoMessage("Trying old meshing code...", "")
        LzTriModel::Outliner::MeshOutline( lOut, lNormalCCW, lTriMesh );

        return;
    }
#endif


#if 0  // UNIT TEST GAZILLION CUTS
    {
        // Load test mesh
        Mesh lFem;
#if 0
        // Normal bone
        lFem.Load(R""(D:\LinkZ\Data\TwInsight\ProstFlex\SurgiTwin_2\resources\L_femur.obj)"");
#else
        // Buggy mesh
        lFem.Load(R""(D:\TwInsight\Dev\surgitwin\test\data\planning\subject_test_capping\resources\L_femur.obj)"");
#endif

        // Set outliner
        LzTriModel::Outliner lOut;
        lOut.Set( lFem );

        // Compute BBox
        const LzGeom::BBox lBBox( lFem.mVertices );

        // Log
        LzLogTimeN("", "------ Computing and meshing a gazillion cuts")

        // Select slicing step
//const double STEP_mm = 2.5;
const double STEP_mm = 1.0;

        // Find widest dim, min and max along this dim
//int dim=0; // *** Found an Undefined orientation but NO self intersections!
//int dim=2;
for( int dim=2 ; dim>=0 ; dim-- )
        {
            // Get Min / Max
            const double lMin = lBBox.Min( dim );
            const double lMax = lBBox.Max( dim );

            // Cut plane equation
            double lCutEq[] = { 0, 0, 0, 0 };
            lCutEq[dim] = -1;

            // Choose a CCW normal
            Vector3D lNorCCW{0, 0, 0};
            lNorCCW.mV[dim] = 1;

            // Slice real nice
            //
            size_t iFile = 0;
//int dim=2; double x=-522.791922485; // inner test FAILS when self inters and leads to huge mis-meshing !!!!
for( double x=lMin ; x<=lMax ; x+= STEP_mm )
            {
                LzLogN("", "Cutting along dim "<<dim<<": x= "<<x)

                // Prepare cut plane
                lCutEq[3] = x;
                Plane3D lCut( lCutEq );

                // Extract outlines
                List<List<Point3D>> lOuts;
                lOut.Cut( lCut, lOuts );

                // Mesh all
                BrowseList( iO, lOuts )
                {
                    // Get outline
                    /*const*/ List<Point3D> & lO = lOuts.GetAt( iO );

////Corrected vertex #352: [ 111.382186839, -106.461660004, -522.791922485
//lO[352] = {111.382186839, -106.461660004, -522.791922485};

                    // Meshing options
                    MeshOutlineOptions lOptions;
//lOptions.mEarlyStopEvaluation = false;

                    // Mesh outline
                    Mesh lMesh;
#if 1
                    // NEW
                    MeshOutline( lMesh, lO, lNorCCW, lOptions );
#elif 0
                    // Shewchuk
                    LzTriModel::Outliner::MeshOutline_Shewchuk( lO, lNorCCW, 10.0, lMesh );

                    Crashes very quickly
#else
                    // OLD
                    LzTriModel::Outliner::MeshOutline( lO, lNorCCW, lMesh );
                    // Fails on Femur -- bad mesh
#endif
                    // Save
                    lMesh.Save(LzServices::StartUpPath()+"/__SLICE_FEMUR/__DEBUG__MeshOutline__Slice_Femur__dim="+std::to_string(dim)+"__"+std::to_string(iFile++)+".obj");
                }
            }
        }

        return;
    }
#endif













}


}
}


