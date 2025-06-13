#include "OutlinerTools.h"
#include <LzGeom/Line3D.h>
#include <LzGeom/RigidTr3D.h>
#include <LzTriModel/Mesh.h>
#include <LzTriModel/MeshTools.h>


namespace LzTriModel
{
namespace OutlinerTools
{
//================================================================================
bool HasSelfIntersections( const List<Point3D> & pOut, double pEpsilon )
{
    // Check closedness
    const bool lIsClosed = Outliner::IsClosed( pOut );

LzLogException("", "*** Recuperer depuis Inferator / mainwindow.cpp !!!!")


    return true;
}

//================================================================================
bool ProjectOnOutline( const Point3D & pPt,
                       const List<Point3D> & pOut,
                       const Vector3D & pDir,
                       Point3D & pProjPt,
                       const Vector3D * ppNormal/*=nullptr*/,
                       std::pair<void*, void*> * ppProjSegPos/*=nullptr*/)
{
    // Deal with the empty case
    if( pOut.Count() == 0 )
        return false;

    // Depending on whether we have a normal or not
    if( ppNormal )
    {
        // Check: a valid outline is either empty or has >= 2 points
        if( pOut.Count() < 2 )
            LzLogException("", "Cannot project point on outline! Outline must have at "
                               "least 2 points.")
    }
    else
    {
        // Check: this check does not ensure that the LSQ plane will be properly computed
        //        as the outline could well be a straight line (and the plane ill defined)
        if( pOut.Count() < 3 )
            LzLogException("", "Cannot project point on outline! Outline must have at "
                               "least 3 points, since no normal was provided.")
    }

    // Find a normal vector to the polygon
    const Vector3D & lNorm = ppNormal ? *ppNormal : Plane3D(pOut).Normal() ;

    // Here: either a normal was provided, or it could be computer (meaning the outline actually defines a plane)

    // Check
    if( fabs(lNorm * pDir) > 1e-6 )
        LzLogException("", "Cannot project point on outline! The projection direction is not within the outline plane.")

    // Compute intersection plane
    const Plane3D lCutPlane( pPt, pDir, lNorm );

    // Here: the cut plane could be computed, meaning the vectors are non null

    // Check all segments
    double lMinDot = -1.0;
    BrowseList( iP0, pOut )
    {
        // Last segment?
        void * iP1 = pOut.NextPos( iP0 );
        if( !iP1 )
            break;

        // Two points
        const Point3D & A = pOut.GetAt( iP0 );
        const Point3D & B = pOut.GetAt( iP1 );

        // Project on segment [A, B]
        Point3D lProjPt;
        double lDot = ProjectOnSegment( pPt, A, B, pDir, lProjPt, lCutPlane );

        // Ignore backward and invalid projections
        if( lDot < 0.0 )
            continue;

        // Deal with the case where pPt projects on itself
        if( lDot == 0.0 )
        {
            // Point projects ON segment
            pProjPt = pPt;

            if( ppProjSegPos )
            {
                ppProjSegPos->first = iP0;
                ppProjSegPos->second = iP1;
            }
            return true;
        }

        // Closer to pPt while on the positive side?
        if( lMinDot<0.0 || lMinDot>lDot )
        {
            lMinDot = lDot;
            pProjPt = lProjPt;

            if( ppProjSegPos )
            {
                ppProjSegPos->first = iP0;
                ppProjSegPos->second = iP1;
            }
        }
    }

    // Have an interserction only if min dot is non negative
    return lMinDot != -1.0;
}

//================================================================================
bool ProjectOnSegments( const Point3D & pPt,
                        const List<std::pair<Point3D, Point3D>> & pSegs,
                        const Vector3D & pDir,
                        Point3D & pProjPt,
                        const Vector3D * ppNormal/*=nullptr*/ )
{
    // Deal with the empty case
    if( pSegs.Count() == 0 )
        return false;

    // Depending on whether we have a normal or not
    if( !ppNormal )
    {
        // Check
        if( pSegs.Count() < 2 )
            LzLogException("", "Cannot project point on segments! Need at least 2 "
                               "segments (at least 3 points), since no normal was provided.")
    }

    // Find a normal vector to the polygon
    const Vector3D & lNorm = ppNormal ? *ppNormal : Plane3D( [&](Plane3D * pThis)
        {
            // Collect all points
            List<Point3D> lAllPts;
            BrowseList( iS, pSegs )
            {
                // Segment
                const std::pair<Point3D, Point3D> & lS = pSegs.GetAt(iS);

                // Two points
                lAllPts.AddTail( lS.first );
                lAllPts.AddTail( lS.second );
            }

            // Compute the temporary LSQ plane
            *pThis = Plane3D(lAllPts);
        } ).Normal() ;

    // Here: either a normal was provided, or it could be computer (meaning the outline actually defines a plane)

    // Check
    if( fabs(lNorm * pDir) > 1e-6 )
        LzLogException("", "Cannot project point on segments! The projection direction is not within the segments plane.")

    // Compute intersection plane
    const Plane3D lCutPlane( pPt, pDir, lNorm );

    // Here: the cut plane could be computed, meaning the vectors are non null

    // Check all segments
    double lMinDot = -1.0;
    BrowseList( iS, pSegs )
    {
        // Segment
        const std::pair<Point3D, Point3D> & lS = pSegs.GetAt(iS);

        // Two points
        const Point3D & A = lS.first;
        const Point3D & B = lS.second;

        // Project on segment [A, B]
        Point3D lProjPt;
        double lDot = ProjectOnSegment( pPt, A, B, pDir, lProjPt, lCutPlane );

        // Ignore backward and invalid projections
        if( lDot < 0.0 )
            continue;

        // Deal with the case where pPt projects on itself
        if( lDot == 0.0 )
        {
            // Point projects ON segment
            pProjPt = pPt;
            return true;
        }

        // Closer to pPt while on the positive side?
        if( lMinDot<0.0 || lMinDot>lDot )
        {
            lMinDot = lDot;
            pProjPt = lProjPt;
        }
    }

    // Have an interserction only if min dot is non negative
    return lMinDot != -1.0;
}

//================================================================================
double ProjectOnSegment( const Point3D & pPt,
                         const Point3D & pA,
                         const Point3D & pB,
                         const Vector3D & pDir,
                         Point3D & pProjPt,
                         const Plane3D & pCutPlane )
{
    // One line
    const Line3D AB( pA, pB );

    // Get intersection
    Point3D I;
    if( pCutPlane.Intersection_NoExc(AB, I) )
    {
        // Here: we have a proper intersection
        // -----

        // Reject if intersection is not on the segment
        if( (I - pA)*(I - pB) > 0.0 )
            return -1.0;

        // Pointing in the right direction
        const double lDot = pDir * (I - pPt);

        // Check dot product
        if( lDot < 0.0 )
        {
            // Ignore intersection if on the wrong side
            return -1.0;
        }
        else
        {
            // Yes, valid intersection
            pProjPt = I;
            return lDot;
        }
    }
    else
    {
        // Here: no proper intersection, the line is parallel to
        // ----- (or included in) the cut plane

        // Check if point is on the line (AB)
        if( AB.DistanceTo(pPt) < 1e-6 )
        {
            // Point is in [A, B]
            if( (pPt - pA)*(pPt - pB) <= 0.0 )
            {
                // Yes, point is on [A, B] and projects onto self
                pProjPt = pPt;
                return 0.0;
            }
            else
            {
                // Point is outside [A, B]

                // Check dot product with A
                const double lDot_A = pDir * (pA - pPt);

                // In all cases, this must be positive
                if( lDot_A < 0.0 )
                    return -1.0;

                // pPt projects 'axially' onto A or B, pick
                // the one with the smallest dot product
                const double lDot_B = pDir * (pB - pPt);
                //
                if( lDot_A < lDot_B )
                {
                    pProjPt = pA;
                    return lDot_A;
                }
                else
                {
                    pProjPt = pB;
                    return lDot_B;
                }
            }
        }
        else
        {
            // Nope, point is not on the line (AB)
            return -1.0;
        }
    }
}

//================================================================================
void SegmentsToMesh( const List<std::pair<Point3D, Point3D>> & pSegs, Mesh & pTo, double pScale/*=1.0*/ )
{
    // Free previous
    pTo.Free();

    // Scale model
    const double S = pScale * 0.5;

    // Create cube
    const Mesh lCube( [&]( Mesh * pThis )
        {
            BBox lBBox(vector<Point3D>{ Point3D{-S, -S, -S}, Point3D{+S, +S, +S} });
            pThis->FromBBox( lBBox );
        } );

    // Cube merger
    auto MergeCubeAt = [&]( const Point3D & pPt )
        {
            Mesh C = lCube;
            C.RigidTransform( {0,0,0, pPt - Point3D(0,0,0)} );
            pTo.Merge( C );
        };

    // Generate mesh
    BrowseList( iS, pSegs )
    {
        // Segment
        const std::pair<Point3D, Point3D> & lS = pSegs.GetAt(iS);

        // Two points
        const Point3D & A = lS.first;
        const Point3D & B = lS.second;

        // Does the segment deserve a cylinder?
        if( A.DistanceTo(B) > 1e-6 )
        {
            // Merge cylinder
            Mesh lCyl;
            LzTriModel::GenCylinder(lCyl, A, S/2, B, S/2, 10, 20 );
            pTo.Merge( lCyl );

            // Merge cubes
            MergeCubeAt( A );
            MergeCubeAt( B );
        }
        else
        {
            // One single cube
            MergeCubeAt( A );
        }
    }
}

//================================================================================
void OutlineToMesh( const List<Point3D> & pOut, Mesh & pTo, double pScale/*=1.0*/ )
{
    // Free previous
    pTo.Free();

    // Check
    if( pOut.Count() < 2 )
        LzLogException("", "Cannot convert outline to mesh! Invalid outline.");

    // Here: the outline has at least 2 points

    // Scale model
    const double S = pScale * 0.5;

    // Create cube
    const Mesh lCube( [&]( Mesh * pThis )
        {
            BBox lBBox(vector<Point3D>{ Point3D{-S, -S, -S}, Point3D{+S, +S, +S} });
            pThis->FromBBox( lBBox );
        } );

    // Create cross
    const Mesh lCross( [&]( Mesh * pThis )
        {
            BBox lBBox(vector<Point3D>{ Point3D{-3*S, -S/2, -S/2}, Point3D{+3*S, +S/2, +S/2} });
            Mesh lBar;
            lBar.FromBBox(lBBox);
            lBar.RigidTransform( RigidTr3D(0,0,45, {0,0,0}) );
            *pThis = lBar;
            lBar.RigidTransform( RigidTr3D(0,0,90, {0,0,0}) );
            pThis->Merge(lBar);
        } );

    // Draw segments (in a closed outline the last point repeats first so no need to go around the list)
    BrowseList( iP0, pOut )
    {
        // Last segment?
        void * iP1 = pOut.NextPos( iP0 );
        if( !iP1 )
            break;

        // Two points
        const Point3D & lV = pOut.GetAt( iP0 );
        const Point3D & lW = pOut.GetAt( iP1 );

        // Check for null segment
        if( lV.DistanceTo(lW) < 1e-6 )
            continue;

        // Merge cylinder
        Mesh lCyl;
        LzTriModel::GenCylinder(lCyl, lV, S/2, lW, S/2, 10, 20 );
        pTo.Merge( lCyl );
    }

//    // Intermediate cubes
//    for( void * iP=pOut.NextPos(pOut.HeadPos()) ; iP!=pOut.TailPos() ; pOut.Next(iP) )
//    {
//        // Read point
//        const Vector3D lToPt = pOut.GetAt(iP) - Point3D{0,0,0};

//        // Add a black pearl
//        Mesh C = lCube;
//        C.RigidTransform( {0,0,0, lToPt} );
//        pTo.Merge( C );
//    }

    // Check whether outline is closed or not
    const bool lIsClosed = Outliner::IsClosed( pOut );

    // Draw end cubes
    if( lIsClosed )
    {
        // Closed outline: just one cross
        Mesh X = lCross;
        X.RigidTransform( {0,0,0, pOut.GetHead() - Point3D(0,0,0)} );
        pTo.Merge( X );
    }
    else
    {
        // Open outline: cubes all the way
        Mesh C = lCube;
        C.RigidTransform( {0,0,0, pOut.GetHead() - Point3D(0,0,0)} );
        pTo.Merge( C );
        //
        C = lCube;
        C.RigidTransform( {0,0,0, pOut.GetTail() - Point3D(0,0,0)} );
        pTo.Merge( C );
    }
}

//================================================================================
void * FindPosOfMaxLength( const List<List<Point3D>> & pOuts )
{
    // Result
    void * iMaxO = nullptr;
    double lMaxLen = 0.0;

    // Check all outlines
    BrowseList( iO, pOuts )
    {
        // Get length of current outline
        const double lLen = Outliner::OutlineLength( pOuts.GetAt(iO) );

        // Update max?
        if( lMaxLen < lLen )
        {
            lMaxLen = lLen;
            iMaxO = iO;
        }
    }

    // All done
    return iMaxO;
}

//================================================================================
void SplitClosedOutline( const List<Point3D> & pOut,
                         const List<void *> & pCutPos,
                         List< List<Point3D> > & pSubOuts )
{
    // Clear previous
    pSubOuts.DelAll();

    // Empty?
    if( pOut.Count() == 0 )
        return;

    // No cuts?
    if( pCutPos.Count() == 0 )
    {
        pSubOuts.AddTail( pOut );
        return;
    }

    // Check
    if( !Outliner::IsClosed(pOut) )
        LzLogException("", "Provided outline is not closed!")

    // Testing if a pos is a cut or not
    std::function<bool(void *)> IsCut = [&]( void * pPos ) -> bool
        {
            BrowseList( iC, pCutPos )
            {
                if( pCutPos.GetAt(iC) == pPos )
                    return true;
            }

            return false;
        };

    // Push current outline
    pSubOuts.AddTail( {} );

    // Process list
    BrowseList( iP, pOut )
    {

LzLogException("", "************ TODO")



    }

}

//================================================================================
void ConvertEdgesToClosedContours( List<EdgeIdx> & pFromEdges,
                                   List< List<EdgeIdx> > & pToContours )
{
    // Setup detonator
    LzServices::Detonator lDet( [&]{ pToContours.DelAll(); },
                                [&]{ pToContours.DelAll(); } );

    // Define Left and Right in an Edge
    auto EdLeft  = []( const EdgeIdx & pEdge ) -> size_t { return pEdge.first; };
    auto EdRight = []( const EdgeIdx & pEdge ) -> size_t { return pEdge.second; };

    // Flip Edge
    auto FlipEd = []( const EdgeIdx & pEdge ) -> EdgeIdx { return {pEdge.second, pEdge.first}; };

    // Define Left and Right in a Strand -- Assuming pStrand has at least one element
    auto StLeft  = [&EdLeft ]( const List<EdgeIdx> & pStrand ) -> size_t { return EdLeft(  pStrand.GetHead() ); };
    auto StRight = [&EdRight]( const List<EdgeIdx> & pStrand ) -> size_t { return EdRight( pStrand.GetTail() ); };

    // Check if contour is closed -- Assuming pStrand has at least one element
    auto IsClosed = [&StLeft, &StRight]( const List<EdgeIdx> & pStrand ) -> bool
        {
            return StLeft(pStrand) == StRight(pStrand);
        };

    // Process all edges
    while( pFromEdges.Count() )
    {
        // Init new strand
        List<EdgeIdx> lStrand;
        lStrand.AddTail( pFromEdges.GetHead() );
        pFromEdges.DelHead();

        // Loop at least once until the contour is closed
        do
        {
            // Find next piece
            void * iE = pFromEdges.HeadPos();
            while( iE )
            {
                // Get edge
                const EdgeIdx & lEd = pFromEdges.GetAt(iE);

                // Does it fit?
                if( EdLeft(lEd) == StRight(lStrand) )
                {
                    // Add piece to the right, and next
                    lStrand.AddTail( lEd );
                    goto FoundOnePiece;
                }
                else
                if( EdLeft(lEd) == StLeft(lStrand) )
                {
                    // Add flipped piece to the left, and next
                    lStrand.AddHead( FlipEd(lEd) );
                    goto FoundOnePiece;
                }
                else
                if( EdRight(lEd) == StLeft(lStrand) )
                {
                    // Add piece to the left, and next
                    lStrand.AddHead( lEd );
                    goto FoundOnePiece;
                }
                if( EdRight(lEd) == StRight(lStrand) )
                {
                    // Add flipped piece to the right, and next
                    lStrand.AddTail( FlipEd(lEd) );
                    goto FoundOnePiece;
                }
                else
                {
                    // Next
                    pFromEdges.Next(iE);
                }
            }

            // Error
            LzLogException("", "Could not find next piece for an open contour!")

        FoundOnePiece:

            // Remove found piece from list
            pFromEdges.DelAtAndNext(iE);
        }
        while( !IsClosed(lStrand) );

        // Stash
        pToContours.AddTail( std::move(lStrand) );
    }

    // All is well
    lDet.defuse();
}
}
}
