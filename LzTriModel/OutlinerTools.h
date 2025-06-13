#pragma once

#include "Outliner.h"


namespace LzTriModel
{
namespace OutlinerTools
{
//********************** RECUP TOUS LES STATICS de Outliner.h
//********************** RECUP TOUS LES STATICS de Outliner.h
//********************** RECUP TOUS LES STATICS de Outliner.h
//********************** RECUP TOUS LES STATICS de Outliner.h
//********************** RECUP TOUS LES STATICS de Outliner.h

/**
 * @brief HasSelfIntersections
 * @param pOut
 * @return
 */
bool HasSelfIntersections( const List<Point3D> & pOut );


/**
 * @brief ProjectOnOutline
 * @param pPt
 * @param pOut
 * @param pDir
 * @param pProj
 * @param ppNormal
 * @return
 */
bool ProjectOnOutline( const Point3D & pPt,
                       const List<Point3D> & pOut,
                       const Vector3D & pDir,
                       Point3D & pProjPt,
                       const Vector3D * ppNormal=nullptr,
                       std::pair<void*, void*> * ppProjSegPos=nullptr);


/**
 * @brief ProjectOnSegments
 * @param pPt
 * @param pSegs
 * @param pDir
 * @param pProjPt
 * @param ppNormal
 * @return
 */
bool ProjectOnSegments( const Point3D & pPt,
                        const List<std::pair<Point3D, Point3D>> & pSegs,
                        const Vector3D & pDir,
                        Point3D & pProjPt,
                        const Vector3D * ppNormal=nullptr );


/**
 * @brief ProjectOnSegment
 * @param pPt
 * @param pA
 * @param pB
 * @param pDir
 * @param pProj
 * @param ppNormal
 * @return
 */
//********************************************************* HIDE IMPLEMENTATION ????
//********************************************************* HIDE IMPLEMENTATION ????
double ProjectOnSegment( const Point3D & pPt,
                         const Point3D & pA,
                         const Point3D & pB,
                         const Vector3D & pDir,
                         Point3D & pProjPt,
                         const Plane3D & pCutPlane );
//********************************************************* HIDE IMPLEMENTATION ????
//********************************************************* HIDE IMPLEMENTATION ????


/**
 * @brief SegmentsToMesh
 * @param pSegs
 * @param pTo
 * @param pScale
 */
void SegmentsToMesh( const List<std::pair<Point3D, Point3D>> & pSegs, Mesh & pTo, double pScale=1.0 );


/**
 * @brief OutlineToMesh
 * @param pOut
 * @param pTo
 */
void OutlineToMesh( const List<Point3D> & pOut, Mesh & pTo, double pScale=1.0 );


/**
 * @brief FindPosOfMaxLength
 * @param pOuts
 * @return
 */
void * FindPosOfMaxLength( const List<List<Point3D>> & pOuts );


/**
 * @brief SplitClosedOutline
 * @param pOut
 * @param pCutPos
 * @param pSubOuts
 */
void SplitClosedOutline( const List<Point3D> & pOut,
                         const List<void *> & pCutPos,
                         List< List<Point3D> > & pSubOuts );
//
void TEST_SplitClosedOutline();


/**
 * @brief ConvertEdgesToClosedContours
 * @param pEdges
 * @param pContours
 */
using EdgeIdx = std::pair<size_t, size_t>;
void ConvertEdgesToClosedContours( List<EdgeIdx> & pFromEdges,
                                   List< List<EdgeIdx> > & pToContours );

}
}
