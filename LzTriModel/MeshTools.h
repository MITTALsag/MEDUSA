#pragma once

#include "Mesh.h"
#include "FastDistComp.h"
#include <LzServices/Vector.h>
//
#include <Eigen/Dense>


namespace LzTriModel
{
using LzServices::Vector;



//void ContractMesh( Mesh & pMesh, const List<size_t> & pMobileVers );



#pragma region "Patch extraction"
/**
 * @brief GetTriangleAwayFromCracks
 * @param pMesh
 * @return
 */
size_t GetTriangleAwayFromCracks( const Mesh & pMesh, const MeshTopology & pTopo );

/**
 *
 ********************* Code from LzTriModel::PatchComputer should be moved here!
 *
 **/
//
//* The extracted patch will most likely contain unused vertices. However if there were no duplicate vertices (or normals)
//* in the source mesh then there won't be any duplicate items in the output mesh either.
//
void GetPatchFrom( const Mesh & pMesh,
                   const MeshTopology & pTopo,
                   size_t pFromVerIdx,
                   const List<size_t> & pHermeticEdgIdx,
                   Mesh & pVerMesh );
#pragma endregion


#pragma region "Cutting"
//class CutFromVertexOptions
//{
////public:
////    CutFromVertexOptions() {}

//public:
//    double mCutTol{ 1e-6 };
//    MeshTopology * mpTopo{ nullptr };
//    Mesh * mpOppMesh{ nullptr };
//};
/**
 * @brief CutFromVertex
 * @param pMesh Mesh is not const as vertices will be snapped to cut plane before applying the cut
 * @param pTopo /!\ MUST be the topology associated with the mesh
 * @param pCutPlane
 * @param pVerIdx
 * @param pVerMesh Mesh containing the vertex
 * @param pOppMesh Opposite mesh
 */
void CutFromVertex( const Mesh & pMesh,
                    const MeshTopology & pTopo,
                    size_t pFromVerIdx,
                    const Plane3D & pCutPlane,
                    Mesh & pVerMesh,
                    Mesh & pOppMesh );
#pragma endregion


#pragma region "Mesh analysis & shape analysis"
/**
 * @brief CheckMeshProperties
 * @param pMesh
 * @param pProp
 * @param ppTopo
 * @return
 *
 * @attention if a topology is provided through ppTopo, its compatibility will be checked
 */
bool CheckMeshProperties( const Mesh & pMesh, uint16_t pProp, const MeshTopology * ppTopo=nullptr );
namespace MeshProperties
{
    static constexpr uint16_t NO_CRACKS       = 1;
    static constexpr uint16_t NO_MULTIEDGES   = 2;
    static constexpr uint16_t NO_UNUSED_VERS  = 4;
    //
    static constexpr uint16_t IS_SINGLE_PATCH = 8;
    static constexpr uint16_t IS_MULTI_PATCH  = 16;
    static constexpr uint16_t IS_METERS       = 32;
    static constexpr uint16_t IS_MILLIMETERS  = 64;
    static constexpr uint16_t IS_SMOOTH_NORM  = 128;
};

/**
 * @brief CheckTrianglesConsistency detects degenerate or inverted triangles by comparing the orientation of their normals to
 * the orientation of the normals in the neighbor triangles
 * @param pMesh
 * @param pMaxAngle_deg
 * @param pMinTriArea
 * @param ppTopo
 * @return
 */
bool CheckTrianglesConsistency( const Mesh & pMesh,
                                double pMaxAngle_deg,
                                double pMinTriArea,
                                bool pListAllErrors,
                                const MeshTopology * ppTopo );

/**
 * @brief The PseudoSymPlaneOptions class
 */
class PseudoSymPlaneOptions
{
public:
    // Enabling subsampling
    void EnableSubSamp( size_t pApproxSubSamps );

    // Max nb of steps, and tolerance to reach fixed point
    size_t mMaxSteps = 100;
    double mTol_mm  = 1e-3;
    double mTol_deg = 1e-3;

    // Mesh patching; needed to set-up FDC
    LzTriModel::Patching mPatching = LzTriModel::Patching::Single;

    // Constrained sym plane normal
    bool mConstrainSymNormal = false;
    Plane3D mConstrainedNorPlane;

    // Subsampling
    bool mSubSample = true;
    size_t mApproxSubSamps = 1000;
};

// Find pseudo-symmetry plane of an object: mesh version
Plane3D FindPseudoSymPlane( const Plane3D & pInitPlane, const Mesh & pMesh, const FastDistComp & pFDC,
                            const PseudoSymPlaneOptions & pOptions=PseudoSymPlaneOptions() );

// Find pseudo-symmetry plane of an object: points cloud version
Plane3D FindPseudoSymPlane( const Plane3D & pInitPlane, const vector<Point3D> & pPts, const Tree3D & pTree3D,
                            const PseudoSymPlaneOptions & pOptions=PseudoSymPlaneOptions() );

// Compute Mean, Max, StdDev for the candidate plane: mesh version
vector<double> ComputeSymScore( const Plane3D & pSymPlane, const Mesh & pMesh,
                                const FastDistComp & pFDC, bool pLogValues );

// Compute Mean, Max, StdDev for the candidate plane: points cloud version
vector<double> ComputeSymScore( const Plane3D & pSymPlane, const vector<Point3D> & pPts,
                                const Tree3D & pTree3D, bool pLogValues );

// Returns true if difference > tolerance in mm or deg between the two planes
bool SmallDiff( const Plane3D & pP0, const Plane3D & pP1,
                double pTol_mm=1e-6, double pTol_deg=1e-6,
                const Point3D & pRefPt=Point3D(0,0,0), bool pLog=false );

// Compute volume of a cut object (open mesh), where pOrigin should be on the cut plane
// to account for the fact that the mesh is (potentially) open
double ComputePartialVolume( const Mesh & pCutObj, const Point3D & pOrigin );

// Bigger value of PVAR indicates more volume for less area (e.g. ANTERIOR side of a vertebra)
// Smaller value of PVAR indicates less volume for more area (e.g. POSTERIOR side of a vertebra)
// return value [0] gives the ratio for +pNormal
// return value [1] gives the ratio for -pNormal
vector<double> ComputePartialVolAreaRatio( const Mesh & pObject,
                                           const Point3D & pOrigin, const Vector3D & pNormal );

// Points cloud version
vector<double> ComputePartialVolAreaRatio( const vector<Point3D> & pVolPts, const vector<Point3D> & pSurfPts,
                                           const Point3D & pOrigin, const Vector3D & pNormal );

// Compute scored line intersections
// Score = dot product = (Inter - pO) * pU
// Intersection points are sorted by INCREASING score values
// If pU is normalized then the score is a linear coordinate along the cut axis
// The algorithm might produce duplicate points due to the tolerance used in computing the
// intersection between the line and the triangles
//
// The points are stored in the INCREASING order of Score. Head = smallest score, Tail = biggest score
//
using ScoredPointsList = List< LzServices::Sortable<Point3D, double> >;
void GetLineIntersections( const Mesh & pMesh,
                           const Point3D & pO, const Vector3D & pU,
                           ScoredPointsList & pSortedInters );
#pragma endregion


#pragma region "Geometry generation"
void GenEllipsoid( Mesh & pTo, const Point3D & pCenter,
                   double pXRadius, double pYRadius, double pZRadius,
                   size_t pMeridianSampling=30, size_t pParallelSampling=30 );

// Create a closed cylinder mesh
void GenCylinder( Mesh & pTo, const Point3D & pCenter1, double pRadius1,
                              const Point3D & pCenter2, double pRadius2,
                  size_t pSegments, size_t pLayers );
#pragma endregion


#pragma region "Error measurement"
void AssessDelta_A_to_B( const Mesh & pA, const Mesh & pB, double & pMeanErr, double & pMaxErr, double & pStdDev,
                         const FastDistComp * ppFDC_B=nullptr, const MeshTopology * ppTopo_B=nullptr );
void AssessBilateralDelta( const Mesh & pA, const Mesh & pB, double & pMeanErr, double & pMaxErr, double & pStdDev,
                           const FastDistComp * ppFDC_A=nullptr, const FastDistComp * ppFDC_B=nullptr,
                           const MeshTopology * ppTopo_A=nullptr, const MeshTopology * ppTopo_B=nullptr );
//#pragma region "Error measurement"
//vector<double> GetMaxErrors( const Mesh & pA, const FastDistComp * ppFDC_A, const Mesh & pB, const FastDistComp * ppFDC_B );
//#pragma endregion
#pragma endregion





#pragma region "ICP fit"
//*** RigidFit pourrait utiliser le principe du FatFit pour eviter de tomber dans les minimas locaux
//*** RigidFit pourrait utiliser le principe du FatFit pour eviter de tomber dans les minimas locaux
//*** RigidFit pourrait utiliser le principe du FatFit pour eviter de tomber dans les minimas locaux
//*** RigidFit pourrait utiliser le principe du FatFit pour eviter de tomber dans les minimas locaux

//TxLogWarningMessage("Checker code ICP....","TxMath::ToolBox");
//check mais apparemment ca marche bien !
//********* retenir la meilleure distance au cours de toutes les iterations (sont elles calculees ?? verifie-t-on que l'erreur est strictement decroissante ????)
//********* retenir la meilleure distance au cours de toutes les iterations (sont elles calculees ?? verifie-t-on que l'erreur est strictement decroissante ????)
//********* retenir la meilleure distance au cours de toutes les iterations (sont elles calculees ?? verifie-t-on que l'erreur est strictement decroissante ????)
//*** Dans Arun loop, projeter sur la surface et ne pas prendre le point le plus proche
//*** ce qui nous oblige a passer un FDC et pas un tree3D !!
RigidTr3D Rigid_ICP( const Mesh & pSrcMesh, const Mesh & pDstMesh, double pMaxEulerRotDeg=0.1, double pMaxTransNorm=0.1, size_t pMaxSteps=300 );
Matrix Rigid_ICP( const Mesh & pSrcMesh, const Mesh & pDstMesh, const Matrix & pAffinePose, double pMaxEulerRotDeg=0.1, double pMaxTransNorm=0.1, size_t pMaxSteps=300 );
#pragma endregion


#pragma region "PCA fit"
Matrix RegisterMinMax( const Mesh & pSrc, const Mesh & pDst, const RigidTr3D & pSrc2Wrld, const RigidTr3D & pDst2Wrld );

// Affine registration: scaling along principal axis ENHANCED (see. TxMath::ToolBox)
// Returns an affine matrix that transforms Src to Dst while minimizing SSD point (src) to surface (dst) distance, within pMaxIter iterations.
// Minimizes an asymetrical energy as only distance from Src to Dst is taken into account in the SSD.
//
//*** Creates a FastDistComp for all destination meshes or uses the one(s) provided by client.
//
// Optional parameter ppLockedScales lists the directions (0, 1, or 2) which will have scales fixed to 1.0 during the iteration
Matrix AffinePCAFitPairedMeshes( Mesh * pSrc, const Mesh * pDst, size_t pMaxIter, double pTolToIdentity, List<size_t> * ppLockedScales=nullptr,
                                 const FastDistComp * ppDstFDC=nullptr, const MeshTopology * ppDstTopo=nullptr );
Matrix AffinePCAFitPairedMeshes( const Vector<Mesh*> & pSrc,
                                 const Vector<const Mesh*> & pDst, size_t pMaxIter, double pTolToIdentity, List<size_t> * ppLockedScales=nullptr,
								 const Vector<const FastDistComp*> * pDstFDCs=nullptr, const Vector<const MeshTopology*> * pDstTopos=nullptr );
//*** TODO
//
//*** To be used from within a loop when the destination data does not move. Uses a FastDistComp created outside.
//
//Matrix AffinePCAFitPairedMeshes( Mesh * pSrc, const FastDistComp * pDst, size_t pMaxIter, double pTolToIdentity );
//Matrix AffinePCAFitPairedMeshes( const Vector<Mesh*> & pSrc,
//								 const Vector<const FastDistComp*> & pDst, size_t pMaxIter, double pTolToIdentity );
//*** TODO

//// Surfacic PCA
//// returns Alias PCA_to_W
//RigidTr3D SurfacicPCA( const Mesh & pObject );

// Compute this meshe's inertial PCA
void TEST_InertiaTensor( const Mesh & pObject, Point3D & pMean, Eigen::Matrix3d & pI );
void InertialPCA( const Mesh & pObject, Point3D & pMean, Vector3D pV[3] );
void TEST_InertialPCA();
#pragma endregion


#pragma region "Elastic fit"
//************************************ MAGIC NUMBERS
//pFitIter = 20;
//pSmoothIter = 5;
//************************************ MAGIC NUMBERS
void ElasticFit( Mesh & pFrom, const Mesh & pTo, size_t pFitIter, size_t pSmoothIter, bool pBilateral,
                 List<size_t> * ppSrcIdx=nullptr, List<Point3D> * ppDstPos=nullptr );
#pragma endregion


#pragma region "Extending from cracks"
// Extends a mesh from a closed crack in a direction given by the extension plane.
// Uses the provided topology or computes its own.
//
// /!\ pMesh can be left in an incoherent state if an exception is raised during the call.
// /!\ It is up to the user to do the proper backup/restore and exception check operations.
//
// - Identifies a set of cracks above the extension plane to extend from
// - The cracks must form a unique, ideally convex, closed contour
// - All crack vertices are projected on a least square plane having the same normal as the extension plane
// - Extends the surface in the direction of the normal of the extension plane
// - Generates pLayersCount layers, each having a pLayerHeight
// - Each crack vertex has its own extension direction computed by sampling the surface's shape underneath
// - The sampling is done below the least square extension plane by pSamplesCount samples spaced every pSampleHeight
// - The final mesh is Laplacian relaxed in pLaplaceRelaxSteps iterations
// - All vertices below (or on) pExtPlane are fixed during relaxation as well as the vertices on the new top crack 
//
void ExtendFromCracks( Mesh & pMesh, const MeshTopology * ppMeshTopo,
					   const Plane3D & pExtPlane,
                       size_t pLayersCount, double pLayerHeight,
                       size_t pSamplesCount, double pSampleHeight,
                       size_t pLaplaceRelaxSteps );

// Mesh repulsion (to avoid conflict between deformed capsule and patient's skin, for example; see AtlasTransfert)
void RepulseMesh( Mesh & pMesh, const Mesh & pRepulsor, double pRepDist, double pRelaxRadius );

// Bilateral mesh repulsion
// At the end, both meshes must be at least pMinCollDist away from one another
class SolveBilatConflictOptions
{
public:
    // Debug save intermediate files
    bool mDebugSave = false;

    // Patching
    Patching mPatching_A = Patching::Single;
    Patching mPatching_B = Patching::Single;

    // Maximal number of iterations for conflict solving
    size_t mMaxIter = 20;

    // Conflict solving rate
    double mSolveRate = 0.5;  // Percentage of conflict solving applied at each displacement

    // Distance above which (0 meaning 'no conflict') the conflict is solve as
    // a one shot displacement of the involved nodes, possibly meaning a bilateral
    // cross over which would result in an unnecessary collision margin of mOneShotMinDist
    double mOneShotMinDst = -0.01;

    // HC Laplacian rexation (smoothing) steps at each conflict solving stage
    size_t mHCSmoothSteps = 5;
    //
    // Warning: too many HC iterations while fixing the boundaries are likely to produce concavities
    //          from otherwise nicely convex meshes, so one single relaxation step to even things out
    //          should be enough. We are also expecting very small collisions.

    // Mesh relaxation radius around the (most recently i.e., at current iteration) displaced nodes
    double mRelaxRadius = 5.0; // mm
};
/**
 * @brief SolveBilatConflict
 * @param pA
 * @param pB
 * @param pOptions
 * @return Iteration index at which the conflict is solved. If 0 then there was no conflict in the first place.
 */
size_t SolveBilatConflict( Mesh & pA, Mesh & pB,
                           const SolveBilatConflictOptions & pOptions=SolveBilatConflictOptions() );
//
// Unit test SolveBilatConflict
void UNIT_TEST__SolveBilatConflict();
#pragma endregion


#pragma region "Mesh refinement"
enum class RefineMode
{
	AbsoluteLength, // Splits edges in half if len(edge) > pThreshold
	//--> Iterations converge as edge lenghts are getting smaller

	RelativeLength  // Splits edges in half if len(edge) > pThreshold * len(shortest edge)
	//--> Iterations may not converge
};
// Lists edges that must be refined. Uses MeshTopology numbering (edge 0 is opposed to V0, ...)
//void RefineTriangle( const Point3D & pV0, const Point3D & pV1, const Point3D & pV2, RefineMode pMode, List<size_t> & pRefEdges );
void RefineTriangle( RefineMode pMode, double pThreshold, const Point3D pVers[3], List<size_t> & pRefEdges );
// Refines whole mesh. Returns number of bad triangles i.e. still violating RefineMode rule.
size_t RefineMesh( RefineMode pMode, double pThreshold, Mesh & pMesh );
#pragma endregion


#pragma region "Mesh simplification"
// Topology preserving cleaner
//////////////////////CrackPolicy : crackscannotbecollapsed, 
//////////////////////CrackPolicy : crackscannotbecollapsed, 
//////////////////////CrackPolicy : crackscannotbecollapsed, 
void TopoClean( Mesh & pMesh, double pVerTol, double pEdgTol, bool pPreserveCracks, bool pCheckMEs );
#pragma endregion


#pragma region "Mesh relaxation"
// Moves vertices to the mean position of all vertices surounding current vertex
// Current vertex is not considered in the computation of new mean position (strictly Laplacian approach)
//
// Applies sliding constraint: signed_dist(v) >= 0 for all vertices
// *** NON: fausse bonne idee.... pVerNeighs: number of neighbor generations involved in computation of new position
// *** desequilibre la relaxation aux bords, et ne lisse pas la surface
//
void SlideVerticesToMeanPos( Mesh & pMesh, const FastDistComp & pOnto, size_t pSteps,
                             /*size_t pVerNeighs=1,*/ const MeshTopology * ppTopo=nullptr, const Vector<bool> * pIsVerFixed=nullptr );


// Laplacian smoothing followed by a smooth retro-projection
//****************************************
//****************************************
//**
//**
//**
//** CA MARCHE !!!!!!
//**
//**
//**
//****************************************
//****************************************
void SmoothSmoothing( Mesh & pMesh, size_t pSteps, /*, size_t pFitIter, size_t pSmoothIter,
                      size_t pTriNeighs=1, */const MeshTopology * ppTopo=nullptr, const Vector<bool> * pIsVerFixed=nullptr );


#pragma region "....failures"
// Mesh smoothing
//
// **** testing idee : pTriNeighs..........
//
void TangentMeshSmoothing( Mesh & pMesh, size_t pSteps, size_t pTriNeighs=1, const MeshTopology * ppTopo=nullptr, const Vector<bool> * pIsVerFixed=nullptr );

// Mesh smoothing
//
// **** testing idee : pTriNeighs..........
//
void TanVerSmooth( Mesh & pMesh/*, const FastDistComp & pOnto*/, size_t pSteps, size_t pVerNeighs=1, const MeshTopology * ppTopo=nullptr, const Vector<bool> * pIsVerFixed=nullptr );
void TanVerSmooth_2( Mesh & pMesh, const FastDistComp & pOnto, size_t pSteps, size_t pVerNeighs=1, const MeshTopology * ppTopo=nullptr, const Vector<bool> * pIsVerFixed=nullptr );

// Mesh normals' smoothing
void SmoothPerVerNormals( Mesh & pMesh, size_t pSteps, size_t pVerNeighs=1, const MeshTopology * ppTopo=nullptr );

#pragma endregion
#pragma endregion


}
