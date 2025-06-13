#pragma once

#include "LzLib_TriModel.h"
#include "Triangle.h"
#include <LzGeom/Point3D.h>
#include <LzGeom/Vector3D.h>
#include <LzGeom/Plane3D.h>
#include <LzGeom/BBox.h>
#include <LzServices/List.h>
#include <vector>


// DLL exports
namespace LzTriModel
{
    class DLL_EXPORT Mesh;
}


namespace LzMath { class Matrix; }


namespace LzTriModel
{
using LzGeom::Point3D;
using LzGeom::Vector3D;
using LzGeom::Line3D;
using LzGeom::Plane3D;
using LzGeom::RigidTr3D;
using LzGeom::BBox;
using LzMath::Matrix;
using LzServices::List;
using LzServices::Vector;
using std::vector;


class MeshTopology;


class Mesh
{
#pragma region "Construction & destruction"
public:
	Mesh();
    Mesh( const Mesh & pMesh );
    Mesh( Mesh && pMesh );
    Mesh( std::function<void(Mesh * pThis)> pInit );
	virtual ~Mesh();

	virtual void Free();
#pragma endregion


#pragma region "Load & save"
public:
	virtual void Load( const std::string & pFName );
	virtual void Save( const std::string & pFName ) const;
protected:
	void LoadXml( const std::string & pFName );
	void SaveXml( const std::string & pFName ) const;
	void LoadObj( const std::string & pFName );
	void LoadEncryptedObj( const std::string & pFName );
	void SaveObj( const std::string & pFName ) const;
	void SaveEncryptedObj( const std::string & pFName ) const;
	void LoadMdl( const std::string & pFName );
	void SaveMdl( const std::string & pFName ) const;
	void LoadStl( const std::string & pFName );
	void SaveStl( const std::string & pFName ) const;
	//void LoadMesh( const std::string & pFName );
	void SaveMesh( const std::string & pFName ) const;
#pragma endregion


#pragma region "Operations"
public:
    // Mesh loaded
    //
    // A mesh is considered as loaded as soon as there are 1 or more vertices.
    // A loaded mesh can thus be a mesh without triangles (points cloud) or a
    // surface mesh (defined by triangles) without normals.
    bool IsLoaded() const { return mVertices.size() != 0; }

	// Costly operation: computes bbox of all vertices
    bool IsMeters( double pMaxLengthInMeters=3.0
    /* no human body larger than 3 meters...
       and most body parts are larger than 3 mm*/ ) const;

    // Returns number of unused vertices and normals
    size_t RemoveUnusedVerticesAndNormals();
//********************************* bool RemoveUnusedTexCoords();
	void RemoveAllTexCoords();

//*************************** RemoveDuplicateNormals( double pMergeTol ) : plus clair, au lieu de dupliquer le code de RemoveDupVers
    // Returns true (i.e., "Finished") if no more triangles/vertices need to be removed
	bool RemoveDuplicateVertices( double pMergeTol=1e-6 );
	bool RemoveDuplicateVerticesAndNormals( double pMergeTol=1e-6 );
//*************************** RemoveDuplicateNormals( double pMergeTol ) : plus clair, au lieu de dupliquer le code de RemoveDupVers

    /**
     * @brief RemoveCollapsedTriangles removes all triangles with duplicate vertices
     * @return Returns the number of erased triangles
     */
    size_t RemoveCollapsedTriangles();

    /**
     * @brief RemoveOppositeTriangles
     * @return Returns the number of erased triangles
     */
    size_t RemoveOppositeTriangles();



    Mesh & operator=( const Mesh & pMesh );
    Mesh & operator=( Mesh && pMesh );
    //
    void Flip( const Plane3D & pMirror );
	void Scale( double pScale ) { Scale( pScale, pScale, pScale ); }
	void Scale( double pScaleX, double pScaleY, double pScaleZ );
	void RigidTransform( const RigidTr3D & pT );
	void AffineTransform( const Matrix & pT );
	void ElasticTransform( const Vector<Vector3D> & pDispPerVer );
    //
    void Merge( const Mesh & pOther );
//    void Merge( Mesh && pOther ); .... to do ?
    void MergeTrianglesOnly( const Mesh & pOther );
    //
    void KeepOnlyTriangles( const vector<size_t> & pTriIdx );
	void SwapTrianglesOrientation();

	// Removes any triangle with vertices on the negative side of the cut plane.
	// Does not remove unused vertices. No new vertices added.
	void UglyCut( const Plane3D & pCut );

	// Removes any triangle with vertices outside of the bbox.
	// Does not remove unused vertices. No new vertices added.
	void UglyCut( const BBox & pCut );

	// Removes the part of the mesh located on the negative side of the cut plane.
	// Does not remove unused vertices. Adds many new vertices on the cut plane.
    void Cut( const Plane3D & pCut, double pSnapRange,
              const vector<Point3D> * pScalpelQuad=nullptr, const Plane3D * pCutBoundary=nullptr );

	// Removes the part of the mesh located inside or outside another mesh.
	// Does not remove unused vertices. Adds many new vertices.
	enum class CutMeshMode { RemoveInside, RemoveOutside };
    void Cut( const Mesh & pCut, double pSnapRange, CutMeshMode pMode/*, LzTriModel::Patching pPatching=LzTriModel::Patching::Split*/ );

#pragma region "Laplacian relax"
	// Moves vertices to the mean position of all vertices surounding current vertex
    void MoveVerticesToMeanPos( size_t pSteps, const MeshTopology * ppTopo=nullptr, const vector<bool> * pIsVerFixed=nullptr );

	// Moves vertices to the mean position of all vertices surounding current vertex.
	// Recovers initial cartesian bbox
    void MoveVerticesToMeanPosAndRecoverBBox( size_t pSteps, const MeshTopology * ppTopo=nullptr, const vector<bool> * pIsVerFixed=nullptr );

	// Moves vertices to the mean position of all vertices surounding current vertex using the HC algorithm (which tends to keep the original shape by moving back the vertices towards their original position)
    void MoveVerticesToMeanPos_HC( size_t pSteps,
                                   const MeshTopology * ppTopo=nullptr,
                                   const vector<bool> * ppIsVerFixed=nullptr,
                                   const vector<const Plane3D *> * ppSlidePlanes=nullptr,
                                   double pAlpha=0.4, double pBeta=0.6 );
#pragma endregion

	// Computes smooth normals
	void SmoothNormalize( const MeshTopology * ppTopo=nullptr );

	// Computes sharp normals
	void SharpNormalize( double pMinSharpAngleDeg, bool pUseAreaPonderation, const MeshTopology * ppTopo=nullptr );

	// Extrusion / intrusion
	void NormalExtrude( double pNormalShift, bool pUseAreaPonderation );
	void NormalExtrude( double pNormalShift_X, double pNormalShift_Y, double pNormalShift_Z, bool pUseAreaPonderation );
    void NormalExtrude( const vector<double> & pHeightMap );
    void NormalExtrude( const vector<double> & pHeightMap, bool pUseAreaPonderation );

	// Converts each triangle to 4 smaller triangles
    void Resample_4xTri( /*const List<size_t> * ppPatch = nullptr */ );

	// Creates tri mesh from a bounding box
	void FromBBox( const BBox & pBBox );


//// Finds loops of nearest vertices between this mesh and another. Explores all attractors starting from all vertices in this mesh.
//void FindAttractors( const Mesh & pMeshB, List< List<size_t> > & pAtoB ) const;
		




//*********** BEFORE USING GetLineIntersections CHECK THE CORRESPONDING FUNCTIONS IN MeshTools
//*********** BEFORE USING GetLineIntersections CHECK THE CORRESPONDING FUNCTIONS IN MeshTools
    // Lists intersections along axis
    // Uses matrix inversion which can have bad conditioning based on a triangle's aspect ratio!
    //
    // **** IDEA: use FDC + tolerance distance to accept / reject intersection points
    // ****       This would be more costly but probably a possible optimisation would be
    // ****       to use the intersecting triangle as a defining search neighborhood
    // ****
    //
    using ScoredPointsList = List< LzServices::Sortable<Point3D, double> >;
    void GetLineIntersections_OLD( const Line3D & pLine,
                                   List<Point3D> & pInters,
                                   ScoredPointsList * pScoredInters=nullptr ) const;

	// Optimization for repeated intersection computation
	void ComputeTriMats( Vector<Matrix> & pTriMats ) const;
    void GetLineIntersections( const Line3D & pLine, List<Point3D> & pInters,
                               const List<size_t> & pTriIdx, const Vector<Matrix> & pTriMats ) const;
//*********** BEFORE USING GetLineIntersections CHECK THE CORRESPONDING FUNCTIONS IN MeshTools
//*********** BEFORE USING GetLineIntersections CHECK THE CORRESPONDING FUNCTIONS IN MeshTools





	// Centroid of all faces ponderated by each face's area
    void GetFacesCentroids( vector<Point3D> & pCentroids, vector<double> & pAreas ) const;
    void GetFacesCentroids( Vector<Point3D> & pCentroids, Vector<double> & pAreas ) const;
    Point3D WeightedCentroid() const;
		
	// Volumetric centroid: centroid of all tetrahedrons, weighted by their respective volume
    void GetVolumetricCentroids( vector<Point3D> & pCentroids, vector<double> & pVolumes ) const;
    void GetVolumetricCentroids( Vector<Point3D> & pCentroids, Vector<double> & pVolumes ) const;
    Point3D VolumetricWeightedCentroid() const;

	// Compute most mirroring plane
    size_t GetBestMirrorPlane( const Vector<Plane3D> & pPlanes ) const;

	// Index of nearest vertex
    size_t NearestVertexIndex( const Point3D & pPt ) const;
	const Point3D & NearestVertex( const Point3D & pPt ) const	{ return mVertices[ NearestVertexIndex(pPt) ]; }

    // Extremum in a given direction
    size_t GetExtremumIndex( const Vector3D & pAlongDir ) const;
    const Point3D & GetExtremum( const Vector3D & pAlongDir ) const { return mVertices[ GetExtremumIndex(pAlongDir) ]; }
#pragma endregion


#pragma region "Properties"
public:
    double Volume() const;
	double Area() const;
	// ProbeInversion implemente un algorithme un peu naze pour tester une inversion de maillage....
	//bool ProbeInversion( DistanceComputer::Patching pPatch ) const;
#pragma endregion


// #pragma region "Draw"
public:
 	void ResetDisplayList() {}
// #if !defined(NO_OPENGL)
//     virtual void Draw();
// //    virtual void Draw( std::function<void(size_t)> pSetVerColor );
// #endif
// #pragma endregion


#pragma region "Data"
public:
    class /*DLL_EXPORT*/ TexCoord2D
	{
	public:
		TexCoord2D( double pS=0, double pT=0 ) : mS(pS), mT(pT) {}

		double mS;
		double mT;
	};

    vector<Point3D> mVertices;
    Vector<TexCoord2D> mTexCoords2D;
    vector<Vector3D> mNormals;
    vector<Triangle> mTriangles;

	// Display list
    List<size_t> mDLs;
#pragma endregion
};
}


