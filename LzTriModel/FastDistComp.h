#pragma once

#include "FastDistComp_Patching.h"
#include "Mesh.h"
#include "MeshTopology.h"
#include <LzGeom/Vector3D.h>
#include <LzGeom/Tree3D.h>
//#include <LzServices/Vector.h>
#include <LzServices/List.h>


namespace LzTriModel
{
using LzGeom::Vector3D;
using LzGeom::Tree3D;
//using LzServices::Vector;
using LzServices::List;
using std::string;


class FastDistComp
{
#pragma region "Construction / destruction"
public:
    FastDistComp();
    FastDistComp( std::function<void(FastDistComp * pThis)> pInit );
//    FastDistComp( FastDistComp && pOther) = default;
    ~FastDistComp();
	virtual void Free();
#pragma endregion


#pragma region "Setting"
public:
	// When no topology available: will build one locally
	// pMesh can be a temp variable: vertex information is saved in Tree3D
    void Set( const Mesh & pMesh,
              Patching pPatching,
              bool pIgnoreNullNormals=false,
              bool pIgnoreMultiedges=false );

	// if a topology is available: saves memory space
	// pMesh can be a temp variable: vertex information is saved in Tree3D
	// pTopo can be a temp variable
    void Set( const Mesh & pMesh,
              const MeshTopology & pTopo,
              Patching pPatching,
              bool pIgnoreNullNormals=false,
              bool pIgnoreMultiedges=false );

protected:
	class Patch;
    void SetPatch( const Mesh & pMesh,
                   const MeshTopology & pTopo,
                   Patch & pToPatch,
                   bool pIgnoreNullNormals,
                   bool pIgnoreMultiedges );
#pragma endregion


#pragma region "Computing"
public:
    Point3D GetVertex( size_t pPatchIdx, size_t pVerIdx ) const;
    Point3D NearestVertex( Tree3D::NearestMode pTMode, const Point3D & pFrom ) const;

    // Low level information for nearest distance
	class SgnDstInfo
	{
	public:
        size_t mPatchIdx;
        size_t mV[3];
		double mW[3];
        enum class Target { A=0, B=1, C=2, AB=3, AC=4, BC=5, ABC };
		Target mTarget;

        size_t CountTargetVers() const;
        void LogProj( const std::vector<Point3D> & pVers, const Point3D & lCheckPoint ) const;
        void Log( const string & pPrefix="" ) const;
	};


    //********* TODO: add signed dist options
    //********* TODO: add signed dist options
    //********* TODO: add signed dist options
    //********* TODO: add signed dist options
    //********* TODO: add signed dist options
    //********* TODO: add signed dist options


	double SignedDistance( Tree3D::NearestMode pTMode, const Point3D & pFrom, Point3D & pProjP, Vector3D & pProjN, SgnDstInfo * ppInfo=nullptr ) const;
	//
    // Shorthand for Signed Distance
    double AccurateSignedDistance( const Point3D & pFrom ) const
    {
        Point3D lUnused_1;
        Vector3D lUnused_2;
        return SignedDistance(Tree3D::NearestMode::Accurate, pFrom, lUnused_1, lUnused_2);
    }
    //
    // Shorthand for Projection
    Point3D AccurateProjection( const Point3D & pFrom ) const
    {
        Point3D lProj;
        Vector3D lUnused_2;
        SignedDistance(Tree3D::NearestMode::Accurate, pFrom, lProj, lUnused_2);
        return lProj;
    }
    //
	double InOutDistance( Tree3D::NearestMode pTMode, const Point3D & pFrom ) const;

    // Returns true iff the point is inside the mesh (signed dist <= 0)
    bool AccurateIsIn( const Point3D & pFrom ) const;
#pragma endregion


#pragma region "Data"
public:
    // Check number of patches
    size_t GetNbPatches() const { return mPatches.size(); }

	// Privileged access to Tree3D for people who know what they are doing
    const Tree3D & GetTreeForPatch( size_t pPatchIdx ) const;

protected:
	// Triple index (to avoid using Lists)
    struct TriIdx { size_t mI[3]; };

	// Data for a single (connex) patch
	class Patch
	{
	public:
		// Keeps track of the position of vertices in the patch's mesh
		Tree3D mTree;

		// Bounding box of the patch
		BBox mBBox;

		// All the data below uses old indices of the input patches i.e. indices before quick-sorting by Tree3D
        vector< List<size_t> > mVer2Tris;
        vector< TriIdx > mTri2Vers;		// Vertices: 0, 1, 2
        vector< TriIdx > mTri2Edges;	// Edges: 12, 02, 01, see specs in MeshTopology::TopTri::mE

		// Normals computed at initialization time
        vector<Vector3D> mVerNors;
        vector<Vector3D> mEdgNors;
        //vector<Vector3D> mTriNors;
	};

	// All patches
    vector<Patch> mPatches;
#pragma endregion
};
}
