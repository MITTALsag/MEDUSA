#pragma once

#include "Triangle.h"
#include "Mesh.h"
#include <LzGeom/Point3D.h>
#include <LzGeom/Vector3D.h>
#include <LzServices/List.h>
#include <vector>
#include <limits>


// DLL exports
namespace LzTriModel
{
    class DLL_EXPORT TopVer;
    class DLL_EXPORT TopTri;
    class DLL_EXPORT TopEdge;
    class DLL_EXPORT MeshTopology;
}


namespace LzTriModel
{
using LzGeom::Point3D;
using LzGeom::Vector3D;
using LzServices::List;
using std::vector;
	

////
//remplacer unsigned int par size_t

/*utiliser*/ static constexpr size_t sUnusedIdx = std::numeric_limits<size_t>::max();



//================================================================================
#pragma region "TopVer"
class TopVer
{
public:
	TopVer();
	~TopVer();

	// If the vertex is unused then both of these lists will remain empty

    List<size_t> mE;
    List<size_t> mT;
};
#pragma endregion


//================================================================================
#pragma region "TopTri"
class TopTri
{
public:
	TopTri( int pV0=-1, int pV1=-1, int pV2=-1 );
	~TopTri();
	int EdgeWithVertex( int pIdxV, int pExcludeEdge=-1 ) const;
	int EdgeWithoutVertex( int pIdxV ) const;
	int VertexOppositeToEdge( int pIdxE, int * pTriIdx=nullptr ) const;
	int OtherEdge( int pExcludeEdge1, int pExcludeEdge2=-1 ) const;
	bool HasSameOrientationAs( const TopTri & pOther ) const;
	bool HasVertex( int pIdxV ) const;

	// Use signed type for safe usage of Other...( -1 ) or Opposite( -1 )
	// mE[0] is opposed to mV[0], etc
	int mV[3];
	int mE[3]; 
};
#pragma endregion


//================================================================================
#pragma region "TopEdge"
class TopEdge
{
public:
	TopEdge( int pV0=-1, int pV1=-1 );
	~TopEdge();
	bool HasVertex( int pIdxV ) const;
	int OtherVer( int pExcludeVer ) const;
	int OtherTri( int pExcludeTri ) const;
    bool IsConnected( const TopEdge & pEdge ) const;

	// Use signed type for safe usage of Other...( -1 )
	int mV[2];
	List<int> mT;
};
#pragma endregion


//================================================================================
#pragma region "MeshTopology"
class MeshTopology
{
#pragma region "Construction & destruction"
public:
    MeshTopology();
    MeshTopology( std::function<void(MeshTopology * pThis)> pInit );
    ~MeshTopology();
		
	// Needs to be called if the associated mesh has been altered (cleaned, etc)
	void Free();

    // Checking compatibility
    bool IsCompatible( const Mesh & pMesh ) const
    {
        return mTopVers.size()==pMesh.mVertices.size()
            && mTopTris.size()==pMesh.mTriangles.size();
    }
#pragma endregion


#pragma region "Analysis"
public:
	// MeshTopology keeps a pointer on pMesh data for commodity. This is not a magnificent design as it should only deal with topology.
	// But helps a lot as MeshTopology serves as a proxy for the mesh in some methods.
	void Set( const Mesh & pMesh );
//bool IsSet() const { return mpVertices != 0; }
bool IsSet() const { return mTopVers.size() != 0; }
    void GetConnectedCracks( List<List<int>> & pTo_ConnectedCracks ) const;
//*********** DEPRECATED
    //void FindPatches( /*bool pMultiEdgesAreBoundaries in this implementation = true!,*/ List< List<size_t> > & pPatches ) const;
    //bool CheckPatchOrientation( const List<size_t> & pPatch ) const;
//*********** DEPRECATED

#pragma endregion


#pragma region "Data"
public:
	// Topology
	const vector<TopVer> & TopVers() const { return mTopVers; }
    const TopVer & TopVers( size_t pIdx ) const;
	//
	const vector<TopTri> & TopTris() const { return mTopTris; }
    const TopTri & TopTris( size_t pIdx ) const;
	//
	const vector<TopEdge> & TopEdges() const { return mTopEdges; }
    const TopEdge & TopEdges( size_t pIdx ) const;
	//
	const vector<int> & Cracks() const { return mCracks; }
	const vector<int> & MultiEdges() const { return mMultiEdges; }

    // Find elements
    bool HasCrackVer( size_t pVer ) const;

protected:
	// Topology
	vector<TopVer> mTopVers;
	vector<TopTri> mTopTris;
	vector<TopEdge> mTopEdges;

	// Secondary
	vector<int> mCracks;
	vector<int> mMultiEdges;
#pragma endregion
};
#pragma endregion
}
