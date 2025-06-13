#pragma once

#include "Mesh.h"
#include "MeshTopology.h"
#include <LzServices/Vector.h>


namespace LzTriModel
{
using LzGeom::Vector3D;
using LzTriModel::Mesh;
using LzTriModel::MeshTopology;
using LzServices::Vector;
//
using std::vector;


#pragma region "EnhancedVertex"
//================================================================================
class EnhancedVertex
{
public:
    EnhancedVertex();
    virtual ~EnhancedVertex() {}

    int mIdx;
    List<int> mTriangles;
    List<int> mNeighbors;
    int mCollapseTgt;
    double mCollapseCost;

    bool IsValid() const { return mIdx != -1 && mCollapseTgt != -1; }
    void Log() const;
};
#pragma endregion


#pragma region "EnhancedTriangle"
//================================================================================
class EnhancedTriangle
{
public:
    EnhancedTriangle() : mIdx(-1), mHasFlipped(false) {}
    virtual ~EnhancedTriangle() {}

    int mIdx;
    Vector3D mNormal; 
    Vector3D mLastNormal;
    bool mHasFlipped;

    bool IsValid() const { return mIdx != -1; }
    bool HasFlipped() const { return mHasFlipped; }
};
#pragma endregion


#pragma region "MeshDecimator"
//================================================================================
class MeshDecimator
{
#pragma region "Construction & destruction"
public:
// TODO : Mettre à jour => Methodes statiques uniquement ?
    MeshDecimator() : mCostMin(-1), mCostMax(-1) {};
	~MeshDecimator();
	void Free();
#pragma endregion
        

#pragma region "Data"
public:
    // Mesh's vertices & triangles enhanced versions
    vector<EnhancedVertex> mEV;
    vector<EnhancedTriangle> mET;

    // Min, Max cost of non preserved edge (thus mCostMax != DBL_MAX)
    double mCostMin;
    double mCostMax;
#pragma endregion


#pragma region "Decimate - GameDeveloper algo 1998"
// Decimate a mesh using edge collapse, see Doc/GameDeveloper - 1998.pdf
// Each vertex is enhanced with information: its best neighbor collapsing candidate and its cost
// See EnhancedVertex class
// Cost is function of edge length and curvature difference between both triangles sharing this edge
public:
    void Decimate_GameDeveloper1998( Mesh & pMesh, const MeshTopology & pTopo, unsigned int pTgtTrianglesNo, bool pPreserveCracks=false, double pMaxNormalAngles_Deg=180 );

protected:
    double ComputeCost( const Mesh & pMesh, const MeshTopology & pTopo, const EnhancedVertex & pVOrg, const EnhancedVertex & pVTgt, bool pPreserseCracks, double pMaxNormalAngles_Deg );
    void ComputeBestCost( const Mesh & pMesh, const MeshTopology & pTopo, EnhancedVertex & pV, bool pPreserseCracks, double pMaxNormalAngles_Deg );
#pragma endregion


#pragma region "Edge Collapse Cost"
public:
    void ComputeEnhancedVertexStructure( const Mesh & pMesh, const MeshTopology & pTopo );
#ifdef USE_CLI
    void DrawCollapseCostEdge( const Mesh & pMesh );
#endif
#pragma endregion
};
#pragma endregion


}
