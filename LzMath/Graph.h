#pragma once

#include "LzLib_Math.h"
#include <LzServices/List.h>


// DLL exports
namespace LzGraph
{
    class DLL_EXPORT_LzMath Node;
    enum class DLL_EXPORT_LzMath PathMode;
    class DLL_EXPORT_LzMath Graph;
}


namespace LzTriModel { class Mesh; class MeshTopology; }


namespace LzGraph
{
using LzServices::List;
using LzTriModel::Mesh;
using LzTriModel::MeshTopology;
using std::vector;


/**
 *
 * TODO: geodesic approximation. From a tri mesh, build dual graph based on triangle centers
 * and per-edge connectivity. Find min len path. Adjust path by sliding the vertices on the
 * edges between successive triangles. Init = points in the middle of the edge. Adjust param
 * position of point on edge as a percentage [0, 1].
 *
 *
 *
 *
 */

//==============================================================================
class Node
{
public:
    // Assumption: both list hold the same number of elements
    List<size_t> mNextVers;
    List<size_t> mNextEds;
    List<double> mLens;

    // Path finding
    void * mPosInQ{ nullptr };
    size_t mPrevVer; // Makes sense only once mPosInQ has become nullptr (i.e., the node has been visited)
    size_t mPrevEd;  // Makes sense only once mPosInQ has become nullptr (i.e., the node has been visited) AND pMode == EdgeIndices mode
    double mDist{ std::numeric_limits<double>::max() };
};


//==============================================================================
enum class PathMode
{
    VertexIndices,
    EdgeIndices
};


//==============================================================================
class Graph
{
public:
    void Init( const Mesh & pMesh, const MeshTopology * ppTopo=nullptr );

    // If pFrom == pTo, then the returned path list is empty
    void FindShortestPath( size_t pFrom, size_t pTo, List<size_t> & pIdx, PathMode pMode );

protected:
    vector<Node> mNodes;
};

//==============================================================================
void DLL_EXPORT_LzMath GetPathFromWaypoints( const vector<size_t> & pWaypoints,
                                             List<size_t> & pPath,
                                             const Mesh & pMesh,
                                             const MeshTopology & pTopo );
}
