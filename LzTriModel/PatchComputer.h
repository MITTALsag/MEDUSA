#pragma once

#include "MeshTopology.h"


namespace LzServices
{
	template<class T> class List;
	template<class T> class Vector;
}


namespace LzTriModel
{
using LzServices::List;
using std::vector;

class Mesh;

class PatchComputer
{
public:
    PatchComputer() = delete; // To enforce a static usage of the class

public:
    //-----------------------------
    // Mesh splitting
    //-----------------------------

    // Ok if not too many hermetic edges... otherwise hermeticity check not super optimal
    static void SplitMesh( const Mesh & pMesh,
                           List<Mesh> & pPatches,
                           bool pRemoveUnusedversNors,
                           const List<size_t> * ppHermeticEdges=nullptr );

    static void SplitMesh( const Mesh & pMesh,
                           const MeshTopology & pTopo,
                           List<Mesh> & pPatches,
                           bool pRemoveUnusedversNors,
                           const List<size_t> * ppHermeticEdges=nullptr );

    static void * FindPosOfMaxArea( const List<Mesh> & pPatches );

    //-----------------------------
    // Patch counting
    //-----------------------------

    // Ok if not too many hermetic edges... otherwise hermeticity check not super optimal
    static size_t CountPatches( const Mesh & pMesh,
                                const List<size_t> * ppHermeticEdges=nullptr );

    static size_t CountPatches( const MeshTopology & pTopo,
                                const List<size_t> * ppHermeticEdges=nullptr );

//// Better suited for cases with lots of hermetic edges. Execution optimal, but requires allocation of a large hermetic<bool> vector
//static void SplitMesh( const Mesh & pMesh, List<Mesh> & pPatches, const Vector<bool> & ppIsEdgeHermetic );
//static void SplitMesh( const Mesh & pMesh, const MeshTopology & pTopo, List<Mesh> & pPatches, const Vector<bool> & ppIsEdgeHermetic );

    //-----------------------------
    // Multi-edge mesh splitting
    //-----------------------------

    // Hermetic edges = multi-edges
    static void SplitMesh_HermeticME( const Mesh & pMesh,
                                      List< List<size_t> > & pPatchIdx );

    static void SplitMesh_HermeticME( const MeshTopology & pTopo,
                                      List< List<size_t> > & pPatchIdx );

    //-----------------------------
    // Patch extraction from idx
    //-----------------------------

    // Patch reconstruction from lists of triangles
    static void ExtractPatch( const Mesh & pFrom,
                              const List<size_t> & pTris,
                              Mesh & pTo,
                              bool pRemoveUnusedversNors );

//static void SplitMesh( const Mesh & pMesh, List<size_t> & pTriIdx );
//static void SplitMesh( const MeshTopology & pTopo, List<size_t> & pTriIdx );

//****************** TO DO
//****************** TO DO
//void SetOrientationInPatch( const List<size_t> & pPatch, size_t pRefTri );
//void SetMostPopularOrientationInPatch( const List<size_t> & pPatch );
//void InvertPatchOrientation( const List<size_t> & pPatch );
//****************** TO DO
//****************** TO DO

protected:
    static void ExtractFrom( vector<bool> & pDone,
                             const MeshTopology & pTopo,
                             size_t pIdx,
                             List<size_t> * ppPatchIdx,
                             const List<size_t> * ppHermeticEdges );

    static void ExtractFrom_HermeticME( vector<bool> & pDone,
                                        const MeshTopology & pTopo,
                                        size_t pIdx,
                                        List<size_t> & pPatchIdx );
};
}
