#pragma once

#include "Coord3D.h"
#include "Point3D.h"
#include "BBox.h"
#include <vector>


namespace LzGeom
{
using std::vector;


class Tree3D
{
#pragma region "Construction / destruction"
public:
    Tree3D();
    Tree3D( const Tree3D & pOther ) { *this = pOther; }// = delete;
    Tree3D( Tree3D && pOther) = default;
    virtual ~Tree3D();
    virtual void Free();
#pragma endregion


#pragma region "Copying"
public:
//    const Tree3D & operator=( const Tree3D & pOther ) = delete;
    const Tree3D & operator=( const Tree3D & /*pOther*/ )
    {
        // Free current
        Free();

        //        implementer le bidule ! ==> shared ptr
        //        implementer le bidule ! ==> shared ptr
        //        implementer le bidule ! ==> shared ptr
        //        implementer le bidule ! ==> shared ptr
        //        implementer le bidule ! ==> shared ptr
        //        implementer le bidule ! ==> shared ptr
LzLogException("", "NOT IMPLEMENTED")
        return *this;
    }
#pragma endregion


#pragma region "Set-up"
public:
    //
    // Calls Node::UpdateTree to recursively subdivide the points cloud into "bunches" where separators are
    // stored at the tree's branching nodes. Terminal nodes (leaves) contain the actual data. No data is stored
    // at branching (non-terminal) nodes.
    //
    // pLogStats:
    //      Log tree stats or not.
    //
    // pMinDiameter:
    //      A bunch of points whose bbox is smaller than pMinDiameter in all 3 dimensions is considered ponctual.
    //      The corresponding node is not subdivided.
    //
    // pMaxLeafCount:
    //      Should be MinSubdiv. Number of nodes below which a node is not subdivided.
    //      If number of points in bunch <= pMaxLeafCount in a node then the node is not subdivided.
    //
    // pMaxTreeDepth:
    //      Maximal depth of the tree. 10 subdivisions will produce 1024 leaves.
    //
    void Create( const vector<Point3D> & pPoints, bool pLogStats=true, double pMinDiameter=0.001, unsigned int pMaxLeafCount=10, unsigned int pMaxTreeDepth=10 );
    void Create( const Vector<Point3D> & pPoints, bool pLogStats=true, double pMinDiameter=0.001, unsigned int pMaxLeafCount=10, unsigned int pMaxTreeDepth=10 );
#pragma endregion


#pragma region "Computation"
public:
    enum class NearestMode { Fast, Accurate };
	// Returns index in old indexation (i.e. the one provided by caller, before quicksort in Tree3D)
    unsigned int GetNearestVertex( const Point3D & pPt, NearestMode pMode ) const;
//double GetDistance( const Point3D & pPt, const Mesh * ppTriModel, Point3D & pBestPoint );

	// Returns a list of old indices of all vertices within range < pRange to point pPt
	// This gets more costly than linear search when range increases (as we get tree filtering of leaves)
	// For small ranges ==> very efficient especially when pRangeIdx is empty... (no intersection)
    void GetRangeVertices( const Point3D & pPt, double pRange, List<size_t> & pRangeIdx ) const;

    // Access to vertices using old indexation (i.e. indexation before quicksort in Tree3D)
    const Point3D & OldPoint3D( size_t pOldIdx ) const;
#pragma endregion


#pragma region "Debug"
#if defined(USING_QT) && !defined(NO_OPENGL)
public:
    void DrawLeaves() const;
    void DrawPlanes() const;
    void Draw2D( double pWidth, double pHeight ) const;
    //void DrawCandidatesTri(const unsigned int pNearest);
//void Draw(const bool & pIsMeters) const { mpTree->Draw(); };
//void Log() const {mpTree->LogDisplay(0, mpPoints);};
//protected:
//  void Draw2D_Recurse( const Node * ppNode, const Point3D & pPos ) const;
#endif
#pragma endregion


#pragma region "Data"
protected:
    friend class Node;

	// Statistics
    unsigned int mDepthMax;
    unsigned int mNbLeaf;
    double mMoyLeaf;
    unsigned int mSizeMax;
    unsigned int mSizeMin;
//    unsigned int mConcurrentLeaves; //******************* debug data only

    // Root
    class Node;
    Node * mpRootNode;

    // QSorted points cloud
    Vector<Point3D> mPoints;
    Vector<unsigned int> mNew2Old;
    Vector<unsigned int> mOld2New;
#pragma endregion


#pragma region "Node"
protected:
    class Node
    {
#pragma region "Construction / Destruction"
    public:
        Node();
        virtual ~Node();
#pragma endregion


#pragma region "Computations"
    public:
        void UpdateTree( Vector<Point3D> & pPoints, Vector<unsigned int> & pIndices,
//Tree3D & pTree, d'abord mettre au point test auto puis refondre ce truc
						 unsigned int & pDepthMax, unsigned int & pNbLeaf, unsigned int & pSizeMax, unsigned int & pSizeMin,
						 int pBegin, int pEnd, unsigned int pCurrDepth,
						 unsigned int pMaxLeafCount, double pMinDiameter, unsigned int pMaxDepth );

        const Node * GetNearestLeaf( const Point3D & pPt ) const;
        unsigned int GetNearestPoint( const Point3D & pPt, const Vector<Point3D> & pPoints ) const;

		// Accurate leaf exploration for nearest vertex
        void GetNextLeaf( unsigned int & pConcurrentLeaves, unsigned int & pNearest, double & pMinDist, const Point3D & pPt, const Vector<Point3D> & pPoints ) const;

		// Leaf exploration for range vertices
        void CollectRangeVertices( const Point3D & pPt, double pRange, List<size_t> & pRangeIdx, const Vector<Point3D> & pPoints, const Vector<unsigned int> & pNew2Old ) const;
#pragma endregion


#pragma region "Display"
        void DrawLeaves( double pMoyLeaf, unsigned int pSizeMax, unsigned int pSizeMin/*, const bool pIsMeters*/ ) const;
        void DrawPlanes() const;
//void LogDisplay(const unsigned int pDepth, const vector<Point3D> * ppPoints) const;
//void Draw() const;
        void DrawTree2D( const Point3D & pPos, double pCurrW, double pOffY ) const;
#pragma endregion


    public:
        // Data
        BBox mBBox;

    protected:
        // Data
        Node * mpNodeLeft;
        Node * mpNodeRight;
        double mSeparator;
        int mIdxStart;
        int mIdxEnd;
        int mAxis;
        unsigned int mSize;

        // Back-office
        int BestAxis( double pMinDiameter ) const;
        double Pivot( const Vector<Point3D> & ppPoints, int & pMaxLeftIdx ) const;
        bool IsChoice( const Point3D & pPt, double pMinDist ) const;
//void DrawHyperplane() const;
    };
#pragma endregion
};
}
