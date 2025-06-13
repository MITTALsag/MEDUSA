#pragma once

#include "Triangle.h"
#include "MeshTopology.h"
#include <LzGeom/Vector3D.h>
#include <LzGeom/Plane3D.h>


namespace LzTriModel
{
using LzGeom::Point3D;
using LzGeom::Vector3D;
using LzGeom::Plane3D;
using LzServices::List;
//
using std::vector;


class Mesh;
//class MeshTopology;


//================================================================================
class OutVer
{
public:
	double mD;
	int mStatus;
	// Status semantics:
	//		-1 if D<-tolerance: vertex is below cut plane
	//		 0 if D in [-tolerance,+tolerance]: vertex is on cut plane
	//		+1 if D>+tolerance: vertex is above cut plane
};


//================================================================================
class OutEdge
{
public:
	bool mHasStrictCut;
	Point3D mCutPoint;
	bool mAlreadyChecked;
};


//================================================================================
class Outliner
{
#pragma region "Construction & destruction"
public:
	Outliner();
	virtual ~Outliner();
	void Free();
#pragma endregion


#pragma region "Outlining"
public:
    // Warning:
    // Outliner keeps pointers on both objects Mesh and MeshTopology
    //
    // There is no need to Set an object again if the object has been transformed
    // but the topology has not changed.
    //
    void Set( const Mesh & pMesh, const MeshTopology & pTopo );
    void Set( const Mesh & pMesh );
    // If last point of outline == first point then the outline is closed, if not the outline is open
	//
	// Convention
	// - Empty outline: 0 points
	// - Non-empty outline: at least 2 points
    // - Closed outline: first point == last point
    //
	void Cut( const Plane3D & pCutPlane, List< List<Point3D> > & pOutlines, double pTolerance=1e-6 );
	
	// Meshes the outline/polygon
	// pOutline = closed outline = circular list of points = (P0, P1, P2, ..., PN, P0)
	static void MeshOutline( const List<Point3D> & pOutline, const Vector3D & pCutNormal, Mesh & pMesh );

//#define COMPILE_SHEWCHUK

#ifdef COMPILE_SHEWCHUK
	// Meshes the outline/polygon using Shewchuk's Triangle library
	// pOutline = closed outline = circular list of points = (P0, P1, P2, ..., PN, P0)
	// pCutNormal must be the outgoing CCW normal of the outline
	// pMaxTriArea is the maximal triangle area allowed in the output mesh
	//
	// Uses: https://www.cs.cmu.edu/~quake/triangle.switch.html
	static void MeshOutline_Shewchuk( const List<Point3D> & pOutline, const Vector3D & pCutNormal, double pMaxTriArea, Mesh & pMesh );
#endif

protected:
    void CutFromVertex( int pVerIdx, List<Point3D> & pOutline );
	void CutFromEdge( int pEdgeIdx, List<Point3D> & pOutline );
#pragma endregion


#pragma region "Tools"
public:
	static bool IsClosed( const List<Point3D> & pOutline );
    enum class Orientation { CW, CCW, Undefined };
	static Orientation GetOrientation( const List<Point3D> & pOutline, const Vector3D & pRef );
    static string OrientationToString( Orientation pOr );
	static double OutlineLength( const List<Point3D> & pOutline );
    static double SignedOutlineArea( const List<Point3D> & pOutline, const Vector3D * ppNormal=nullptr );
	static Point3D OutlineWeightedCentroid( const List<Point3D> & pOutline );
	static void LogOutlinesInfo( const List< List<Point3D> > & pOutlines, const std::string & pTag="Outlines" );
	static void OutlineIntersections( const List<Point3D> & pOutline, const Plane3D & pCutPlane, List<Point3D> & pInters, double pTolerance=1e-6 );
	static void OutlineIntersections( const List< List<Point3D> > & pOutlines, const Plane3D & pCutPlane, List< List<Point3D> > & pInters, double pTolerance=1e-6 );


    //
//    static List<Point3D> ConvexifyOutline( const List<Point3D> & pClosedOut, const Vector3D * ppNormal=nullptr );
    static List<Point3D> ConvexifyOutline( /*const*/ List<Point3D> /*&*/ pClosedOut, const Vector3D * ppNormal=nullptr );





    // Signed distance to closed outline. Last point repeats the first point in outline to mark closedness.
	// Returns the projection of the point on the outline.
	// Optionally a normal to the outline plane can be specified, which will save the trouble of doing a PCA of the points cloud.
	// The orientation of this optional vector is unimportant as in/out is determined by solid angle logic.
    static double SgnDistToClosedOut( const List<Point3D> & pOutline, const Point3D & pPt,
                                      Point3D & pProjP/*, Vector3D & pProjN*/, const Vector3D * ppNormal=nullptr );
    static double SgnDistToClosedOuts( const List< List<Point3D> > & pOutlines, const Point3D & pPt,
                                       Point3D & pProjP/*, Vector3D & pProjN*/, const Vector3D * ppNormal=nullptr );

	// Outline intersection and union
	enum class Operator { Union, Inter, Sub, AllInOne };
	class CSG_debug_info
	{
	public:
		void Free()
		{
            mNods.clear();
			mCuts.DelAll();
		
			mStrips_A.DelAll();
			mStrips_B.DelAll();

			mFilStr_A.DelAll();
			mFilStr_B.DelAll();
		}

	public:
        vector<Point3D> mNods;
        List<size_t> mCuts;

        List< List<size_t> > mStrips_A;
        List< List<size_t> > mStrips_B;

        List< List<size_t> > mFilStr_A;
        List< List<size_t> > mFilStr_B;
	};

	// Outlines must be closed (repeating last point) and lie on the XY plane.
	static void CSG( Operator pOp,
					 const List<Point3D> & pA, const List<Point3D> & pB,
					 List< List<Point3D> > & pC,
					 double pTolerance=1e-6,
					 CSG_debug_info * ppDebugInfo=nullptr );

    // Offsetting
    static void NormalOffsetOutline( List<Point3D> & pOutline, double pOffset );

    // Cropping
    //
    // *** TODO !!
    //
//    static void Crop( const List<Point3D> & pOut, List< List<Point3D> > & pTo, std::function<bool(const Point3D & pPt)> pKeep );
    static void Crop( const List<Point3D> & pOutline, const Plane3D & pCutPlane, List< List<Point3D> > & pRes, double pTol=1e-6 );

    // Resampling
    static Point3D FindCurviPoint( const List<Point3D> & pOutline, double pS );
    //
    // Uses curvilinear coord to resample outline. Works well if the outline has a smooth shape but will
    // produce uneven sampling if the outline has regions with strong curvature variations.
    // This method DOES NOT retain the initial nodes.
    static void UniformResample( const List<Point3D> & pOutline, size_t pNbSegs, List<Point3D> & pRes );
#pragma endregion


#pragma region "Data"
protected:
	const Mesh * mpMesh;
	const MeshTopology * mpTopo;
    MeshTopology mTopo; // If no topology has been provided

	vector<OutVer> mOutVers;
	vector<bool> mOutTrisChecked;
	vector<OutEdge> mOutEdges;
#pragma endregion
};
}
