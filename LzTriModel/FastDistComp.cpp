#include "FastDistComp.h"
#include "PatchComputer.h"


namespace LzTriModel
{
#pragma region "Construction / destruction"
//================================================================================
FastDistComp::FastDistComp()
{
}

//================================================================================
FastDistComp::FastDistComp( std::function<void(FastDistComp * pThis)> pInit )
{
    // Init
    pInit(this);
}

//================================================================================
FastDistComp::~FastDistComp()
{
	Free();
}

//================================================================================
void FastDistComp::Free()
{
    mPatches.clear();
}
#pragma endregion


#pragma region "Setting"
//================================================================================
void FastDistComp::Set( const Mesh & pMesh,
                        Patching pPatching,
                        bool pIgnoreNullNormals/*=false*/,
                        bool pIgnoreMultiedges/*=false*/ )
{
	MeshTopology lTopo;
	lTopo.Set( pMesh );

    Set( pMesh, lTopo, pPatching, pIgnoreNullNormals, pIgnoreMultiedges );
}

//================================================================================
void FastDistComp::Set( const Mesh & pMesh,
                        const MeshTopology & pTopo,
                        Patching pPatching,
                        bool pIgnoreNullNormals/*=false*/,
                        bool pIgnoreMultiedges/*=false*/ )
{
    // Check
    if( !pTopo.IsCompatible(pMesh) )
        LzLogException("", "Provided mesh and topology are not compatible!")

	// Check
	if( !pTopo.TopTris().size() )
        LzLogException("", "Cannot set topology: no triangles in mesh!")

	// Check
    if( !pIgnoreMultiedges && pTopo.MultiEdges().size() )
        LzLogException("", "Found "<<pTopo.MultiEdges().size()<<" multiedges in mesh! Signed distance computation may get funky...")

try
{
	// Free previous data
	Free();

	//------------------------------
	// Analyze patches
	//------------------------------

	if( pPatching == Patching::Split )
	{
		// Compute patches
		List<Mesh> lPatches;
        PatchComputer::SplitMesh( pMesh, pTopo, lPatches, true/* Remove unused vertices and normals */ );

		// Set patches
        mPatches.resize( lPatches.Count() );

        // Treat single patch case specifically to reuse available topology
        if( lPatches.Count() == 1 )
        {
            SetPatch( lPatches.GetHead(),
                      pTopo,
                      mPatches[0],
                      pIgnoreNullNormals,
                      pIgnoreMultiedges );
        }
        else
        {
            int iPatch = 0;
            BrowseList( iP, lPatches )
            {
                MeshTopology lPatchTopo;
                lPatchTopo.Set( lPatches.GetAt(iP) );

                SetPatch( lPatches.GetAt(iP),
                          lPatchTopo,
                          mPatches[iPatch],
                          pIgnoreNullNormals,
                          pIgnoreMultiedges );
                iPatch++;
            }
        }
	}
	else
	{
		// Exception will be raised if there are any unused vertices and normals
		// in the original mesh

		// Set a single patch
        mPatches.resize( 1 );
        SetPatch( pMesh,
                  pTopo,
                  mPatches[0],
                  pIgnoreNullNormals,
                  pIgnoreMultiedges );
	}
}
catch ( const std::exception & e )
{
    LzLogErr("", e.what())
    Free();
    LzLogException("", "Could not set distance computer!");
}
}

//================================================================================
void FastDistComp::SetPatch( const Mesh & pMesh,
                             const MeshTopology & pTopo,
                             Patch & pToPatch,
                             bool pIgnoreNullNormals,
                             bool pIgnoreMultiedges )
{
    // Chack if the patch has any unused vertices
    for( size_t v=0 ; v<pTopo.TopVers().size() ; v++ )
	{
		// Unused vertex?
		if( pTopo.TopVers(v).mT.Count() == 0 )
            LzLogException("", "Cannot set patch! Vertex "<<v<<" is not used by any triangle. Clean mesh before trying again.")
	}

	// Init tree 3D
	pToPatch.mTree.Create( pMesh.mVertices );

	// BBox
	pToPatch.mBBox.Reset();
	pToPatch.mBBox.Update( pMesh.mVertices );


	//-------------------------------------
#pragma region "Temp normals at triangles"
	//-------------------------------------

	// Patch data
    pToPatch.mTri2Vers.resize( pTopo.TopTris().size() );
    pToPatch.mTri2Edges.resize( pTopo.TopTris().size() );

	// Temp data
    Vector<Vector3D> lTriangleNormals( pTopo.TopTris().size() );
    for( size_t t=0 ; t<pTopo.TopTris().size() ; t++ )
	{
		const TopTri & lTri = pTopo.TopTris()[ t ];

		// Triangle vertices
		for( int i=0 ; i<3 ; i++ )
			pToPatch.mTri2Vers[t].mI[i] = lTri.mV[i];

		// Triangle edges
		for( int i=0 ; i<3 ; i++ )
			pToPatch.mTri2Edges[t].mI[i] = lTri.mE[i];

		// Triangle normals
		const Point3D & lA = pMesh.mVertices[ lTri.mV[0] ];
		const Point3D & lB = pMesh.mVertices[ lTri.mV[1] ];
		const Point3D & lC = pMesh.mVertices[ lTri.mV[2] ];

		lTriangleNormals[t] = (lB - lA)^(lC - lA);

		// Normalize normal
        _LzLogTry
			lTriangleNormals[t].Normalize();
        _LzLogCatchAndThrow("Could not compute triangle normal at triangle "<<t<<": A= "<<lA.ToString()<<", B= "<<lB.ToString()<<", C= "<<lC.ToString()
                            <<"! (patch has "<<pMesh.mVertices.size()<<" vertice(s) and "<<pMesh.mTriangles.size()<<" triangle(s))")
	}
#pragma endregion


	//-------------------------------------
#pragma region "Pseudo-normals at vertices"
	//-------------------------------------

	// Patch data
    pToPatch.mVer2Tris.resize( pTopo.TopVers().size() );
    pToPatch.mVerNors.resize( pTopo.TopVers().size() );
    for( size_t v=0 ; v<pTopo.TopVers().size() ; v++ )
	{
		// Vertices to triangles
		pToPatch.mVer2Tris[v] = pTopo.TopVers()[v].mT;
		
		//
		// All vertices in patch are used
		//

		// Vertex index, position & topology
		const Point3D & lVer = pMesh.mVertices[ v ];
		const TopVer & lTopVer = pTopo.TopVers()[ v ];

		// Reset pseudo normal
		pToPatch.mVerNors[v] = Vector3D(0,0,0);

        // Loop over adjacent triangles and accumulate
		BrowseList( iT, lTopVer.mT )
		{
			// Triangle index
            size_t t = lTopVer.mT.GetAt(iT);

			// Find the two opposite vertices in this triangle
			const TopEdge & lOppTopEdge = pTopo.TopEdges()[ pTopo.TopTris()[t].EdgeWithoutVertex(v) ];
			const Point3D & lOppV0 = pMesh.mVertices[ lOppTopEdge.mV[0] ];
			const Point3D & lOppV1 = pMesh.mVertices[ lOppTopEdge.mV[1] ];

			// Compute angle in [0,Pi] and update pseudo normal
			double lAngleRad = (lOppV0-lVer).ShortPositiveAngleTo( lOppV1-lVer );
			pToPatch.mVerNors[v] += lAngleRad * lTriangleNormals[ t ];
		}

		// Normalize pseudo-normal
        try
        {
			pToPatch.mVerNors[v].Normalize();
        }
        catch( ... )
        {
            // Should null vertex normals be silent or not?
            if( pIgnoreNullNormals )
            {
                // Log
                LzLogErr("", "Could not compute vertex pseudo-normal at vertex "<<v<<"! (patch has "<<pMesh.mVertices.size()<<" vertice(s) and "<<pMesh.mTriangles.size()<<" triangle(s))")

                // Set default vector
                pToPatch.mVerNors[v] = Vector3D(1,0,0);
            }
            else
            {
                // Do not ignore! Throw!
                LzLogException("", "Could not compute vertex pseudo-normal at vertex "<<v<<"! (patch has "<<pMesh.mVertices.size()<<" vertice(s) and "<<pMesh.mTriangles.size()<<" triangle(s))")
            }
        }
	}
#pragma endregion


	//-------------------------------------
#pragma region "Pseudo-normals at edges"
	//-------------------------------------

    pToPatch.mEdgNors.resize( pTopo.TopEdges().size() );
    for( size_t e=0 ; e<pToPatch.mEdgNors.size() ; e++ )
	{
		// Edge topology
		const TopEdge & lTopEdge = pTopo.TopEdges()[ e ];

		// Compute edge pseudo normal
		if( lTopEdge.mT.Count() == 1 )
		{
			// Crack edge
			pToPatch.mEdgNors[e] = lTriangleNormals[ lTopEdge.mT.GetHead() ];
		}
		else
		if( lTopEdge.mT.Count() == 2 )
		{
			// Manifold edge
			int lTri1 = lTopEdge.mT.GetHead();
			int lTri2 = lTopEdge.mT.GetTail();
			pToPatch.mEdgNors[e] = lTriangleNormals[lTri1] + lTriangleNormals[lTri2];

			// Normalize pseudo-normal
            try
            {
                pToPatch.mEdgNors[e].Normalize();
            }
            catch(...)
            {
                if( pIgnoreNullNormals )
                {
                    // Set a dummy normal
                    pToPatch.mEdgNors[e] = Vector3D(1, 0, 0);
                }
                else
                {
                    LzLogException("", "Could not compute edge pseudo-normal between triangles "<<lTri1<<" and "<<lTri2<<"! (patch has "
                                       <<pMesh.mVertices.size()<<" vertice(s) and "<<pMesh.mTriangles.size()<<" triangle(s))")
                }
            }
		}
		else
        if( pIgnoreMultiedges )
        {
            // Set a dummy normal
            pToPatch.mEdgNors[e] = Vector3D(1, 0, 0);
        }
        else
        {
            // Not ignoring multi-edges... crash
            LzLogException("", "Cannot compute pseudo-normal for edge "<<e<<"! Edge has "<<lTopEdge.mT.Count()<<" triangles.");
        }
#pragma endregion
	}
}
#pragma endregion


#pragma region "Computing"
//================================================================================
Point3D FastDistComp::GetVertex( size_t pPatchIdx, size_t pVerIdx ) const
{
	// Check
    if( pPatchIdx >= mPatches.size() )
        LzLogException("", "Invalid patch index= "<<pPatchIdx<<"! We only have "<<mPatches.size()<<" patch(es).")

	return mPatches[ pPatchIdx ].mTree.OldPoint3D( pVerIdx );
}

//================================================================================
Point3D FastDistComp::NearestVertex( Tree3D::NearestMode pTMode, const Point3D & pFrom ) const
{
	// Check
    if( !mPatches.size() )
        LzLogException("", "No patches in distance computer!");

	// Minimal distance for all patches
	double lMinDist = +std::numeric_limits<double>::max();
	Point3D lBestVer;

	// Check all patches
    for( size_t p=0 ; p<mPatches.size() ; p++ )
	{
		// Find nearest point in this patch
        const Patch & lPatch = mPatches[p];
        const Point3D & lPatchVer = lPatch.mTree.OldPoint3D( lPatch.mTree.GetNearestVertex(pFrom, pTMode) );

		// Have a better distance?
		double lDist = pFrom.DistanceTo( lPatchVer );
		if( lMinDist > lDist )
		{
			lMinDist = lDist;
			lBestVer = lPatchVer;
		}
	}

	return lBestVer;
}

//================================================================================
size_t FastDistComp::SgnDstInfo::CountTargetVers() const
{
    switch( mTarget )
    {
        case FastDistComp::SgnDstInfo::Target::A:
        case FastDistComp::SgnDstInfo::Target::B:
        case FastDistComp::SgnDstInfo::Target::C:
            return 1;

        case FastDistComp::SgnDstInfo::Target::AB:
        case FastDistComp::SgnDstInfo::Target::AC:
        case FastDistComp::SgnDstInfo::Target::BC:
            return 2;

        case FastDistComp::SgnDstInfo::Target::ABC:
            return 3;

        default:
            LzLogException("", "Unexpected target!")
    }
}

//================================================================================
void FastDistComp::SgnDstInfo::LogProj( const std::vector<Point3D> & pVers,
                                        const Point3D & lCheckPoint ) const
{
    Point3D lP(0, 0, 0);
    for( int w=0 ; w<3 ; w++ )
    for( int i=0 ; i<3 ; i++ )
    {
        lP.mV[i] += mW[w] * pVers[ mV[w] ].mV[i];
    }

    const double lErr = lP.DistanceTo(lCheckPoint);
    LzLogMsg("", "Delta: P= "<<lErr)

    if( lErr > 1e-12 )
        LzLogException("", "*** Found a big difference!");
}

//================================================================================
void FastDistComp::SgnDstInfo::Log( const string & pPrefix/*=""*/ ) const
{
    LzLogNode("", "SgnDstInfo "<<pPrefix)
    LzLogM("", "Patch idx= "<<mPatchIdx)
    LzLogM("", "V= "<<mV[0]<<"; "<<mV[1]<<"; "<<mV[2])
    LzLogM("", "W= "<<mW[0]<<"; "<<mW[1]<<"; "<<mW[2])
    LzLogM("", "Target= "<<static_cast<int>( mTarget ))
}

//================================================================================
static double SignedDistanceToTri( const Point3D & pP, Point3D & pProjP, Vector3D & pProjN,
                                   const Point3D & pA, const Point3D & pB, const Point3D & pC,
                                   const Vector3D & pNA, const Vector3D & pNB, const Vector3D & pNC,
                                   const Vector3D & pNAB, const Vector3D & pNAC, const Vector3D & pNBC,
                                   FastDistComp::SgnDstInfo * ppInfo=nullptr )
{
    // Output signed distance
    double lDist;

    // Triangle vectors
    const Vector3D AB = pB - pA;
    const Vector3D AC = pC - pA;

    // Normal vector pointing towards positive half-space, assuming (ABC) is CCW
    const Vector3D N = AB^AC;
    const double N2 = N * N;


    //------------------------------
#pragma region "Projection on plane (ABC)"
    //------------------------------

    // P1 = P projected on plane
    Point3D P1;
    Vector3D AP1;
    {
        // Projection of P on ABC plane
        Vector3D AP = pP - pA;
        double k = ( AP * N ) / N2;
        P1 = pP - k*N;

        // Check if P1 is in triangle
        AP1 = P1 - pA;
        double x = (AP1 ^  AC) * N / N2;
        double y = (AB  ^ AP1) * N / N2;

        // AP1 = x*AB + y*AC
        if( x>=0 && y>=0 && x+y<=1 )
        {
            // Projection params
            pProjP = P1;
            pProjN = N;
            pProjN.Normalize();

            // Signed distance computed using triangle normal
            double lAbsDist = pP.DistanceTo(P1);
            double lSgn = N * ( pP - P1 );
            lDist = lSgn<0 ? -lAbsDist : +lAbsDist ;

            // Fill info
            if( ppInfo )
            {
                ppInfo->mW[0] = 1.0 - x - y;
                ppInfo->mW[1] = x;
                ppInfo->mW[2] = y;

                ppInfo->mTarget = FastDistComp::SgnDstInfo::Target::ABC;
            }

            return lDist;
        }
    }
#pragma endregion


    //-------------------------------------------------------------------
    // P1 is out of the triangle : compute projections on triangle edges
    //-------------------------------------------------------------------

    // pProjN = normal at projection used to compute the sign


    //-------------------------------
#pragma region "Distance to segment AB"
    //-------------------------------
    {
        double AB2 = AB * AB;
        double t = ( AP1 * AB ) / AB2;

        if( t <= 0 )
        {
            pProjP = pA;
            pProjN = pNA;

            // Fill info
            if( ppInfo )
            {
                ppInfo->mW[0] = 1.0;
                ppInfo->mW[1] = 0.0;
                ppInfo->mW[2] = 0.0;

                ppInfo->mTarget = FastDistComp::SgnDstInfo::Target::A;
            }
        }
        else
        if( t >= 1 )
        {
            pProjP = pB;
            pProjN = pNB;

            // Fill info
            if( ppInfo )
            {
                ppInfo->mW[0] = 0.0;
                ppInfo->mW[1] = 1.0;
                ppInfo->mW[2] = 0.0;

                ppInfo->mTarget = FastDistComp::SgnDstInfo::Target::B;
            }
        }
        else
        {
            pProjP = pA + t * AB;

            // Pseudo normal of the edge AB
            pProjN = pNAB;

            // Fill info
            if( ppInfo )
            {
                ppInfo->mW[0] = 1.0 - t;
                ppInfo->mW[1] = t;
                ppInfo->mW[2] = 0.0;

                ppInfo->mTarget = FastDistComp::SgnDstInfo::Target::AB;
            }
        }

        // Set first candidate value: will be challenged by distances below
        lDist = pProjP.DistanceTo( pP );

        // Pass through: Sign computation at the end
    }
#pragma endregion


    // Sub-info
    //
    // Only weight and target information are updated here
    //
    double lInfoW[3];
    FastDistComp::SgnDstInfo::Target lInfoTarget;
    //
    Point3D P2;
    Vector3D P2N;


    //-------------------------------
#pragma region "Distance to segment AC"
    //-------------------------------
    {
        double AC2 = AC * AC;
        double t = ( AP1 * AC ) / AC2;

        if( t <= 0 )
        {
            P2 = pA;
            P2N = pNA;

            // Fill tmp info
            if( ppInfo )
            {
                lInfoW[0] = 1.0;
                lInfoW[1] = 0.0;
                lInfoW[2] = 0.0;

                lInfoTarget = FastDistComp::SgnDstInfo::Target::A;
            }
        }
        else
        if( t >= 1 )
        {
            P2 = pC;
            P2N = pNC;

            // Fill tmp info
            if( ppInfo )
            {
                lInfoW[0] = 0.0;
                lInfoW[1] = 0.0;
                lInfoW[2] = 1.0;

                lInfoTarget = FastDistComp::SgnDstInfo::Target::C;
            }
        }
        else
        {
            P2 = pA + t * AC;

            // Pseudo normal of the edge AC
            P2N = pNAC;

            // Fill tmp info
            if( ppInfo )
            {
                lInfoW[0] = 1.0 - t;
                lInfoW[1] = 0.0;
                lInfoW[2] = t;

                lInfoTarget = FastDistComp::SgnDstInfo::Target::AC;
            }
        }

        double d = P2.DistanceTo( pP );

        // Better ?
        if( d < lDist )
        {
            // Update
            pProjP = P2;
            pProjN = P2N;
            if( ppInfo )
            {
                for( int w=0 ; w<3 ; w++ )
                    ppInfo->mW[w] = lInfoW[w];

                ppInfo->mTarget = lInfoTarget;
            }

            lDist = d;
        }
    }
#pragma endregion


    //-------------------------------
#pragma region "Distance to segment BC"
    //-------------------------------
    {
        // Triangle vector
        Vector3D BC = pC - pB;

        double BC2 = BC * BC;
        double t = ( (P1-pB) * BC ) / BC2;

        if( t <= 0 )
        {
            P2 = pB;
            P2N = pNB;

            // Fill tmp info
            if( ppInfo )
            {
                lInfoW[0] = 0.0;
                lInfoW[1] = 1.0;
                lInfoW[2] = 0.0;

                lInfoTarget = FastDistComp::SgnDstInfo::Target::B;
            }
        }
        else
        if( t >= 1 )
        {
            P2 = pC;
            P2N = pNC;

            // Fill tmp info
            if( ppInfo )
            {
                lInfoW[0] = 0.0;
                lInfoW[1] = 0.0;
                lInfoW[2] = 1.0;

                lInfoTarget = FastDistComp::SgnDstInfo::Target::C;
            }
        }
        else
        {
            P2 = pB + t * BC;

            // Pseudo normal of the edge BC
            P2N = pNBC;

            // Fill tmp info
            if( ppInfo )
            {
                lInfoW[0] = 0.0;
                lInfoW[1] = 1.0 - t;
                lInfoW[2] = t;

                lInfoTarget = FastDistComp::SgnDstInfo::Target::BC;
            }
        }

        double d = P2.DistanceTo( pP );

        // Better ?
        if( d < lDist )
        {
            // Update
            pProjP = P2;
            pProjN = P2N;
            if( ppInfo )
            {
                for( int w=0 ; w<3 ; w++ )
                    ppInfo->mW[w] = lInfoW[w];

                ppInfo->mTarget = lInfoTarget;
            }

            lDist = d;
        }
    }
#pragma endregion


    //-------------------------------
    // Orientation
    //-------------------------------

    // Change distance sign
    if( pProjN * ( pP - pProjP ) < 0 )
        lDist = -lDist;

    return lDist;
}

//================================================================================
double FastDistComp::SignedDistance( Tree3D::NearestMode pTMode,
                                     const Point3D & pFrom,
                                     Point3D & pProjP,
                                     Vector3D & pProjN,
                                     SgnDstInfo * ppInfo/*=nullptr*/ ) const
{
	// Check
    if( !mPatches.size() )
        LzLogException("", "No patches in distance computer!");

	// Minimal distance for all patches
	double lMinSignedDistance = +std::numeric_limits<double>::max();
    for( size_t p=0 ; p<mPatches.size() ; p++ )
	{
        /*****************************************************************
         * The computation is optimized assuming that patches are closed
         * surfaces or at least "independent" surfaces in the sense that
         * the object is not represented by a stitching of the patches
         * but as their union.
         *
         * Multi-patch, CAD-type surfaces should be treated as Single
         * patch, i.e. one patch enclosing all connex components.
         *
         *****************************************************************/

		// The patch
        const Patch & lPatch = mPatches[p];

        // If have a negative min distance (we are already inside a patch)
        // AND point is outside the bbox of the current patch ==> skip current patch
        // ('false' flag for pIncludingBoundary to reject points even at the surface of another patch)
        if( lMinSignedDistance<0 && !lPatch.mBBox.PointIsIn(pFrom, false) )
			continue;

		// Compute projection on patch
        double lPatchSignedDist = 0.0;
		Point3D lPatchProjP(-666, -666, -666); // In case point is not initialized by subsequent computations
		Vector3D lPatchProjN(-666.666, -666.666, -666.666);
		//
		SgnDstInfo lPatchInfo;
		if( ppInfo )
            lPatchInfo.mPatchIdx = p;

        // Find the nearest projection point on the patch.
        {
			// Assuming that there is no surface crossing between the point and its projection,
			// then orientation will be given by considering the signed distance to the projected point.
			//
			double lPatchMinAbsDist = +std::numeric_limits<double>::max();
            const size_t lNearIdx = lPatch.mTree.GetNearestVertex( pFrom, pTMode );
            const List<size_t> & lTris = lPatch.mVer2Tris[ lNearIdx ];

			// Check
			if( lTris.Count() == 0 )
                LzLogException("", "Cannot compute signed distance to patch #"<<p<<"! Vertex "<<lNearIdx<<" is not used by any triangle. Clean mesh before trying again.")

//***** The approximation/optimization is HERE!
//***** The approximation/optimization is HERE!
//***** The approximation/optimization is HERE!

			// Find nearest projection in neighborhood of nearest vertex
			SgnDstInfo lWeightTargetInfo;
			BrowseList( iT, lTris )
			{
#pragma region "Triangle data"
				// Indices
                const size_t lTriIdx = lTris.GetAt(iT);
				//
				const TriIdx & lTriVers = lPatch.mTri2Vers[ lTriIdx ];
                const size_t I0 = lTriVers.mI[ 0 ];
                const size_t I1 = lTriVers.mI[ 1 ];
                const size_t I2 = lTriVers.mI[ 2 ];

				// Vertices
				const Point3D & lA = lPatch.mTree.OldPoint3D( I0 );
				const Point3D & lB = lPatch.mTree.OldPoint3D( I1 );
				const Point3D & lC = lPatch.mTree.OldPoint3D( I2 );

				// Vertex normals
				const Vector3D & lNorA = lPatch.mVerNors[ I0 ];
				const Vector3D & lNorB = lPatch.mVerNors[ I1 ];
				const Vector3D & lNorC = lPatch.mVerNors[ I2 ];
	
				// Edge normals
				const TriIdx & lTriEdgs = lPatch.mTri2Edges[ lTriIdx ];
				const Vector3D & lNor_AB = lPatch.mEdgNors[ lTriEdgs.mI[2] ];
				const Vector3D & lNor_AC = lPatch.mEdgNors[ lTriEdgs.mI[1] ];
				const Vector3D & lNor_BC = lPatch.mEdgNors[ lTriEdgs.mI[0] ];
#pragma endregion

				// Calc my distance and projection
				Point3D lProjP;
				Vector3D lProjN;
				double lMyDist = SignedDistanceToTri( pFrom, lProjP, lProjN,
													  lA, lB, lC,
													  lNorA, lNorB, lNorC,
													  lNor_AB, lNor_AC, lNor_BC,
                                                      ppInfo ? &lWeightTargetInfo/* Only weight and target information are updated here */
                                                             : nullptr );
                // Take absolute distance
				double lMyAbsDist = lMyDist<0 ? -lMyDist : +lMyDist ;

				// Found a nearest projection?
				if( lPatchMinAbsDist > lMyAbsDist )
				{
					// Keep projection
					lPatchMinAbsDist = lMyAbsDist;
					lPatchSignedDist = lMyDist;
					lPatchProjP = lProjP;
					lPatchProjN = lProjN;

					// Fill remaining info about projection
					if( ppInfo )
					{
						// Fill patch info
						lPatchInfo.mV[0] = I0;
						lPatchInfo.mV[1] = I1;
						lPatchInfo.mV[2] = I2;

						// Recover rest of the data from TriInfo
						for( int w=0 ; w<3 ; w++ )
							lPatchInfo.mW[w] = lWeightTargetInfo.mW[w];
						//
						lPatchInfo.mTarget = lWeightTargetInfo.mTarget;
					}
				}
			}
		}

		// Find the smallest distance for all patches
		//
		// SplitPatches define a UNION of domains. This computation gives an exact result if patches are not intersecting
		// each other (separate surfaces) but is approximates the distance if the surfaces do overlap. In this case a negative
		// distance is returned if the point is inside one of the surfaces but the value of the distance is not always correct.
		//
		if( lMinSignedDistance > lPatchSignedDist )
		{
			lMinSignedDistance = lPatchSignedDist;
			pProjP = lPatchProjP;
			pProjN = lPatchProjN;

			// Store interpolation information, if needed
			if( ppInfo )
				*ppInfo = lPatchInfo;
		}
	}

	return lMinSignedDistance;
}

//================================================================================
double FastDistComp::InOutDistance( Tree3D::NearestMode pTMode, const Point3D & pFrom ) const
{
	// Check
    if( !mPatches.size() )
        LzLogException("", "No patches in distance computer!");

    for( size_t iPatchIdx=0 ; iPatchIdx<mPatches.size() ; iPatchIdx++ )
	{
		// The patch
		const Patch & lPatch = mPatches[iPatchIdx];

		// If outside patch, keep current distance
		if( !lPatch.mBBox.PointIsIn(pFrom,true) )
			continue;

		// Compute projection on patch
		double lPatchMinAbsDist = +std::numeric_limits<double>::max();
		double lPatchSignedDist;
		{
			// Find the nearest projection point on the patch.
            const size_t lNearIdx = lPatch.mTree.GetNearestVertex( pFrom, pTMode );
            const List<size_t> & lTris = lPatch.mVer2Tris[ lNearIdx ];

			// Find nearest projection in neighborhood of nearest vertex
			BrowseList( iT, lTris )
			{
#pragma region "Triangle data"
				const int lTriIdx = lTris.GetAt(iT);

				// Vertices
				const TriIdx & lTriVers = lPatch.mTri2Vers[ lTriIdx ];
				const Point3D & lA = lPatch.mTree.OldPoint3D( lTriVers.mI[0] );
				const Point3D & lB = lPatch.mTree.OldPoint3D( lTriVers.mI[1] );
				const Point3D & lC = lPatch.mTree.OldPoint3D( lTriVers.mI[2] );

				// Vertex normals
				const Vector3D & lNorA = lPatch.mVerNors[ lTriVers.mI[0] ];
				const Vector3D & lNorB = lPatch.mVerNors[ lTriVers.mI[1] ];
				const Vector3D & lNorC = lPatch.mVerNors[ lTriVers.mI[2] ];
					
				// Edge normals
				const TriIdx & lTriEdgs = lPatch.mTri2Edges[ lTriIdx ];
				const Vector3D & lNor_AB = lPatch.mEdgNors[ lTriEdgs.mI[2] ];
				const Vector3D & lNor_AC = lPatch.mEdgNors[ lTriEdgs.mI[1] ];
				const Vector3D & lNor_BC = lPatch.mEdgNors[ lTriEdgs.mI[0] ];
#pragma endregion

				// Calc my distance and projection
                Vector3D tmp1;
                Point3D tmp2;
				double lMyDist = SignedDistanceToTri( pFrom, tmp2, tmp1,
                                                      lA, lB, lC,
                                                      lNorA, lNorB, lNorC,
                                                      lNor_AB, lNor_AC, lNor_BC );

				double lMyAbsDist = lMyDist<0 ? -lMyDist : +lMyDist ;

				// Found a nearest projection?
				if( lPatchMinAbsDist > lMyAbsDist )
				{
					lPatchMinAbsDist = lMyAbsDist;
					lPatchSignedDist = lMyDist;
				}
			}
		}

		// Inside a patch?
		if( lPatchSignedDist < 0 )
			return -1;
	}

	return +1;
}

//================================================================================
static double SignedDistanceToTri_NoProj( const Point3D & pP,
                                          const Point3D & pA, const Point3D & pB, const Point3D & pC,
                                          const Vector3D & pNA, const Vector3D & pNB, const Vector3D & pNC,
                                          const Vector3D & pNAB, const Vector3D & pNAC, const Vector3D & pNBC )
{
    // Output signed distance
    double lDist;

    // Triangle vectors
    const Vector3D AB = pB - pA;
    const Vector3D AC = pC - pA;

    // Normal vector pointing towards positive half-space, assuming (ABC) is CCW
    const Vector3D N = AB^AC;
    const double N2 = N * N;


    //------------------------------
#pragma region "Projection on plane (ABC)"
    //------------------------------

    // P1 = P projected on plane
    Point3D P1;
    Vector3D AP1;
    {
        // Projection of P on ABC plane
        Vector3D AP = pP - pA;
        double k = ( AP * N ) / N2;
        P1 = pP - k*N;

        // Check if P1 is in triangle
        AP1 = P1 - pA;
        double x = (AP1 ^  AC) * N / N2;
        double y = (AB  ^ AP1) * N / N2;

        // AP1 = x*AB + y*AC
        if( x>=0 && y>=0 && x+y<=1 )
        {
            // Signed distance computed using triangle normal
            double lAbsDist = pP.DistanceTo(P1);
            double lSgn = N * ( pP - P1 );
            lDist = lSgn<0 ? -lAbsDist : +lAbsDist ;

            return lDist;
        }
    }
#pragma endregion


    //-------------------------------------------------------------------
    // P1 is out of the triangle : compute projections on triangle edges
    //-------------------------------------------------------------------

    // pProjN = normal at projection used to compute the sign
    Point3D pProjP;
    Vector3D pProjN;


    //-------------------------------
#pragma region "Distance to segment AB"
    //-------------------------------
    {
        double AB2 = AB * AB;
        double t = ( AP1 * AB ) / AB2;

        if( t <= 0 )
        {
            pProjP = pA;
            pProjN = pNA;
        }
        else
        if( t >= 1 )
        {
            pProjP = pB;
            pProjN = pNB;
        }
        else
        {
            pProjP = pA + t * AB;

            // Pseudo normal of the edge AB
            pProjN = pNAB;
        }

        // Set first candidate value: will be challenged by distances below
        lDist = pProjP.DistanceTo( pP );

        // Pass through: Sign computation at the end
    }
#pragma endregion


    // Sub-info
    Point3D P2;
    Vector3D P2N;


    //-------------------------------
#pragma region "Distance to segment AC"
    //-------------------------------
    {
        double AC2 = AC * AC;
        double t = ( AP1 * AC ) / AC2;

        if( t <= 0 )
        {
            P2 = pA;
            P2N = pNA;
        }
        else
        if( t >= 1 )
        {
            P2 = pC;
            P2N = pNC;
        }
        else
        {
            P2 = pA + t * AC;

            // Pseudo normal of the edge AC
            P2N = pNAC;
        }

        double d = P2.DistanceTo( pP );

        // Better ?
        if( d < lDist )
        {
            // Update
            pProjP = P2;
            pProjN = P2N;

            lDist = d;
        }
    }
#pragma endregion


    //-------------------------------
#pragma region "Distance to segment BC"
    //-------------------------------
    {
        // Triangle vector
        Vector3D BC = pC - pB;

        double BC2 = BC * BC;
        double t = ( (P1-pB) * BC ) / BC2;

        if( t <= 0 )
        {
            P2 = pB;
            P2N = pNB;
        }
        else
        if( t >= 1 )
        {
            P2 = pC;
            P2N = pNC;
        }
        else
        {
            P2 = pB + t * BC;

            // Pseudo normal of the edge BC
            P2N = pNBC;
        }

        double d = P2.DistanceTo( pP );

        // Better ?
        if( d < lDist )
        {
            // Update
            pProjP = P2;
            pProjN = P2N;

            lDist = d;
        }
    }
#pragma endregion


    //-------------------------------
    // Orientation
    //-------------------------------

    // Change distance sign
    if( pProjN * ( pP - pProjP ) < 0 )
        lDist = -lDist;

    return lDist;
}

//================================================================================
bool FastDistComp::AccurateIsIn( const Point3D & pPt ) const
{
    // Check
    if( !mPatches.size() )
        LzLogException("", "No patches in distance computer!");

    // Use accurate Tree3 mode
    const LzGeom::Tree3D::NearestMode lTreeMode = LzGeom::Tree3D::NearestMode::Accurate;

    // Count patches
    const size_t lNbPatches = mPatches.size();

    //-------------------------
    // Compute 'in' flags
    //-------------------------

    // Is the point in the bbox of a patch?
    vector<bool> lIsInPatchBB;

    // Is the point in any bbox?
    bool lIsInAnyBB = false;

    // Compute flags
    lIsInPatchBB.resize(lNbPatches, false);
    for( size_t p=0 ; p<lNbPatches ; p++ )
    {
        // Include boundary so that points on the surface of a patch don't get rejected
        if( mPatches[p].mBBox.PointIsIn(pPt, /*pIncludingBoundary=*/true) )
        {
            lIsInPatchBB[p] = true;
            lIsInAnyBB = true;
        }
    }

    // Not contained in any patch
    if( !lIsInAnyBB )
        return false;

    //-------------------------
    // Point is in 1+ patch
    //-------------------------

    for( size_t p=0 ; p<lNbPatches ; p++ )
    {
        // Skip if not contained in a patch
        if( !lIsInPatchBB[p] )
            continue;

        // Get the patch
        const Patch & lPatch = mPatches[p];

        // Get the nearest triangles
        const size_t lNearIdx = lPatch.mTree.GetNearestVertex(pPt, lTreeMode);
        const List<size_t> & lTris = lPatch.mVer2Tris[ lNearIdx ];

        // Check
        if( lTris.Count() == 0 )
            LzLogException("", "Cannot compute signed distance to patch #"<<p<<"! Vertex "<<lNearIdx<<" is not used by any triangle. Clean mesh before trying again.")

        // Check the point's distance to this patch
        double lPatchMinAbsDist = +std::numeric_limits<double>::max();
        double lPatchSignedDist = 0.0;
        BrowseList( iT, lTris )
        {
#pragma region "Triangle data"
            // Indices
            const size_t lTriIdx = lTris.GetAt(iT);
            //
            const TriIdx & lTriVers = lPatch.mTri2Vers[ lTriIdx ];
            const size_t I0 = lTriVers.mI[ 0 ];
            const size_t I1 = lTriVers.mI[ 1 ];
            const size_t I2 = lTriVers.mI[ 2 ];

            // Vertices
            const Point3D & lA = lPatch.mTree.OldPoint3D( I0 );
            const Point3D & lB = lPatch.mTree.OldPoint3D( I1 );
            const Point3D & lC = lPatch.mTree.OldPoint3D( I2 );

            // Vertex normals
            const Vector3D & lNorA = lPatch.mVerNors[ I0 ];
            const Vector3D & lNorB = lPatch.mVerNors[ I1 ];
            const Vector3D & lNorC = lPatch.mVerNors[ I2 ];

            // Edge normals
            const TriIdx & lTriEdgs = lPatch.mTri2Edges[ lTriIdx ];
            const Vector3D & lNor_AB = lPatch.mEdgNors[ lTriEdgs.mI[2] ];
            const Vector3D & lNor_AC = lPatch.mEdgNors[ lTriEdgs.mI[1] ];
            const Vector3D & lNor_BC = lPatch.mEdgNors[ lTriEdgs.mI[0] ];
#pragma endregion

            // Calc my distance and projection
            double lMyDist = SignedDistanceToTri_NoProj( pPt,
                                                         lA, lB, lC,
                                                         lNorA, lNorB, lNorC,
                                                         lNor_AB, lNor_AC, lNor_BC );
            // Compute abs dist
            double lMyAbsDist = lMyDist<0 ? -lMyDist : +lMyDist ;

            // Found a nearest projection?
            if( lPatchMinAbsDist > lMyAbsDist )
            {
                // Keep projection
                lPatchMinAbsDist = lMyAbsDist;
                lPatchSignedDist = lMyDist;
            }
        }

        // Is it inside this patch?
        if( lPatchSignedDist <= 0.0 )
            return true;
    }

    // The point is not inside any patch
    return false;
}

//================================================================================
const Tree3D & FastDistComp::GetTreeForPatch( size_t pPatchIdx ) const
{
    // Check
    if( pPatchIdx >= mPatches.size() )
        LzLogException("", "Cannot access to patch #"<<pPatchIdx<<"! Patch count= "<<mPatches.size()<<".")

    // Return tree
    return mPatches[pPatchIdx].mTree;
}
#pragma endregion

}
