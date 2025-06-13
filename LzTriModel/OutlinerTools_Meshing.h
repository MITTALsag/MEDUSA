#pragma once

#include "Mesh.h"
#include <map>



namespace LzTriModel
{
namespace OutlinerTools
{
/**
 * @brief The MeshOutlineOptions class
 */
class MeshOutlineOptions
{
    //-------------------
    // External options
    //-------------------
public:
    // Speed optimization: stop costly evaluation of splits as soon as
    // the min dist of a split is GREATER than the StopMinDist below
    // AND the split is an inner split (of course).
    bool mEarlyStopEvaluation{ true };
    double mEarlyStopMinDist{ 0.5 };

    // Verbosity
    bool mLog{ false };

    // Self intersections epsilon
    double mSelfInterEpsilon{ 1e-10 };

    //-------------------
    // Internal options
    //-------------------
public:
    // Debug recursive call id
    size_t mCallId{ 0 };

    // Does the outline have self-intersections?
    // This will condition the behavior when no "inner" splits can be found.
    bool mHasSelfInters;

    // Timing info
    string mJobName;
//**** Recup et comm' idees mail (fall back dist outline : il suffit de garder le critere de l'angle et ignorer la projection etc)

//**** scoring heuristics : above min, internal, near half split (time out if not found near half split) ... which links to order of evaluation
};


/**
 * @brief MeshOutline
 * @param pTo
 * @param pOut
 * @param pNormalCCW
 * @param pOptions
 */
void MeshOutline( Mesh & pTo,
                  const List<Point3D> & pOut,
                  const Vector3D & pNormalCCW,
                  const MeshOutlineOptions & pOptions=MeshOutlineOptions() );


/**
 * @brief Unit test helper functions
 */
void Test_MinSegToSeg();
void Test_HasSelfIntersections();
void Test_FlipFlopper();


/**
 * @brief Test_MeshOutline
 */
void Test_MeshOutline();

}
}
