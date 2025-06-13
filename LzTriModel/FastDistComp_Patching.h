#pragma once


namespace LzTriModel
{
//****************************************************************************
//
// Split: object=union of closed meshes. Distance computations are optimized.
// ------
//
// Single: we KNOW that the mesh contains one single closed (or opened patch)
// -------
//         OR
//
//         object=stitching of open meshes to form an approximate close
//         surface.
//
//
// /!\ Tree3D representation of Single surfaces can lead to approximations
//     near "stitches" as only a fraction of the "nearest triangle fan" can
//     be explored.
//
//****************************************************************************
enum class Patching { Split, Single };
}

