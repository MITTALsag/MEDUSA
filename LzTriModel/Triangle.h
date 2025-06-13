#pragma once

#include "LzLib_TriModel.h"
#include <string>


// DLL exports
namespace LzTriModel
{
    class DLL_EXPORT Triangle;
}


namespace LzTriModel
{
using std::string;


class Triangle
{
public:	
	// Construction / destruction
#define I_TO_SZ(I) static_cast<size_t>(I)
    //
	Triangle();
//    // Legacy constructors
//    Triangle( int pV0, int pV1, int pV2,
//              int pN0, int pN1, int pN2 ) : Triangle( I_TO_SZ(pV0), I_TO_SZ(pV1), I_TO_SZ(pV2),
//                                                      I_TO_SZ(pN0), I_TO_SZ(pN1), I_TO_SZ(pN2) ) {}
//    Triangle( int pV0, int pV1, int pV2,
//              int pT0, int pT1, int pT2,
//              int pN0, int pN1, int pN2 ) : Triangle( I_TO_SZ(pV0), I_TO_SZ(pV1), I_TO_SZ(pV2),
//                                                      I_TO_SZ(pT0), I_TO_SZ(pT1), I_TO_SZ(pT2),
//                                                      I_TO_SZ(pN0), I_TO_SZ(pN1), I_TO_SZ(pN2) ) {}
    // New kid on the block
    Triangle( size_t pV0, size_t pV1, size_t pV2,
              size_t pN0, size_t pN1, size_t pN2 );
    Triangle( size_t pV0, size_t pV1, size_t pV2,
              size_t pT0, size_t pT1, size_t pT2,
              size_t pN0, size_t pN1, size_t pN2 );
    virtual ~Triangle();

	// Check
	bool HasRepeatedVers() const;

	// Log
    string ToString() const;

	// Data
    size_t mIdxV[3];	// Vertices
    size_t mIdxT[3];	// Texture
    size_t mIdxN[3];	// Normals
};
}
