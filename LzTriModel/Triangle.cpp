#include "Triangle.h"
#include <LzServices/LzLog.h>
#include <LzMath/ToolBox.h>


namespace LzTriModel
{
//================================================================================
Triangle::Triangle()
{
	mIdxV[0] = mIdxV[1] = mIdxV[2] = -1; // Max unsigned int
	mIdxT[0] = mIdxT[1] = mIdxT[2] = -1; // Max unsigned int
	mIdxN[0] = mIdxN[1] = mIdxN[2] = -1; // Max unsigned int
}

//================================================================================
//Triangle::Triangle( int pV0, int pV1, int pV2,
//                    int pN0, int pN1, int pN2 )
Triangle::Triangle( size_t pV0, size_t pV1, size_t pV2,
                    size_t pN0, size_t pN1, size_t pN2 )
{
	mIdxV[0] = pV0;
	mIdxV[1] = pV1;
	mIdxV[2] = pV2;

    mIdxT[0] = -1;
    mIdxT[1] = -1;
    mIdxT[2] = -1;

    mIdxN[0] = pN0;
	mIdxN[1] = pN1;
	mIdxN[2] = pN2;
}

//================================================================================
//Triangle::Triangle( int pV0, int pV1, int pV2,
//                    int pT0, int pT1, int pT2,
//                    int pN0, int pN1, int pN2 )
Triangle::Triangle( size_t pV0, size_t pV1, size_t pV2,
                    size_t pT0, size_t pT1, size_t pT2,
                    size_t pN0, size_t pN1, size_t pN2 )
{
	mIdxV[0] = pV0;
	mIdxV[1] = pV1;
	mIdxV[2] = pV2;

	mIdxT[0] = pT0;
	mIdxT[1] = pT1;
	mIdxT[2] = pT2;

	mIdxN[0] = pN0;
	mIdxN[1] = pN1;
	mIdxN[2] = pN2;
}

//================================================================================
Triangle::~Triangle()
{
}

//================================================================================
bool Triangle::HasRepeatedVers() const
{
	if( mIdxV[0] == mIdxV[1]
	 || mIdxV[1] == mIdxV[2]
	 || mIdxV[2] == mIdxV[0] )
	{
		return true;
	}

	return false;
}

//================================================================================
std::string Triangle::ToString() const
{
	std::stringstream lStrStr;
	lStrStr << "Triangle( Ver=[" << mIdxV[0] << ", " << mIdxV[1] << ", " << mIdxV[2]
			<< "], Tex=[" << mIdxT[0] << ", " << mIdxT[1] << ", " << mIdxT[2]
			<< "], Nor=[" << mIdxN[0] << ", " << mIdxN[1] << ", " << mIdxN[2] << "] )";
	
	return lStrStr.str();
}

}
