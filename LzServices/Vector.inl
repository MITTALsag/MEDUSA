#include <LzServices/LzLog.h>
#include <sstream>


namespace LzServices
{
using namespace std;


#pragma region "Construction & destruction"
//================================================================================
template<class T> Vector<T>::Vector() : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
}

//================================================================================
template<class T> Vector<T>::Vector( size_t pSize ) : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
    Resize( pSize );
}

//================================================================================
template<class T> Vector<T>::Vector( size_t pSize, const T & pT ) : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
    Assign( pSize, pT );
}

//================================================================================
template<class T> Vector<T>::Vector( size_t pSize, const T * ppT ) : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
    Resize( pSize );
    for( size_t i=0 ; i<pSize ; i++ )
        mpBuffer[i] = ppT[i];
}

//================================================================================
template<class T> Vector<T>::Vector( const Vector<T> & pOther ) : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
    *this = pOther;
}

//================================================================================
template<class T> Vector<T>::Vector( const std::vector<T> & pOther ) : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
    Resize( pOther.size() );
    for( size_t i=0 ; i<mSize ; i++ )
		mpBuffer[i] = pOther[i];
}

//================================================================================
template<class T> Vector<T>::Vector( const std::initializer_list<T> & pList ) : mSize( 0 ), mCapacity( 0 ), mpBuffer( nullptr )
{
    Resize( (size_t)pList.size() );

    size_t iNdex = 0;
	for( T iTem : pList )
		mpBuffer[iNdex++] = iTem;
}

//================================================================================
template<class T> Vector<T>::~Vector()
{
    Free();
}

//================================================================================
template<class T> void Vector<T>::Free()
{
    delete [] mpBuffer;
    mpBuffer = nullptr;
    mSize = 0;
    mCapacity = 0;
}
#pragma endregion


#pragma region "Operations"
//================================================================================
template<class T> void Vector<T>::Resize( size_t pNewSize )
{
    // Adjust capacity
    if( pNewSize > mCapacity )
        SetCapacity( pNewSize );

    // Commit size
    mSize = pNewSize;
}

//================================================================================
template<class T> void Vector<T>::Resize( size_t pNewSize, const T & pNewT )
{
    // Adjust capacity
    if( pNewSize > mCapacity )
        SetCapacity( pNewSize );

	// Fill with new Ts
    for( size_t i=mSize ; i<pNewSize ; i++ )
		mpBuffer[i] = pNewT;

    // Commit size
    mSize = pNewSize;
}

//================================================================================
template<class T> void Vector<T>::Assign( size_t pNewSize, const T & pT )
{
    // Reserve memory
    Resize( pNewSize );

    // Set all
    SetAll( pT );
}

//================================================================================
template<class T> void Vector<T>::SetAll( const T & pT )
{
    // Set all
    for( size_t i = 0 ; i < mSize ; i++ )
        mpBuffer[i] = pT;
}

//================================================================================
template<class T> void Vector<T>::PushBack( const T & pT )
{
    static const size_t sEXTRA_SPACE = 20;

    // Adjust capacity
    if( mSize + 1 > mCapacity )
        SetCapacity( mCapacity + sEXTRA_SPACE );

    // Add new element
    mpBuffer[ mSize ] = pT;
    mSize++;
}

//================================================================================
template<class T> void Vector<T>::PushBack( const Vector<T> & pTs )
{
#if 0
    for( size_t i = 0 ; i < pTs.Size() ; i++ )
        PushBack( pTs[i] );
#else
    if( mSize + pTs.mSize > mCapacity )
        SetCapacity( mSize + pTs.mSize );

    for( size_t i = 0 ; i < pTs.Size() ; i++ )
        mpBuffer[mSize + i] = pTs[i];

    mSize += pTs.mSize;
#endif
}

//================================================================================
template<class T> void Vector<T>::PopBack()
{
    if( mSize )
        mSize--;
}

//================================================================================
template<class T> const Vector<T> & Vector<T>::operator=( const Vector<T> & pOther )
{
    Resize( pOther.mSize );

    for( size_t i = 0 ; i < mSize ; i++ )
        mpBuffer[i] = pOther.mpBuffer[i];

    return *this;
}

//================================================================================
template<class T> double Vector<T>::operator*( const Vector<T> & pOther ) const
{
    // Check
    if( mSize==0 || mSize!=pOther.mSize )
        LzLogException("", "Cannot take contract (take scalar product of) vectors! My size= "<<mSize<<", other size= "<<pOther.mSize<<".")

    // Compute contraction (scalar product)
    double lScal = 0.0;
    for( size_t i=0 ; i<mSize ; i++ )
        lScal += mpBuffer[i] * pOther.mpBuffer[i];

    return lScal;
}

//================================================================================
//template<class T> const Vector<T> Vector<T>::operator-() const
//{
//
//    toudou si besoin
//
//}

//================================================================================
template<class T> T & Vector<T>::operator[]( size_t pIdx )
{
    if( pIdx >= mSize )
        LzLogException( "", "Cannot access element index= " << pIdx << " in Vector (" << mSize << " element(s))!" );

    return mpBuffer[pIdx];
}

//================================================================================
//template<class T> T::reference_wrapper Vector<T>::operator[]( size_t pIdx )
//{
//    if( pIdx >= mSize )
//        TxLogException( "Cannot access element index= "
//                        << pIdx
//                        << " in Vector ("
//                        << mSize
//                        << " element(s))!" );
//
//    return mpBuffer[pIdx];
//}

//================================================================================
template<class T> const T & Vector<T>::operator[]( size_t pIdx ) const
{
    if( pIdx >= mSize )
      LzLogException( "", "Cannot access element index " << pIdx << " in Vector (" << mSize << " element(s)!" );

    return mpBuffer[pIdx];
}

//================================================================================
template<class T> T & Vector<T>::GetBack()
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot get last element in empty Vector!" );

    return mpBuffer[mSize - 1];
}

//================================================================================
template<class T> const T & Vector<T>::GetBack() const
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot get last element in empty Vector!" );

    return mpBuffer[mSize - 1];
}

//================================================================================
template<class T> void Vector<T>::ToStdVector( vector<T> & pTo ) const
{
    if( !mSize )
    {
        pTo.clear();
    }
    else
    {
        pTo.resize( mSize );
        for( size_t i=0 ; i<mSize ; i++ )
            pTo[i] = mpBuffer[i];
    }
}

//================================================================================
template<class T> void Vector<T>::SetCapacity( size_t pNewCapacity )
{
    // Avoid useless reallocation
    if( pNewCapacity == mCapacity )
        return;

    // Alloc new space
    T * lpNewBuffer = nullptr;
    if( pNewCapacity )
        lpNewBuffer = new T [pNewCapacity];

    // Set new size and capacity
    mCapacity = pNewCapacity;
    mSize = mCapacity >= mSize ? mSize : mCapacity ;

    // Copy previous elements
    for( size_t i = 0 ; i < mSize ; i++ )
        lpNewBuffer[i] = mpBuffer[i];

    // Delete previous container
    delete [] mpBuffer;

    // Commit capacity
    mpBuffer = lpNewBuffer;
}
#pragma endregion


#pragma region "Arithmetic"
//================================================================================
template<class T> T Vector<T>::Min() const
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot take min of empty vector!" );

    // Find min
//T lMin = T::MaxValue;
T lMin = std::numeric_limits<T>::max();
    for( size_t i = 0 ; i < mSize ; i++ )
    {
        if( lMin > mpBuffer[i] )
            lMin = mpBuffer[i];
    }

    return lMin;
}

//================================================================================
template<class T> T Vector<T>::Max() const
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot take max of empty vector!" );

    // Find max
//T lMax = T::MinValue;
T lMax = std::numeric_limits<T>::lowest();
   for( size_t i = 0 ; i < mSize ; i++ )
    {
        if( lMax < mpBuffer[i] )
            lMax = mpBuffer[i];
    }

    return lMax;
}

//================================================================================
template<class T> double Vector<T>::Mean() const
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot take mean of empty vector!" );

	// Compute mean
	double lSum = 0.0;
    for( size_t i=0 ; i<mSize ; i++ )
		lSum += mpBuffer[i];

    return lSum / (double)mSize;
}

//================================================================================
template<class T> size_t Vector<T>::MinIdx() const
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot take min index of empty vector!" );

    // Find min index
//T lMin = T::MaxValue;
T lMin = std::numeric_limits<T>::max();
    size_t lMinIdx = 0;
    for( size_t i = 0 ; i < mSize ; i++ )
    {
        if( lMin > mpBuffer[i] )
        {
            lMin = mpBuffer[i];
            lMinIdx = i;
        }
    }

    return lMinIdx;
}

//================================================================================
template<class T> size_t Vector<T>::MaxIdx() const
{
    // Check
    if( !mSize )
        LzLogException( "", "Cannot take max index of empty vector!" );

    // Find max index
	T lMax = std::numeric_limits<T>::lowest();
    size_t lMaxIdx = 0;
    for( size_t i=0 ; i<mSize ; i++ )
    {
        if( lMax < mpBuffer[i] )
        {
            lMax = mpBuffer[i];
            lMaxIdx = i;
        }
    }

    return lMaxIdx;
}
#pragma endregion


#pragma region "Debug"
//================================================================================
template<class T> void Vector<T>::Log( const string & pPrefix/*=""*/ ) const
{
    LzLogMsg( "", pPrefix << "Size= " << mSize << " (capacity= " << mCapacity << ")" );
}

//================================================================================
template<class T> std::string Vector<T>::ToString0( const string & pPrefix/*=""*/ ) const
{
	std::stringstream lStr;
	lStr << pPrefix;
    for( size_t i=0 ; i<mSize ; i++ )
	{
		lStr << mpBuffer[i];

		if( i+1 < mSize )
			lStr << ", ";
	}

    return lStr.str();
}

//================================================================================
template <class T> std::ostream & operator<<( std::ostream & pOut, const Vector<T> & pVec )
{
	pOut << pVec.ToString0();
	return pOut;
}


//================================================================================
template<class T> std::string Vector<T>::ToString1( const string & pPrefix/*=""*/ ) const
{
	std::stringstream lStr;
	lStr << "---------------" << std::endl;
	lStr << pPrefix << "Size= " << mSize << " (capacity= " << mCapacity << ")" << std::endl;
    for( size_t i=0 ; i<mSize ; i++ )
		lStr << "[" << i << "]= " << mpBuffer[i] << std::endl;
	lStr << "---------------" << std::endl;

    return lStr.str();
}

//================================================================================
template<class T> std::string Vector<T>::ToString2( const string & pPrefix/*=""*/ ) const
{
	std::stringstream lStr;
    lStr << "---------------" << endl;
	lStr << pPrefix << "Size= " << mSize << " (capacity= " << mCapacity << ")" << endl;
    for( size_t i=0 ; i<mSize ; i++ )
		lStr << "[" << i << "]= " << mpBuffer[i].ToString() << endl;
	lStr << "---------------" << endl;

    return lStr.str();
}
#pragma endregion
}
