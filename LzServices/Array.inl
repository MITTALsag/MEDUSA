#include "LzLog.h"


namespace LzServices
{
#pragma region "Construction & destruction"
//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( const string & pName/*=""*/ )
 : mName(pName)
{}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( std::function<void(Array<T, IDX_T> * pThis)> pInit, const string & pName/*=""*/ )
 : mName(pName)
{
    // Init
    pInit(this);
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( const vector<size_t> & pSize, const string & pName/*=""*/ )
 : mName(pName)
{
    Resize( pSize );
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( const vector<size_t> & pSize, const T & pT, const string & pName/*=""*/ )
 : mName(pName)
{
    Resize( pSize );
    SetAll( pT );
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( const Array<T, IDX_T> & pOther, const string & pName/*=""*/ )
 : mName(pName)
{
    // Rename only if no explicit name
    if( mName == "" )
        mName = "Copy_of("+pOther.Name()+")_ctor";

    // Copy data
    Resize( pOther.mSize );
    for( size_t c=0 ; c<mCount ; c++ )
        mpBuffer[c] = pOther.mpBuffer[c];
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( Array<T, IDX_T> && pOther, const string & pName/*=""*/ )
 : mName(pName)
{
    /* https://stackoverflow.com/questions/3106110/what-is-move-semantics
     *
     * To summarize, the copy constructor makes a deep copy, because the source must remain untouched.
     * The move constructor, on the other hand, can just copy the pointer and then set the pointer in the source to null.
     * It is okay to "nullify" the source object in this manner, because the client has no way of inspecting the object again.
     */

    // Rename only if no explicit name
    if( mName == "" )
        mName = "Moved_from("+pOther.Name()+")_ctor";

    // Move data
    mpBuffer = pOther.mpBuffer;
    mCount = pOther.mCount;
    mSize = pOther.mSize;

    // Nullify r-value source object
    pOther.mpBuffer = nullptr;
    pOther.Free();
}

//================================================================================
template<class T>
size_t CountFromSize( const vector<T> & pSize )
{
    size_t lCount = 1;
    for( size_t d=0 ; d<pSize.size() ; d++ )
        lCount *= pSize[d];

    return lCount;
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( const vector<size_t> & pSize,
                        const std::initializer_list<T> & pFor_i_For_j, const string & pName/*=""*/ )
 : mName(pName)
{
    // Check initializers and size adequacy before allocating memory
    /*
     * Destructor of current class not called if exception thrown from constructor!
     * Use std::auto_ptr<T> to prevent leakage
     * https://stackoverflow.com/questions/188693/is-the-destructor-called-if-the-constructor-throws-an-exception
     */
    if( CountFromSize(pSize) != pFor_i_For_j.size() )
        LzLogException("", "Size mismatch between Array size and initializers! "<<CountFromSize(pSize)<<" != "<<pFor_i_For_j.size()<<".")

    // Alloc array
    Resize( pSize );

    // Assign values
    size_t i = 0;
    for( const T & iTem : pFor_i_For_j )
        mpBuffer[i++] = iTem;
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::Array( const vector<size_t> & pSize,
                        const vector<T> & pFor_i_For_j, const string & pName/*=""*/ )
 : mName(pName)
{
    // Check initializers and size adequacy before allocating memory
    /*
     * Destructor of current class not called if exception thrown from constructor!
     * Use std::auto_ptr<T> to prevent leakage
     * https://stackoverflow.com/questions/188693/is-the-destructor-called-if-the-constructor-throws-an-exception
     */
    if( CountFromSize(pSize) != pFor_i_For_j.size() )
        LzLogException("", "Size mismatch between Array size and initializers! "<<CountFromSize(pSize)<<" != "<<pFor_i_For_j.size()<<".")

    // Alloc array
    Resize( pSize );

    // Assign values
    size_t i = 0;
    for( const T & iTem : pFor_i_For_j )
        mpBuffer[i++] = iTem;
}

//================================================================================
template<class T, class IDX_T>
Array<T, IDX_T>::~Array()
{
    if( mCount )
        Free();
}

//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::Free()
{
    // Log
    LzLogTimeN("", "Array<T, IDX_T>::free( "<<Name()<<" ), "<<mCount<<" element(s).")

    // Release memory
    delete [] mpBuffer;

    // Mark as empty
    // Name remains unchanged
    mpBuffer = nullptr;
    mCount = 0;
    mSize.clear();
}
#pragma endregion


#pragma region "Operations"
//================================================================================
template<class T, class IDX_T>
const string & Array<T, IDX_T>::Name() const
{
    static const string sUnnamed = "<unnamed>";
    return mName.length() ? mName : sUnnamed;
}

//================================================================================
template<class T, class IDX_T>
size_t Array<T, IDX_T>::Size( size_t pD ) const
{
    // Check
    if( pD >= mSize.size() )
        LzLogException("", "Invalid dimension: "<<pD<<"! Array has "<<mSize.size()<<" dimensions.")

    return mSize[pD];
}

//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::Resize( const vector<size_t> & pNewSize )
{
    // Log
    LzLogN("", "Resizing '"<<Name()<<"'")

    // Empty: no size specified?
    if( pNewSize.size() == 0 )
    {
        Free();
        return;
    }

    // Count elements
    size_t lNewCount = CountFromSize( pNewSize );

    // Empty: one of the sizes is 0?
    if( lNewCount == 0 )
    {
        Free();
        return;
    }

    // Here: NewCount != 0

    // Alloc
    try
    {
        // Need to reallocate memory?
        if( lNewCount != mCount )
        {
            // NewCount is different, and a buffer was allocated: clear the buffer
            if( mCount )
                Free();
// Log
LzLogM("", "Array<T, IDX_T>::Resize: trying to allocate "<<lNewCount<<" element(s): "<<lNewCount*sizeof(T)<<" bytes.")

            // Here: mpBuffer == nullptr

            // Alloc new buffer
            mpBuffer = new T[ lNewCount ];

            // If alloc succeded, mpBuffer != nullptr
// Log
LzLogM("", "Array<T, IDX_T>::Resize:          allocated "<<lNewCount<<" element(s): "<<lNewCount*sizeof(T)<<" bytes.")

            // Set new count
            mCount = lNewCount;
        }

        // Set new size
        mSize = pNewSize;
    }
    catch( const std::exception & pExc )
    {
        Free();

        LzLogException("", "Could not allocate array of "<<lNewCount*sizeof( T )<<" bytes! Exception: '"<<pExc.what()<<"'.")
    }
    catch( ... )
    {
        Free();

        LzLogException("", "Could not allocate array of "<<lNewCount*sizeof( T )<<" bytes! Unknown exception.")
    }
}

//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::SetAll( const T & pT )
{
    for( size_t i=0 ; i<mCount ; i++ )
        mpBuffer[i] = pT;
}

//================================================================================
template<class T, class IDX_T>
const Array<T, IDX_T> & Array<T, IDX_T>::operator=( const Array<T, IDX_T> & pOther )
{            
    // Rename only if no explicit name
    if( mName == "" )
        mName = "Copy_of("+pOther.Name()+")_assign";

    // Copy data
    Resize( pOther.mSize );
    for( size_t i=0 ; i<mCount ; i++ )
        mpBuffer[i] = pOther.mpBuffer[i];

    return *this;
}

//================================================================================
template<class T, class IDX_T>
const Array<T, IDX_T> & Array<T, IDX_T>::operator=( Array<T, IDX_T> && pOther )
{
    /* https://stackoverflow.com/questions/3106110/what-is-move-semantics
     *
     * To summarize, the copy constructor makes a deep copy, because the source must remain untouched.
     * The move constructor, on the other hand, can just copy the pointer and then set the pointer in the source to null.
     * It is okay to "nullify" the source object in this manner, because the client has no way of inspecting the object again.
     */

    // Rename only if no explicit name
    if( mName == "" )
        mName = "Moved_from("+pOther.Name()+")_assign";

    // Move data
    mpBuffer = pOther.mpBuffer;
    mCount = pOther.mCount;
    mSize = pOther.mSize;

    // Nullify r-value source object
    pOther.mpBuffer = nullptr;
    pOther.Free();

    return *this;
}
#pragma endregion


#pragma region "Comparison"
//================================================================================
template<class T, class IDX_T>
bool Array<T, IDX_T>::operator==( const Array<T, IDX_T> & pOther ) const
{
    // One of the arrays is empty?
    if( IsEmpty() )
        return pOther.IsEmpty();

    // This Array is not empty

    // Compare dims (other dim = 0, it pOther is empty)
    if( Dims() != pOther.Dims() )
        return false;
    //
    for( size_t d=0 ; d<Dims() ; d++ )
    {
        if( Size(d) != pOther.Size(d) )
            return false;
    }

    // Compare elements
    for( size_t c=0 ; c<Count() ; c++ )
    {
        if( Buffer()[c] != pOther.Buffer()[c] )
            return false;
    }

    return true;
}
#pragma endregion


#pragma region "Iteration"
//================================================================================
template<class T, class IDX_T>
vector<IDX_T> Array<T, IDX_T>::InitIdx() const
{
    vector<IDX_T> lIdx;
    lIdx.assign( mSize.size(), 0 );
    return lIdx;
}

//================================================================================
template<class T, class IDX_T>
bool Array<T, IDX_T>::CheckIdx( const vector<IDX_T> & pIdx ) const
{
    // Check
    if( mSize.size() != pIdx.size() )
        LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: "<<pIdx.size()<<".")

    // Empty?
    if( !mSize.size() )
        return false;

    // Non-empty: check idx.
    for( size_t d=0 ; d<mSize.size() ; d++ )
    {
        if( pIdx[d] >= mSize[d] )
            return false;
    }

    // Everything ok
    return true;
}

//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::NextIdx( vector<IDX_T> & pIdx ) const
{
    // Check
    if( mSize.size() != pIdx.size() )
        LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: "<<pIdx.size()<<".")

    // Check
    if( !mSize.size() )
        LzLogException("", "Cannot iterate through an empty array!")

	// Possible to iterate over an empty array: returns false immediately.
    for( int d=mSize.size()-1 ; d>=0 ; d-- )
    {
        if( pIdx[d] < mSize[d] - 1 )
        {
            pIdx[d]++;
            return;
        }
        else
            pIdx[d] = 0;
    }

    // Overflow first index to mark out of bounds
    pIdx[0] = mSize[0];
}

//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::CounterToIdx( size_t pCounter, vector<IDX_T> & pIdx ) const
{
    // Check
    if( pCounter >= Count() )
        LzLogException("", "Counter overflow! Elements count in array= "<<Count()<<".")

    pIdx.resize( Dims() );
    size_t d = Dims() - 1;

    while( true )
    {
        // Compute index in dimension
        pIdx[d] = pCounter % mSize[d];

        // Processed last dim?
        if( d == 0 )
            break;

        // Pop dimension
        pCounter = (pCounter - pIdx[d]) /  mSize[d];

        // Decrease dim
        d--;
    }
}

//================================================================================
template<class T, class IDX_T>
vector<IDX_T> Array<T, IDX_T>::CounterToIdx( size_t pCounter ) const
{
    vector<IDX_T> lIdx;
    CounterToIdx( pCounter, lIdx );
    return lIdx;
}

//================================================================================
template<class T, class IDX_T>
size_t Array<T, IDX_T>::IdxToCounter( IDX_T pIdx0 ) const
{
    // Check
    if( mSize.size() != 1 ) LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: 1." )

    // Check
    if( pIdx0 >= mSize[0] ) LzLogException("", "Dimension 0, index "<<pIdx0<<" is out of range! Max= "<<(mSize[0] - 1)<<".")

    // Compute index = pI0
    return pIdx0;
}

//================================================================================
template<class T, class IDX_T>
size_t Array<T, IDX_T>::IdxToCounter( IDX_T pIdx0, IDX_T pIdx1 ) const
{
    // Check
    if( mSize.size() != 2 ) LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: 2." )

    // Check
    if( pIdx0 >= mSize[0] ) LzLogException("", "Dimension 0, index "<<pIdx0<<" is out of range! Max= "<<(mSize[0] - 1)<<".")
    if( pIdx1 >= mSize[1] ) LzLogException("", "Dimension 1, index "<<pIdx1<<" is out of range! Max= "<<(mSize[1] - 1)<<".")

    // Compute index = pI1 + mSize[1]*pI0
    return pIdx1 + mSize[1]*pIdx0;
}

//================================================================================
template<class T, class IDX_T>
size_t Array<T, IDX_T>::IdxToCounter( IDX_T pIdx0, IDX_T pIdx1, IDX_T pIdx2 ) const
{
    // Check
    if( mSize.size() != 3 ) LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: 3." )

    // Check
    if( pIdx0 >= mSize[0] ) LzLogException("", "Dimension 0, index "<<pIdx0<<" is out of range! Max= "<<(mSize[0] - 1)<<".")
    if( pIdx1 >= mSize[1] ) LzLogException("", "Dimension 1, index "<<pIdx1<<" is out of range! Max= "<<(mSize[1] - 1)<<".")
    if( pIdx2 >= mSize[2] ) LzLogException("", "Dimension 2, index "<<pIdx2<<" is out of range! Max= "<<(mSize[2] - 1)<<".")

    // Compute index = pI2 + mSize[2]*( pI1 + mSize[1]*pI0 )
    return pIdx2 + mSize[2]*( pIdx1 + mSize[1]*pIdx0 );
}

//================================================================================
template<class T, class IDX_T>
size_t Array<T, IDX_T>::IdxToCounter( const vector<IDX_T> & pIdx ) const
{
    // Check
    if( mSize.size() != pIdx.size() )
        LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: "<<pIdx.size()<<"." )

    // Compute index = pI2 + mSize[2]*( pI1 + mSize[1]*pI0 )
    size_t lIdx = 0;
    for( size_t d=0 ; d<mSize.size() ; d++ )
    {
        // Check
        if( pIdx[d] >= mSize[d] )
            LzLogException("", "Dimension "<<d<<", index "<<pIdx[d]<<" is out of range! Max= "<<(mSize[d] - 1)<<".")

        // Accumulate
        lIdx *= mSize[d];
        lIdx += pIdx[d];
    }

    return lIdx;
}

//================================================================================
template<class T, class IDX_T>
T & Array<T, IDX_T>::operator()( IDX_T pIdx0 )
{
    // Check
    if( mSize.size() != 1 ) LzLogException("", "Array has "<<mSize.size()<<" dimensions, not 1!")

    // Check
    if( pIdx0 >= mSize[0] ) LzLogException("", "Index "<<pIdx0<<" is out of range! Max= "<<( mSize[0] - 1 )<<".")

    return mpBuffer[pIdx0];
}

//================================================================================
template<class T, class IDX_T>
T & Array<T, IDX_T>::operator()( IDX_T pIdx0, IDX_T pIdx1 )
{
    // Check
    if( mSize.size() != 2 ) LzLogException("", "Array has "<<mSize.size()<<" dimensions, not 2!")

    // Check
    if( pIdx0 >= mSize[0] ) LzLogException("", "Index "<<pIdx0<<" is out of range! Max= "<<( mSize[0] - 1 )<<".")
    if( pIdx1 >= mSize[1] ) LzLogException("", "Index "<<pIdx1<<" is out of range! Max= "<<( mSize[1] - 1 )<<".")

    return mpBuffer[pIdx1 + mSize[1] * pIdx0];
}

//================================================================================
template<class T, class IDX_T>
T & Array<T, IDX_T>::operator()( IDX_T pIdx0, IDX_T pIdx1, IDX_T pIdx2 )
{
    // Check
    if( mSize.size() != 3 ) LzLogException("", "Array has "<<mSize.size()<<" dimensions, not 3!")

    // Check
    if( pIdx0 >= mSize[0] ) LzLogException("", "Index 0: "<<pIdx0<<" is out of range! Max= "<<( mSize[0] - 1 )<<".")
    if( pIdx1 >= mSize[1] ) LzLogException("", "Index 1: "<<pIdx1<<" is out of range! Max= "<<( mSize[1] - 1 )<<".")
    if( pIdx2 >= mSize[2] ) LzLogException("", "Index 2: "<<pIdx2<<" is out of range! Max= "<<( mSize[2] - 1 )<<".")

    return mpBuffer[pIdx2 + mSize[2] * ( pIdx1 + mSize[1] * pIdx0 )];
}

//================================================================================
template<class T, class IDX_T>
T & Array<T, IDX_T>::operator()( const vector<IDX_T> & pIdx )
{
    // Check
    if( !mSize.size() )
        LzLogException("", "Cannot access an element of an empty array!")

    return mpBuffer[ IdxToCounter(pIdx) ];
}

//================================================================================
template<class T, class IDX_T>
const T & Array<T, IDX_T>::operator()( IDX_T pIdx0 ) const
{
    return const_cast<Array<T, IDX_T>*>(this)->operator()( pIdx0 );
}

//================================================================================
template<class T, class IDX_T>
const T & Array<T, IDX_T>::operator()( IDX_T pIdx0, IDX_T pIdx1 ) const
{
    return const_cast<Array<T, IDX_T>*>(this)->operator()( pIdx0, pIdx1 );
}

//================================================================================
template<class T, class IDX_T>
const T & Array<T, IDX_T>::operator()( IDX_T pIdx0, IDX_T pIdx1, IDX_T pIdx2 ) const
{
    return const_cast<Array<T, IDX_T>*>(this)->operator()( pIdx0, pIdx1, pIdx2 );
}

//================================================================================
template<class T, class IDX_T>
const T & Array<T, IDX_T>::operator()( const vector<IDX_T> & pIdx ) const
{
    return const_cast<Array<T, IDX_T>*>(this)->operator()( pIdx );
}
#pragma endregion


#pragma region "Bounds"
//================================================================================
template<class T, class IDX_T>
T Array<T, IDX_T>::GetMin() const
{
    // Check
    if( !mCount )
        LzLogException("", "Cannot take min of empty array!")

    // Find min
    T lMin = +std::numeric_limits<T>::max();
    for( size_t i=0 ; i<mCount ; i++ )
    {
        if( lMin > mpBuffer[i] )
            lMin = mpBuffer[i];
    }

    return lMin;
}

//================================================================================
template<class T, class IDX_T>
T Array<T, IDX_T>::GetMax() const
{
    // Check
    if( !mCount )
        LzLogException("", "Cannot take max of empty array!")

    // Find max
    T lMax = -std::numeric_limits<T>::max();
    for( size_t i=0 ; i<mCount ; i++ )
    {
        if( lMax < mpBuffer[i] )
            lMax = mpBuffer[i];
    }

    return lMax;
}
#pragma endregion


#pragma region "Transformations"
//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::SubArrayTo( const vector<int> & pIdx, Array<T, IDX_T> & pTo ) const
{
    // Check
    if( mSize.size() != pIdx.size() )
        LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: "<<pIdx.size()<<".")

    // Count new dimensions
    vector<size_t> lToSize;
    for( size_t d=0 ; d<mSize.size() ; d++ )
    {
        // Free dim?
        if( pIdx[d] < 0 )
        {
            // Recover size at dim d
            lToSize.push_back( mSize[d] );
        }
        else
        // Imposed dim
        {
            // Check
            if( pIdx[d] >= mSize[d] )
                LzLogException("", "Imposed index "<<pIdx[d]<<", at dim "<<d<<", exceeds size= "<<mSize[d]<<"!")
        }
    }

    // Resize recipient
    pTo.Resize( lToSize );

    // Sub
    vector<IDX_T> lFromIdx = InitIdx();
    for( vector<IDX_T> iToIdx=pTo.InitIdx() ; pTo.CheckIdx(iToIdx) ; pTo.NextIdx(iToIdx) )
    {
        // Compute 'from' index
        int j = 0;
        for( size_t d=0 ; d<mSize.size() ; d++ )
        {
            // Free dim?
            if( pIdx[d] < 0 )
            {
                // Read index from iterator
                lFromIdx[d] = iToIdx[j++];
            }
            else
            // Imposed dim
            {
                // Read index from imposed indices
                lFromIdx[d] = static_cast<IDX_T>( pIdx[d] );
            }
        }

        // Copy value from *this to pTo
        pTo( iToIdx ) = (*this)( lFromIdx );
    }
}
#pragma endregion


#pragma region "Debug"
//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::Log( const string & pTag/*=""*/ ) const
{
    // Use provided tag or Array's name
    const string & lTag = pTag=="" ? Name() : pTag ;

    // Check
    if( !Count() )
    {
        LzLogM("", lTag<<": "<<IdxToString(mSize)<<" = EMPTY ARRAY" )
        return;
    }

    // Header
    LzLogN("", lTag<<" - Size= "<<IdxToString(mSize) )

    // Log elements
#if 1
    for( size_t c=0 ; c<Count() ; c++ )
        LzLogM("", lTag<<" ["<<c<<"] "<<IdxToString( CounterToIdx(c) )<<"= "<<mpBuffer[c])
#else
    for( vector<IDX_T> i=InitIdx() ; CheckIdx(i) ; NextIdx(i) )
        LzLogM("", lTag<<" "<<IdxToString(i)<<"= "<<mpBuffer[ IdxToCounter(i) ])
#endif
}

//================================================================================
template<class T, class IDX_T>
void Array<T, IDX_T>::LogInfo( const string & pTag/*=""*/ ) const
{
    // Use provided tag or Array's name
    const string & lTag = pTag=="" ? Name() : pTag ;

    // Log
    LzLogM("", "Array '"<<lTag<<"': "<<mSize.size()<<" dimension(s), size= "<<IdxToString( mSize )<<", "<<mCount<<" element(s).")
}

//================================================================================
template<class T, class IDX_T>
template<class I>
string Array<T, IDX_T>::IdxToString( const vector<I> & pIdx ) const
{
    // Check
    if( mSize.size() != pIdx.size() )
        LzLogException("", "Dimensions mismatch! Array: "<<mSize.size()<<", index: "<<pIdx.size()<<".")

    std::stringstream lStrStr;
    lStrStr << "[";
    for( size_t d=0 ; d<mSize.size() ; d++ )
    {
        lStrStr << " " << static_cast<size_t>( pIdx[d] );
        if( d < mSize.size() - 1 )
            lStrStr << ",";
    }
    lStrStr << " ]";

    return lStrStr.str();
}
#pragma endregion
}
