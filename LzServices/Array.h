#pragma once

#include <functional>
#include <vector>
#include <string>


namespace LzServices
{
using std::string;
using std::vector;


/**
 *
 * T is the type of the elements.
 * IDX_T is the type of the element coordinates (indices).
 *
 * Rationale: In some situations, where many element coordinates (vectors of indices) have to be stored,
 * memory space can be saved by downgrading IDX_T from full blown size_t down to a more economic format
 *
 * Element and dimensions counting uses size_t since
 * 1) using IDX_T would not provide significant RAM savings
 * 2) with narrow IDX_T, element count can easily overflow.
 *    Example: IDX_T = uchar (0--255); 2 dims; size(dim 0)=size(dim 1)=100 => count = 10,000: overflows IDX_T, aka uchar!
 *
 */
template<class T, class IDX_T=size_t> class Array
{
#pragma region "Construction & destruction"
public:
    Array( const string & pName="" );
    Array( std::function<void(Array<T, IDX_T> * pThis)> pInit, const string & pName="" );
    //
    Array( const vector<size_t> & pSize, const string & pName="" );
    Array( const vector<size_t> & pSize, const T & pT, const string & pName="" );
    //
    Array( const Array<T, IDX_T> & pOther, const string & pName="" );
    Array( Array<T, IDX_T> && pOther, const string & pName="" );
    //
    Array( const vector<size_t> & pSize,
           const std::initializer_list<T> & pFor_i_For_j, const string & pName="" );
    Array( const vector<size_t> & pSize,
           const vector<T> & pFor_i_For_j, const string & pName="" );
    //
    ~Array();
    void Free();
#pragma endregion


// You Ain't Gonna Need It
//#pragma region "Conservation & resurection"
//public:
//    void Save( const string & pFileName ) const;
//    void Load( const string & pFileName );
//#pragma endregion
// You Ain't Gonna Need It


#pragma region "Operations"
public:
    const string & Name() const;
    bool IsEmpty() const { return mCount == 0; }
    size_t Dims() const { return mSize.size(); }
    size_t Size( size_t pD ) const;
    const vector<size_t> & SizeVec() const { return mSize; }
    /*
     * Array::resize will not reallocate memory if resizing results in the
     * *same* elements count. In that case, the elements of Array will not
     * be reinitialized by resize however the contents might be rearranged
     * (if pNewSize provides a different layout of the elements);
     *
     */
    void Resize( const vector<size_t> & pNewSize );
    void SetAll( const T & pT );
    const Array<T, IDX_T> & operator=( const Array<T, IDX_T> & pOther );
    const Array<T, IDX_T> & operator=( Array<T, IDX_T> && pOther );
    size_t Count() const { return mCount; }
    //
    // Access by pointer to buffer
    //
    /**
     * @brief Buffer
     * @return copy of pointer to CONST elements
     */
    T const * Buffer() const { return mpBuffer; }
    /**
     * @brief Buffer
     * @return copy of pointer to NON-CONST elements
     */
    T * Buffer() { return mpBuffer; }
#pragma endregion


#pragma region "Comparison"
public:
    /**
     * @brief operator ==
     * @param pOther
     * @return true iff both arrays are empty OR both arrays are non-empty
     * AND they have the same sizes AND all elements are identical.
     * @warning 2 arrays can be == AND HAVE different names.
     */
    bool operator==( const Array<T, IDX_T> & pOther ) const;
    bool operator!=( const Array<T, IDX_T> & pOther ) const
    {
        return !( *this == pOther );
    }
#pragma endregion


#pragma region "Iteration"
public:
    // Iterator
    // Iterating order = for i, for j, ... (in 2D: row major)
    vector<IDX_T> InitIdx() const;
    bool CheckIdx( const vector<IDX_T> & pIdx ) const;
    void NextIdx( vector<IDX_T> & pIdx ) const;

    // Faster access (no vector construction)
    /*
     * ... to be implemented if needed
     *
     */
    // Slower generic access (via vector)
    /**
     * @brief CounterToIdx Implements the following memory mapping.
     *        index = pI2 + mSize[2]*( pI1 + mSize[1]*pI0 )
     *     or index = k + mSize[2]*( j + mSize[1]*i )
     *     which is the most optimal (from the cache point of view)
     *
     *     layout in memory for loops in the following order:
     *         for( i ... )
     *         for( j ... )
     *         for( k ... )
     *         {
     *             do_something_with_element( i, j, k )
     *         }
     *
     *     indeed, i is the least varying index, and k is the most
     *     varying, which ensures better locality.
     *
     * @param pCounter
     * @param pIdx
     */
    void CounterToIdx( size_t pCounter, vector<IDX_T> & pIdx ) const;
    vector<IDX_T> CounterToIdx( size_t pCounter ) const;

    // Faster access (no vector construction)
    size_t IdxToCounter( IDX_T pIdx0 ) const;
    size_t IdxToCounter( IDX_T pIdx0, IDX_T pIdx1 ) const;
    size_t IdxToCounter( IDX_T pIdx0, IDX_T pIdx1, IDX_T pIdx2 ) const;
    //
    // Slower generic access (via vector)
    size_t IdxToCounter( const vector<IDX_T> & pIdx ) const;

    /**
     *   /!\ NEVER keep references/pointers to Array elements: risk of reallocation!
     */

    // Faster access (no vector construction)
    T & operator()( IDX_T pIdx0 );
    T & operator()( IDX_T pIdx0, IDX_T pIdx1 );
    T & operator()( IDX_T pIdx0, IDX_T pIdx1, IDX_T pIdx2 );
	//
    // Slower generic access (via vector)
    T & operator()( const vector<IDX_T> & pIdx );
    //
    // Faster access (no vector construction)
    const T & operator()( IDX_T pIdx0 ) const;
    const T & operator()( IDX_T pIdx0, IDX_T pIdx1 ) const;
    const T & operator()( IDX_T pIdx0, IDX_T pIdx1, IDX_T pIdx2 ) const;
    //
    // Slower generic access (via vector)
    const T & operator()( const vector<IDX_T> & pIdx ) const;
#pragma endregion


#pragma region "Bounds"
public:
    /**
     * @brief getMin
     * @return Min value recomputed by going through all array elements.
     */
    T GetMin() const;

    /**
     * @brief getMax
     * @return Max value recomputed by going through all array elements.
     */
    T GetMax() const;
#pragma endregion


#pragma region "Transformations"
public:
    /**
     * @brief SubArrayTo extracts a sub array by projecting the source array onto imposed dims
     * @param pIdx: free dims are marked as -1 (any negative int), imposed dims must be specified
     * @param pTo
     *
     * E.g. to extract a 2d slice j=2 from a 3d array 3x4x5,
     *      pIdx = { -1, 2, -1 } will produce a 2d array 3x5
     */
    void SubArrayTo( const vector<int> & pIdx, Array<T, IDX_T> & pTo ) const;
#pragma endregion


#pragma region "Debug"
public:
    void Log( const string & pTag="" ) const;
    void LogInfo( const string & pTag="" ) const;
    template<class I> string IdxToString( const vector<I> & pIdx ) const;
#pragma endregion


#pragma region "Data"
protected:
    string mName;
    T * mpBuffer{ nullptr };
    size_t mCount{ 0 };
    vector<size_t> mSize;
#pragma endregion
};
}

#include "Array.inl"
