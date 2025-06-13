#pragma once

#include <string>
#include <vector>


namespace LzServices
{
using std::string;
using std::vector;


template<class T> class Vector
{
#pragma region "Construction & destruction"
public:
    Vector();
    Vector( size_t pSize );
    Vector( size_t pSize, const T & pT );
    Vector( size_t pSize, const T * ppT );
    Vector( const Vector<T> & pOther );
    Vector( const vector<T> & pOther );
    Vector( const std::initializer_list<T> & pList );
    ~Vector();
    void Free();
#pragma endregion


#pragma region "Operations"
public:
    size_t Size() const { return mSize; }
    void CropCapacityToSize() { SetCapacity( mSize ); }
	//
	// Copies content to the new memory location
    void Resize( size_t pNewSize );
	//
	// Copies content to the new memory location, extra slots init with pNewT
    void Resize( size_t pNewSize, const T & pNewT );
	//
    void Assign( size_t pNewSize, const T & pT );
    void SetAll( const T & pT );
    void PushBack( const T & pT );
    void PushBack( const Vector<T> & pTs );
    void PopBack();
    const Vector<T> & operator=( const Vector<T> & pOther );
    double operator*( const Vector<T> & pOther ) const;
//    const Vector<T> operator-() const;
    /***** /!\ NEVER keep references/pointers to Vector elements: risk of reallocation! *****/
    T & operator[]( size_t pIdx );
// hmm... comment ajouter move semantics pour pouvoir stocker des threads dans un Vector: T::reference_wrapper operator[]( size_t pIdx );
    const T & operator[]( size_t pIdx ) const;
    //
    T & GetBack();
    const T & GetBack() const;
    //
    T const * Buffer() const { return mpBuffer; }
    T * Buffer() { return mpBuffer; }
    /***** /!\ NEVER keep references/pointers to Vector elements: risk of reallocation! *****/

    ////--------------------------------------------------------------------------------
    ////--------------------------------------------------------------------------------
    //static void GetSubVector( size_t pStart, size_t pCount, const Vector<size_t> & pV, Vector<size_t> & pTo )
    //{
    //  // Check
    //  if( pStart>=pV.Size() || pStart+pCount>pV.Size() )
    //      TxLogException("Cannot extract sub-vector: start= "+pStart+", count="+pCount+"! Only "+pV.Size()+" element(s) in source vector.");
    //
    //  // Extract values
    //  pTo.Resize( pCount );
    //  for( size_t i=0 ; i<pCount ; i++ )
    //      pTo[i] = pV[pStart + i];
    //}

    // Convert to std::vector
    void ToStdVector( vector<T> & pTo ) const;

protected:
    void SetCapacity( size_t pNewCapacity );
#pragma endregion


#pragma region "Arithmetic"
public:
    T Min() const;
    T Max() const;
    double Mean() const;
    size_t MinIdx() const;
    size_t MaxIdx() const;
#pragma endregion


#pragma region "Debug"
public:
    void Log( const string & pPrefix="" ) const;
    string ToString0( const string & pPrefix="" ) const;
    string ToString1( const string & pPrefix="" ) const;
    string ToString2( const string & pPrefix="" ) const;
#pragma endregion


#pragma region "Data"
protected:
    size_t mSize;
    size_t mCapacity; // Must be always: mCapacity >= mSize
    T * mpBuffer;
#pragma endregion
};
}

#include "Vector.inl"
