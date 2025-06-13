#pragma once

#include <LzServices/List.h>
#include <vector>


namespace LzServices
{
using std::vector;

class HashPos
{
public:
    unsigned int mIdx;
    void * mPos;
};

template <class T> class HashTable
{
#pragma region "Construction & destruction"
public:
    HashTable<T>( unsigned int pSize = 1000 );
    ~HashTable<T>();
    void SetSize( unsigned int pSize );
//protected:
    void Free();
#pragma endregion


#pragma region "Key generators"
public:
    static unsigned int HashKey( int pFromInteger );
    static unsigned int HashKey( const std::string & pFromString );
#pragma endregion


#pragma region "Usage"
public:
    void Insert( unsigned int pKey, T pVal )
	{
		if( mSize == 0 )
            LzLogException("", "Hash table size not set!")
	
		mpLists[pKey % mSize].AddTail( pVal );
	}
    const List<T> & ListOf( unsigned int pKey ) const
	{
		if( mSize == 0 )
            LzLogException("", "Hash table size not set!")
	
		return mpLists[pKey % mSize];
	}
    List<T> & ListOf( unsigned int pKey )
	{
		if( mSize == 0 )
            LzLogException("", "Hash table size not set!")
	
		return mpLists[pKey % mSize];
	}
    void FlushToVector( vector<T> & pToVector ) const;
    void DelAll();
#pragma endregion


#pragma region "Browsing"
public:
    HashPos HeadPos() const;
    void Next( HashPos & pPos ) const;
    T & GetAt( const HashPos & pPos );
    const T & GetAt( const HashPos & pPos ) const;
#define BrowseHashTable( iA, TableA ) for( LzServices::HashPos iA=(TableA).HeadPos() ; iA.mPos ; (TableA).Next(iA) )
#pragma endregion


#pragma region "Debug"
public:
    void Log( const std::string & pTag = "HashTable" ) const;
#pragma endregion


#pragma region "Data"
protected:
    unsigned int mSize;
    List<T> * mpLists;
#pragma endregion
};
}

#include "HashTable.inl"
