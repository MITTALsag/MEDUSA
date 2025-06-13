#pragma once

#include <LzServices/Vector.h>
#include <LzServices/LzLog.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>


// MACROS
#define BrowseList( iA, ListA ) for( void *iA=(ListA).HeadPos() ; iA ; (ListA).Next(iA) )
#define BrowseTwoLists( iA, ListA, iB, ListB ) for( void *iA=(ListA).HeadPos(), *iB=(ListB).HeadPos() ; iA && iB ; (ListA).Next(iA), (ListB).Next(iB) )
#define BrowseThreeLists( iA, ListA, iB, ListB, iC, ListC ) for( void *iA=(ListA).HeadPos(), *iB=(ListB).HeadPos(), *iC=(ListC).HeadPos() ; iA && iB && iC ; (ListA).Next(iA), (ListB).Next(iB), (ListC).Next(iC) )
//
#define EsworbList( iA, ListA ) for( void *iA=(ListA).TailPos() ; iA ; (ListA).Prev(iA) )
#define EsworbTwoLists( iA, ListA, iB, ListB ) for( void *iA=(ListA).TailPos(), *iB=(ListB).TailPos() ; iA && iB ; (ListA).Prev(iA), (ListB).Prev(iB) )
#define EsworbThreeLists( iA, ListA, iB, ListB, iC, ListC ) for( void *iA=(ListA).TailPos(), *iB=(ListB).TailPos(), *iC=(ListC).TailPos() ; iA && iB && iC ; (ListA).Prev(iA), (ListB).Prev(iB), (ListC).Prev(iC) )


namespace LzServices
{
using std::vector;


//==============================================================================
template <class T, class Score> class Sortable
{
public:
    Sortable( T pT={}, Score pScore={} )
    : mT(pT), mScore(pScore)
	{}

	// Operators for insertion in sorted list
    bool operator>=( const Sortable<T,Score> & pOther ) { return mScore >= pOther.mScore; }
    bool operator<=( const Sortable<T,Score> & pOther ) { return mScore <= pOther.mScore; }
    bool operator!=( const Sortable<T,Score> & pOther ) { return (mScore != pOther.mScore)
                                                              || (mT != pOther.mT); }
	// Data
	T mT;
	Score mScore;
};
//
template <class T, class Score> std::ostream & operator<<( std::ostream & pOut, const Sortable<T,Score> & pSort )
{
    pOut << pSort.mT << " (Score= " << pSort.mScore << ")";
	return pOut;
}


//==============================================================================
template <class T> class List
{
#pragma region "Construction & destruction"
public:
    List<T>()
    : m_pHead( nullptr ), m_pTail( nullptr ), m_Count( 0 )
    {}

    List<T>( const std::initializer_list<T> & pList )
    : m_pHead( nullptr ), m_pTail( nullptr ), m_Count( 0 )
	{
		// Check
		if( pList.size() == 0 )
			return;

		// Stash
		for( auto & iT : pList )
			AddTail( iT );
	}

    List<T>( const List<T> & pOther )
    : m_pHead( nullptr ), m_pTail( nullptr ), m_Count( 0 )
    {
        *this = pOther;
    }

    List<T>( List<T> && pOther )
    : m_pHead( nullptr ), m_pTail( nullptr ), m_Count( 0 )
    {
        *this = std::move(pOther);
    }

    ~List<T>() { DelAll(); }

    // Copy and move assign
    List<T> & operator=( const List<T> & pL );
    List<T> & operator=( List<T> && pL );
#pragma endregion


#pragma region "Iteration"
public:
    void * HeadPos() const { return m_pHead; }
    void * TailPos() const { return m_pTail; }
    //
    void * NextPos( const void * const pPtr ) const { return pPtr ? (( ListElt<T>* )pPtr)->m_pNext : nullptr ; }
    void * PrevPos( const void * const pPtr ) const { return pPtr ? (( ListElt<T>* )pPtr)->m_pPrev : nullptr ; }
    //
    void * CircularNextPos( const void * const pPtr ) const { void * lPos = NextPos(pPtr); return lPos ? lPos : m_pHead ; } //recoder avec la logique de if nullptr then nullptr et checker qui utilise in fine ces trucs
    void * CircularPrevPos( const void * const pPtr ) const { void * lPos = PrevPos(pPtr); return lPos ? lPos : m_pTail ; }
    //
    void Next( void * & pPtr ) const { if( pPtr ) pPtr = (( ListElt<T>* )pPtr)->m_pNext; }
    void Prev( void * & pPtr ) const { if( pPtr ) pPtr = (( ListElt<T>* )pPtr)->m_pPrev; }
    size_t Count() const { return m_Count; }
#pragma endregion


#pragma region "Access"
public:
    // Access refs
    T & GetAt( void * pPtr );
    T & GetAtAndNext( void * & pPtr );
    T & GetTail() { return GetAt( m_pTail ); }
    T & GetHead() { return GetAt( m_pHead ); }
	T & GetAfterHead() { return GetAt( NextPos(m_pHead) ); }
	T & GetAfterAfterHead() { return GetAt( NextPos(NextPos(m_pHead)) ); }

    // Access const refs
    const T & GetAt( const void * pPtr ) const;
    const T & GetAtAndNext( /*const*/ void * & pPtr ) const;
    const T & GetAtAndPrev( /*const*/ void * & pPtr ) const;
/*
const char c = 'x';    // 1
char *p1;			   // 2
const char **p2 = &p1; // 3
*p2 = &c;			   // 4
*p1 = 'X';			   // 5

In line 3, we assign a char ** to a const char **. (The compiler should complain.)
In line 4, we assign a const char * to a const char *; this is clearly legal.
In line 5, we modify what a char * points to--this is supposed to be legal.
However, p1 ends up pointing to c, which is const.
This came about in line 4, because *p2 was really p1.
This was set up in line 3, which is an assignment of a form that is disallowed,
and this is exactly why line 3 is disallowed.
*/
    const T & GetTail() const { return GetAt( m_pTail ); }
    const T & GetHead() const { return GetAt( m_pHead ); }
	const T & GetAfterHead() const { return GetAt( NextPos(m_pHead) ); }
	const T & GetAfterAfterHead() const { return GetAt( NextPos(NextPos(m_pHead)) ); }
#pragma endregion


#pragma region "Search"
public:
    void * FindIdx( size_t pIdx ) const;
    void * FindElt( const T & pT ) const;
    void * FetchPos( void * pStartPos, int pCount );
    T & operator[]( size_t pIdx );
    const T & operator[]( size_t pIdx ) const;

	// Sorted lists
	void * FindInIncrList( const T & pT ) const;
#pragma endregion


#pragma region "Operations"
public:
    void * AddBefore( void * pPtr, const T & pT );
    void * AddBefore( void * pPtr, T && pT );
    //
    void * AddAfter( void * pPtr, const T & pT );
    void * AddAfter( void * pPtr, T && pT );
    //
    void * AddHead( const T & pT );
    void * AddHead( T && pT );
    //
    void * AddTail( const T & pT );
    void * AddTail( T && pT );
    //
    void * Append( const List<T> & pL );
    void * Append( List<T> && pL );
    //
    void * Prepend( const List<T> & pL );
    void * Prepend( const List<T> && pL )
    { LzLogException("", "*** NOT IMPLEMENTED") }
    //
	enum class AddIntoReturnPolicy { NullIfNoAdd, PosOfIdentical };
    void * AddIntoIncrList( const T & pT, bool pUnique, AddIntoReturnPolicy pPolicy=AddIntoReturnPolicy::NullIfNoAdd );
    void * AddIntoDecrList( const T & pT, bool pUnique, AddIntoReturnPolicy pPolicy=AddIntoReturnPolicy::NullIfNoAdd  );
	//
    void AddIntoIncrList( const List<T> & pL, bool pUnique );              // This list will be sorted, pL can be unsorted
    void AddIntoDecrList( const List<T> & pL, bool pUnique );              // This list will be sorted, pL can be unsorted
    void MergeWithIncrList( const List<T> & pL, bool pUniqueFromBtoA );    // Both lists have to be sorted
    void MergeUniqueIncrLists( const List<T> & pL );                       // Both lists have to be sorted and have unique elements. Uniqueness is preserved.
    //
    enum class SubstractMode { SubstractAll, SubstractCount };
    void SubstractIncrLists( const List<T> & pOther, SubstractMode pMode ); // Both lists have to be sorted
    //
    List<T> InterIncrLists( const List<T> & pOther, bool pUniqueInInter ) const;    // Both lists have to be sorted
    T Min() const;
    T Max() const;
	double Sum() const;
	double Mean() const;
	//
	void Reverse();
    void RotateToHead( void * pNewHeadPos ); // pNewHeadPos must be a valid position in the list otherwise an infinite loop occurs
#pragma endregion


#pragma region "Moving list elements"
	void * DetachAt( void * iPos )
	{
        LzLogException("", "*** NOT IMPLEMENTED")
//        return nullptr;
	}

	void AttachAfter( void * iPos )
	{
        LzLogException("", "*** NOT IMPLEMENTED")
	}

	void AttachBefore( void * iPos )
	{
        LzLogException("", "*** NOT IMPLEMENTED")
    }
#pragma endregion


#pragma region "Comparison"
public:
	// Two lists are equal iff they contain the exact same elements in the exact same order. A list is not considered as a 'set' as it is ordered.
    bool operator==( const List<T> & pL ) const; // /!\ Careful with repeated elements. Equality between lists can be tricky.
    bool operator!=( const List<T> & pL ) const;
    bool IsDirectShiftOf( const List<T> & pL ) const;
    bool IsReverseShiftOf( const List<T> & pL ) const;
	bool Has( const T & pT ) const;
#pragma endregion


#pragma region "Conversion to/from container"
public:
    void ToArray( T * pArray ) const;
    void FromArray( const T * pArray, size_t pCount );
    void FromArrayReverse( const T * pArray, size_t pCount );
    //
    void ToVector( vector<T> & pVector ) const;
    void ToVector( Vector<T> & pVector ) const;
    void FromVector( const vector<T> & pVector );
    void FromVector( const Vector<T> & pVector );
    // See also definitions below:
    //   template <class T> cli::array<T> ^ListToArray( const List<T> & pLst )
    //   template <class T, int N> void FromArray( List<T> & pLst, cli::array<T,N> ^pArr )
#pragma endregion


#pragma region "Deletion"
public:
    void DelAll();
    void DelAt( void * Ptr );
    void DelFromToTail( void * pPtr );
    void DelElts( const T & p_T );
    void DelAtAndNext( void * & Ptr );
    void DelHead() { DelAt( m_pHead ); }
    void DelTail() { DelAt( m_pTail );  }
    void FindAndDel( const T & pT, bool pExcIfNotFound );
#pragma endregion


#pragma region "Debug"
public:
    std::string ListToString() const;
    std::string ListOfListToString() const;
#pragma endregion


#pragma region "Data"
protected:
    //..............................................................................
    template <class T2> class ListElt
    {
    public:
        ListElt( const T2 & pT )
        : m_T( pT ), m_pNext( nullptr ), m_pPrev( nullptr )
        {}

        ListElt( T2 && pT )
        : m_T( std::move(pT) ), m_pNext( nullptr ), m_pPrev( nullptr )
        {}

        T2 m_T;

        ListElt * m_pNext;
        ListElt * m_pPrev;
    };
    //..............................................................................

    ListElt<T> * m_pHead;
    ListElt<T> * m_pTail;
    size_t m_Count;
#pragma endregion
};


#pragma region "Construction & destruction"
//==============================================================================
template <class T> List<T> & List<T>::operator=( const List<T> & pL )
{
    // Check
    if( &pL == this )
        return *this;

    // Assign
    DelAll();
    for( void * Pos = pL.HeadPos() ; Pos ; pL.Next( Pos ) )
        AddTail( pL.GetAt( Pos ) );

    return *this;
}

//==============================================================================
template <class T> List<T> & List<T>::operator=( List<T> && pL )
{
//// Check that actually moving stuff
//LzLogM("", "We got to move these refrigerators, we gotta move these color TVs!")

    // Check
    if( &pL == this )
        return *this;

    // Assign
    DelAll();

    // Steal data from other list
    m_pHead = pL.m_pHead;
    pL.m_pHead = nullptr;
    //
    m_pTail = pL.m_pTail;
    pL.m_pTail = nullptr;
    //
    m_Count = pL.m_Count;
    pL.m_Count = 0;

    return *this;
}
#pragma endregion


#pragma region "Access"
//==============================================================================
template <class T> T & List<T>::GetAt( void * Ptr )
{
    if( Ptr )
        return static_cast<ListElt<T>*>(Ptr)->m_T;
    else
        LzLogException("", "Accessing null position pointer!")
}

//==============================================================================
template <class T> T & List<T>::GetAtAndNext( void * & Ptr )
{
    if( Ptr )
    {
        T & lT = ((ListElt<T> *)Ptr )->m_T;
        Ptr = ((ListElt<T> *)Ptr)->m_pNext;
        return lT;
    }
    else
        LzLogException("", "Accessing null position pointer!")
}

//==============================================================================
template <class T> const T & List<T>::GetAt( const void * Ptr ) const
{
    if( Ptr )
        return static_cast<const ListElt<T> *>(Ptr)->m_T;
    else
        LzLogException("", "Accessing null position pointer!")
}

//==============================================================================
template <class T> const T & List<T>::GetAtAndNext( /*const*/ void * & pPtr ) const
{
    if( pPtr )
    {
        const T & lT = static_cast<const ListElt<T> *>(pPtr)->m_T;
        pPtr         = static_cast<const ListElt<T> *>(pPtr)->m_pNext;
        return lT;
    }
    else
        LzLogException("", "Accessing null position pointer!")
}

//==============================================================================
template <class T> const T & List<T>::GetAtAndPrev( /*const*/ void * & pPtr ) const
{
    if( pPtr )
    {
        const T & lT = static_cast<const ListElt<T> *>(pPtr)->m_T;
        pPtr         = static_cast<const ListElt<T> *>(pPtr)->m_pPrev;
        return lT;
    }
    else
        LzLogException("", "Accessing null position pointer!")
}
#pragma endregion


#pragma region "Search"
//==============================================================================
template <class T> void * List<T>::FindIdx( size_t Idx ) const
{
    // Geeky integer division by 2: m_Count >> 1
    // Start from head if idx < count / 2
    // Otherwise start from tail
    if( Idx < (m_Count >> 1) )
    {
        size_t I = 0;

        // From start
        for( void * Pos=m_pHead ; Pos ; Next(Pos) )
        {
            if( I == Idx )
                return Pos;

            I++;
        }
    }
    else
    {
        size_t I = m_Count - 1;

        // From end
        for( void * Pos=m_pTail ; Pos ; Prev(Pos) )
        {
            if( I == Idx )
                return Pos;

            I--;
        }
    }

    return nullptr;
}

//==============================================================================
template <class T> void * List<T>::FindElt( const T & pT ) const
{
    for( void * iT = m_pHead ; iT ; Next( iT ) )
    {
        if( GetAt( iT ) == pT )
            return iT;
    }

    return nullptr;
}

//==============================================================================
template <class T> void * List<T>::FetchPos( void * pStartPos, int pCount )
{
    if( !pStartPos )
        return nullptr;

    if( pCount < 0 )
    {
        while( pCount )
        {
            pStartPos = static_cast<ListElt<T>*>(pStartPos)->mpPrev;
            if( !pStartPos )
                return nullptr;

            pCount++;
        }
    }
    else
    if( pCount > 0 )
    {
        while( pCount )
        {
            pStartPos = static_cast<ListElt<T>*>(pStartPos)->mpNext;
            if( !pStartPos )
                return nullptr;

            pCount--;
        }
    }

    return pStartPos;
}

//==============================================================================
template <class T> T & List<T>::operator[]( size_t pIdx )
{
    void * lPos = FindIdx( pIdx );
    if( lPos )
        return GetAt( lPos );

    LzLogException("", "Invalid index in operator[" << pIdx << "]! Count= "<<m_Count<<".")
}

//==============================================================================
template <class T> const T & List<T>::operator[]( size_t pIdx ) const
{
//LzLogException("", "********* BREAKPOINT")





    void * lPos = FindIdx( pIdx );
    if( lPos )
        return GetAt( lPos );
    else

    LzLogException( "", "Invalid index in operator[" << pIdx << "]" )
}

//==============================================================================
template <class T> void * List<T>::FindInIncrList( const T & pT ) const
{
    BrowseList( i, *this )
    {
        const T & lT = GetAt(i);

        // Found a bigger element
        if( lT > pT )
            return nullptr;

        // Found same element
        if( lT == pT )
            return i;
    }

    // Everyone in the list is smaller and different
    return nullptr;
}
#pragma endregion


#pragma region "Operations"
//==============================================================================
template <class T> void * List<T>::AddBefore( void * Ptr, const T & p_T )
{
    if( Ptr )
    {
        ListElt<T> * l_pNew = new ListElt<T>( p_T );

        l_pNew->m_pNext = ( ListElt<T>* )Ptr;
        l_pNew->m_pPrev = ( ( ListElt<T>* )Ptr )->m_pPrev;

        l_pNew->m_pNext->m_pPrev = l_pNew;
        if( l_pNew->m_pPrev )
            l_pNew->m_pPrev->m_pNext = l_pNew;

        m_Count++;

        // New head ?
        if( !l_pNew->m_pPrev )
            m_pHead = l_pNew;

        return l_pNew;
    }
    else
        return nullptr;
}

//==============================================================================
template <class T> void * List<T>::AddBefore( void * Ptr, T && p_T )
{
    if( Ptr )
    {
        ListElt<T> * l_pNew = new ListElt<T>( std::move(p_T) );

        l_pNew->m_pNext = ( ListElt<T>* )Ptr;
        l_pNew->m_pPrev = ( ( ListElt<T>* )Ptr )->m_pPrev;

        l_pNew->m_pNext->m_pPrev = l_pNew;
        if( l_pNew->m_pPrev )
            l_pNew->m_pPrev->m_pNext = l_pNew;

        m_Count++;

        // New head ?
        if( !l_pNew->m_pPrev )
            m_pHead = l_pNew;

        return l_pNew;
    }
    else
        return nullptr;
}

//==============================================================================
template <class T> void * List<T>::AddAfter( void * Ptr, const T & p_T )
{
    if( Ptr )
    {
        ListElt<T> * l_pNew = new ListElt<T>( p_T );

        l_pNew->m_pNext = ( ( ListElt<T>* )Ptr )->m_pNext;
        l_pNew->m_pPrev = ( ListElt<T>* )Ptr;

        l_pNew->m_pPrev->m_pNext = l_pNew;
        if( l_pNew->m_pNext )
            l_pNew->m_pNext->m_pPrev = l_pNew;

        m_Count++;

        // New tail ?
        if( !l_pNew->m_pNext )
            m_pTail = l_pNew;

        return l_pNew;
    }
    else
        return nullptr;
}

//==============================================================================
template <class T> void * List<T>::AddAfter( void * Ptr, T && p_T )
{
    if( Ptr )
    {
        ListElt<T> * l_pNew = new ListElt<T>( std::move(p_T) );

        l_pNew->m_pNext = ( ( ListElt<T>* )Ptr )->m_pNext;
        l_pNew->m_pPrev = ( ListElt<T>* )Ptr;

        l_pNew->m_pPrev->m_pNext = l_pNew;
        if( l_pNew->m_pNext )
            l_pNew->m_pNext->m_pPrev = l_pNew;

        m_Count++;

        // New tail ?
        if( !l_pNew->m_pNext )
            m_pTail = l_pNew;

        return l_pNew;
    }
    else
        return nullptr;
}

//==============================================================================
template <class T> void * List<T>::AddHead( const T & p_T )
{
    if( m_Count )
    {
        ListElt<T> * l_pNew = new ListElt<T>( p_T );

        l_pNew->m_pNext = m_pHead;
        m_pHead->m_pPrev = l_pNew;

        m_Count++;

        m_pHead = l_pNew;
    }
    else
    {
        m_pHead = m_pTail = new ListElt<T>( p_T );

        m_Count = 1;
    }

    return m_pHead;
}

//==============================================================================
template <class T> void * List<T>::AddHead( T && p_T )
{
    if( m_Count )
    {
        ListElt<T> * l_pNew = new ListElt<T>( std::move(p_T) );

        l_pNew->m_pNext = m_pHead;
        m_pHead->m_pPrev = l_pNew;

        m_Count++;

        m_pHead = l_pNew;
    }
    else
    {
        m_pHead = m_pTail = new ListElt<T>( p_T );

        m_Count = 1;
    }

    return m_pHead;
}

//==============================================================================
template <class T> void * List<T>::AddTail( const T & p_T )
{
    ListElt<T> * l_pNew = new ListElt<T>( p_T );

    if( m_Count )
    {
        l_pNew->m_pPrev = m_pTail;
        m_pTail->m_pNext = l_pNew;

        m_Count++;

        m_pTail = l_pNew;
    }
    else
    {
        m_pHead = m_pTail = l_pNew;

        m_Count = 1;
    }

    return m_pTail;
}

//==============================================================================
template <class T> void * List<T>::AddTail( T && p_T )
{
    ListElt<T> * l_pNew = new ListElt<T>( std::move(p_T) );

    if( m_Count )
    {
        l_pNew->m_pPrev = m_pTail;
        m_pTail->m_pNext = l_pNew;

        m_Count++;

        m_pTail = l_pNew;
    }
    else
    {
        m_pHead = m_pTail = l_pNew;

        m_Count = 1;
    }

    return m_pTail;
}

//==============================================================================
template <class T> void * List<T>::Append( const List<T> & pL )
{
    // Check
    if( this == &pL )
        return m_pTail;

    // Append values
    for( void * iPos=pL.HeadPos() ; iPos ; pL.Next( iPos ) )
        AddTail( pL.GetAt( iPos ) );

    return m_pTail;
}

//==============================================================================
template <class T> void * List<T>::Append( List<T> && pL )
{
    // Check
    if( this == &pL )
        return m_pTail;

    // Check
    if( pL.Count() == 0 )
        return m_pTail;

    //
    // Here: pL is non-empty
    //

    // This list is empty?
    if( m_Count == 0 )
    {
        // Yes: steal elements
        m_pHead = pL.m_pHead;
        m_pTail = pL.m_pTail;
        m_Count += pL.m_Count;
    }
    else
    {
        // Non-empty: attach

        // Link
        m_pTail->m_pNext = pL.m_pHead;
        pL.m_pHead->m_pPrev = m_pTail;

        // Move tail
        m_pTail = pL.m_pTail;

        // Update count
        m_Count += pL.m_Count;
    }

    // Nullifiy other
    pL.m_pHead = nullptr;
    pL.m_pTail = nullptr;
    pL.m_Count = 0;

    return m_pTail;
}

//==============================================================================
template <class T> void * List<T>::Prepend( const List<T> & pL )
{
    for( void * iPos=pL.TailPos() ; iPos ; pL.Prev( iPos ) )
        AddHead( pL.GetAt( iPos ) );

    return m_pHead;
}

//==============================================================================
template <class T> void * List<T>::AddIntoIncrList( const T & p_T, bool p_Unique, AddIntoReturnPolicy pPolicy/*=AddIntoReturnPolicy::NullIfNoAdd*/ )
{
    // Find element higher or equal
    void * iPos = m_pHead;
    while( iPos )
    {
        if( GetAt( iPos ) >= p_T )
            break;

        Next( iPos );
    }

    if( !iPos )
    {
        // No higher element found : this one is the highest => tail
        return AddTail( p_T );
    }
    else if( GetAt( iPos ) != p_T )
    {
        // The element found is not equal, hence greater => add before
        return AddBefore( iPos, p_T );
    }
    else
    {
        // The element found is equal to the one considered : duplicate or not
        if( p_Unique )
		{
			if( pPolicy == AddIntoReturnPolicy::NullIfNoAdd )
				return nullptr;
			else
			if( pPolicy == AddIntoReturnPolicy::PosOfIdentical )
				return iPos;
			else
                LzLogException("", "Unknown AddIntoReturnPolicy!")
		}
        else
            return AddBefore( iPos, p_T );
    }
}

//==============================================================================
template <class T> void * List<T>::AddIntoDecrList( const T & p_T, bool p_Unique, AddIntoReturnPolicy pPolicy/*=AddIntoReturnPolicy::NullIfNoAdd*/ )
{
    // Find element smaller or equal
    void * iPos = m_pHead;
    while( iPos )
    {
        if( GetAt( iPos ) <= p_T )
            break;

        Next( iPos );
    }

	if( !iPos )
    {
        // No smaller element found : this one is the smallest => tail
        return AddTail( p_T );
    }
    else
    if( GetAt( iPos ) != p_T )
    {
        // The element found is not equal, hence smaller => add before
        return AddBefore( iPos, p_T );
    }
    else
    {
        // The element found is equal to the one considered : duplicate or not
        if( p_Unique )
		{
			if( pPolicy == AddIntoReturnPolicy::NullIfNoAdd )
				return nullptr;
			else
			if( pPolicy == AddIntoReturnPolicy::PosOfIdentical )
				return iPos;
			else
                LzLogException("", "Unknown AddIntoReturnPolicy!")
		}
        else
            return AddBefore( iPos, p_T );
    }
}

//==============================================================================
template <class T> void List<T>::AddIntoIncrList( const List<T> & p_L, bool p_Unique )
{
    for( void * Pos = p_L.HeadPos() ; Pos ; p_L.Next( Pos ) )
        AddIntoIncrList( p_L.GetAt( Pos ), p_Unique );
}

//==============================================================================
template <class T> void List<T>::AddIntoDecrList( const List<T> & p_L, bool p_Unique )
{
    for( void * Pos = p_L.HeadPos() ; Pos ; p_L.Next( Pos ) )
        AddIntoDecrList( p_L.GetAt( Pos ), p_Unique );
}

//==============================================================================
template <class T> void List<T>::MergeWithIncrList( const List<T> & pL, bool pUniqueFromBtoA )
{
    void * iPosA = m_pHead;
    void * iPosB = pL.m_pHead;

    while( iPosA && iPosB )
    {
        const T & lA = GetAt( iPosA );
        const T & lB = pL.GetAt( iPosB );

        if( lA > lB )
        {
            AddBefore( iPosA, lB );
            pL.Next( iPosB );
        }
        else if( lA == lB )
        {
            // Allow multiple instances?
            if( !pUniqueFromBtoA )
                AddBefore( iPosA, lB );

			// Do not move pointer for A. To preserve uniqueness in A if B has duplicate elements: {1, 3, 5} + {2, 3, 3, 8} --> 1 >> 2 >> 3 >> 5 >> 8
			//Next( iPosA );
			//
			// If you DO move pointer for A: {1, 3, 5} + {2, 3, 3, 8} --> 1 >> 2 >> 3 >> 3 >> 5 >> 8
			//
            pL.Next( iPosB );
        }
        else
        {
            // lA < lB
            Next( iPosA );
        }
    }

    while( iPosB )
    {
        AddTail( pL.GetAt( iPosB ) );
        pL.Next( iPosB );
    }
}

//==============================================================================
template <class T> void List<T>::MergeUniqueIncrLists( const List<T> & pL )
{
    void * iPosA = m_pHead;
    void * iPosB = pL.m_pHead;

    while( iPosA && iPosB )
    {
        const T & lA = GetAt( iPosA );
        const T & lB = pL.GetAt( iPosB );

        if( lA > lB )
        {
			// iPosA points to a position where the smaller B value must be inserted BEFORE
            AddBefore( iPosA, lB );
            pL.Next( iPosB );
        }
        else if( lA == lB )
        {
			// Move both pointers
			Next( iPosA );
            pL.Next( iPosB );
        }
        else
        {
            // lA < lB
            Next( iPosA );
        }
    }

    while( iPosB )
    {
        AddTail( pL.GetAt( iPosB ) );
        pL.Next( iPosB );
    }
}

//==============================================================================
template <class T> void List<T>::SubstractIncrLists( const List<T> & pOther, /*List*/SubstractMode pMode )
{
    // Check both lists are not empty
    if( !Count() || !pOther.Count() )
        return;

    // Obviously non-intersecting
    if( GetTail() < pOther.GetHead() || pOther.GetTail() < GetHead() )
        return;

    // Compute substraction
    void * iPosA = m_pHead;
    void * iPosB = pOther.m_pHead;

    while( iPosA && iPosB )
    {
        const T & lA = GetAt( iPosA );
        const T & lB = pOther.GetAt( iPosB );

        if( lA > lB )
            pOther.Next( iPosB );
        else if( lA == lB )
        {
            DelAtAndNext( iPosA );

            if( pMode == /*List*/SubstractMode::SubstractCount )
                pOther.Next( iPosB );
        }
        else
        {
            // lA < lB
            Next( iPosA );
        }
    }
}

//==============================================================================
template <class T> List<T> List<T>::InterIncrLists( const List<T> & pOther, bool pUniqueInInter ) const
{
    List<T> lInter;

    // One of the two lists is empty
    if( !Count() || !pOther.Count() )
        return lInter;

    // Obviously non-intersecting
    if( GetTail() < pOther.GetHead() || pOther.GetTail() < GetHead() )
        return lInter;

    // Compute intersection
    void * iPosA = m_pHead;
    void * iPosB = pOther.m_pHead;

    while( iPosA && iPosB )
    {
        const T & lA = GetAt( iPosA );
        const T & lB = pOther.GetAt( iPosB );

        if( lA > lB )
            pOther.Next( iPosB );
        else if( lA == lB )
        {
            if( !pUniqueInInter || !lInter.Count() || lInter.GetTail() != lA )
                lInter.AddTail( lA );

            Next( iPosA );
            pOther.Next( iPosB );
        }
        else
        {
            // lA < lB
            Next( iPosA );
        }
    }

    return lInter;
}

//==============================================================================
template<class T> T List<T>::Min() const
{
    // Check
    if( !m_Count )
        LzLogException("", "Cannot take min of empty list!")

    // Find min
    T lMin = T::MaxValue;
    BrowseList( i, *this )
    {
        T lVal = GetAt( i );

        if( lMin > lVal )
            lMin = lVal;
    }

    return lMin;
}

//==============================================================================
template<class T> T List<T>::Max() const
{
    // Check
    if( !m_Count )
        LzLogException("", "Cannot take max of empty list!")

    // Find max
    T lMax = T::MinValue;
    BrowseList( i, *this )
    {
        T lVal = GetAt( i );

        if( lMax < lVal )
            lMax = lVal;
    }

    return lMax;
}

//==============================================================================
template<class T> double List<T>::Sum() const
{
    double lSum = 0;

    // Compute sum
    BrowseList( i, *this )
        lSum += GetAt( i );

	// Return sum
    return lSum;
}

//==============================================================================
template<class T> double List<T>::Mean() const
{
    // Check
    if( !m_Count )
        LzLogException("", "Cannot compute mean of empty list!")

	// Return mean
    return Sum() / (double)m_Count;
}

//==============================================================================
template<class T> void List<T>::Reverse()
{
	// Swap local orientation
	ListElt<T> * iPos = m_pHead;
	while( iPos )
	{
		ListElt<T> * lpLE = iPos;
		ListElt<T> * iNextPos = lpLE->m_pNext;

		// Swap
		ListElt<T> * lTmp = lpLE->m_pNext;
		lpLE->m_pNext = lpLE->m_pPrev;
		lpLE->m_pPrev = lTmp;

		// Move pointer
		iPos = iNextPos;
	}
	
	// Swap global orientation
	ListElt<T> * lTmp = m_pHead;
	m_pHead = m_pTail;
	m_pTail = lTmp;
}

//==============================================================================
// iPos must be a valid position in the list otherwise an infinite loop occurs
template<class T> void List<T>::RotateToHead( void * pNewHeadPos )
{
    // Check
    if( m_Count == 0 )
        LzLogException("", "Cannot rotate empty list!")

    // Check
    if( m_Count == 1 )
    {
        // Double check
        if( m_pHead == pNewHeadPos )
            return;
        else
            LzLogException("", "This singleton list does not contain the provided position!")
    }

    //*** List has at least 2 elements

    // Rotate
    while( m_pHead != pNewHeadPos )
    {
        // Set tail
        m_pTail->m_pNext = m_pHead;
        m_pTail->m_pNext->m_pPrev = m_pTail;

        // Set new head
        m_pHead = m_pTail->m_pNext->m_pNext;
        m_pTail->m_pNext->m_pNext = nullptr;

        // Move tail
        m_pTail = m_pTail->m_pNext;

        // Set head
        m_pHead->m_pPrev = nullptr;
    }
}
#pragma endregion


#pragma region "Comparison"
//==============================================================================
template <class T> bool List<T>::operator==( const List<T> & pL ) const
{
    if( Count() != pL.Count() )
        return false;

    BrowseTwoLists( iA, *this, iB, pL )
    {
        if( GetAt( iA ) != pL.GetAt( iB ) )
            return false;
    }

    return true;
}

//==============================================================================
template <class T> bool List<T>::operator!=( const List<T> & pL ) const
{
    if( Count() != pL.Count() )
        return true;

    BrowseTwoLists( iA, *this, iB, pL )
    {
        if( GetAt( iA ) != pL.GetAt( iB ) )
            return true;
    }

    return false;
}

//==============================================================================
template <class T> bool List<T>::IsDirectShiftOf( const List<T> & pL ) const
{
    // Check
    if( Count() != pL.Count() )
        return false;

    // Empty lists are similar
    if( !Count() )
        return true;

    // For every possible starting point
    for( void * iStartPos = HeadPos() ; iStartPos ; Next( iStartPos ) )
    {
        void * iPosA = iStartPos;
        for( void * iPosB = pL.HeadPos() ; iPosB ; pL.Next( iPosB ) ) // Direct order
        {
            // Test list elements
            if( GetAt( iPosA ) != pL.GetAt( iPosB ) )
                goto Next;

            // Move my position
            Next( iPosA );
            if( !iPosA )
                iPosA = HeadPos();
        }

        // The lists are identical for this shift
        return true;

Next: /* test next shift... */
        ;
    }

    // The lists are not similar
    return false;
}

//==============================================================================
template <class T> bool List<T>::IsReverseShiftOf( const List<T> & pL ) const
{
    // Check
    if( Count() != pL.Count() )
        return false;

    // Empty lists are similar
    if( !Count() )
        return true;

    // For every possible starting point
    for( void * iStartPos = HeadPos() ; iStartPos ; Next( iStartPos ) )
    {
        void * iPosA = iStartPos;
        for( void * iPosB = pL.TailPos() ; iPosB ; pL.Prev( iPosB ) ) // Reverse order
        {
            // Test list elements
            if( GetAt( iPosA ) != pL.GetAt( iPosB ) )
                goto Next;

            // Move my position
            Next( iPosA );
            if( !iPosA )
                iPosA = HeadPos();
        }

        // The lists are identical for this shift
        return true;

Next: /* test next shift... */
        ;
    }

    // The lists are not similar
    return false;
}

//==============================================================================
template <class T> bool List<T>::Has( const T & pT ) const
{
    for( void * i=HeadPos() ; i ; Next(i) )
	{
		if( GetAt(i) == pT )
			return true;
	}

	return false;
}
#pragma endregion


#pragma region "Conversion to/from container"
//==============================================================================
template <class T> void List<T>::ToArray( T * pArray ) const
{
    int i = 0;

    for( void * iT = m_pHead ; iT ; Next( iT ) )
        pArray[i++] = GetAt( iT );
}

//==============================================================================
template <class T> void List<T>::FromArray( const T * pArray, size_t pCount )
{
    DelAll();

    for( size_t i=0 ; i<pCount ; i++ )
        AddTail( pArray[i] );
}

//==============================================================================
template <class T> void List<T>::FromArrayReverse( const T * pArray, size_t pCount )
{
    DelAll();

    for( size_t i=0 ; i<pCount ; i++ )
        AddHead( pArray[i] );
}

//==============================================================================
template <class T> void List<T>::FromVector( const vector<T> & pVector )
{
    DelAll();

    for( size_t i=0 ; i<pVector.size() ; i++ )
        AddTail( pVector[i] );
}

//==============================================================================
template <class T> void List<T>::FromVector( const Vector<T> & pVector )
{
    DelAll();

    for( size_t i=0 ; i<pVector.Size() ; i++ )
        AddTail( pVector[i] );
}

//==============================================================================
template <class T> void List<T>::ToVector( vector<T> & pVector ) const
{
    pVector.resize( Count() );
    size_t i = 0;
    BrowseList( iT, *this )
        pVector[i++] = GetAt( iT );
}

//==============================================================================
template <class T> void List<T>::ToVector( Vector<T> & pVector ) const
{
    pVector.Resize( Count() );
    size_t i = 0;
    BrowseList( iT, *this )
        pVector[i++] = GetAt( iT );
}
#pragma endregion


#pragma region "Deletion"
//==============================================================================
template <class T> void List<T>::DelAll()
{
    ListElt<T> * l_pTarget = m_pHead;
    while( l_pTarget )
    {
        ListElt<T> * l_pNext = l_pTarget->m_pNext;
        delete l_pTarget;
        l_pTarget = l_pNext;
    }

    m_pHead = m_pTail = nullptr;
    m_Count = 0;
}

//==============================================================================
template <class T> void List<T>::DelAt( void * Ptr )
{
    if( Ptr )
    {
        ListElt<T> * l_pTarget = ( ListElt<T>* )Ptr;

        if( l_pTarget->m_pPrev )
            l_pTarget->m_pPrev->m_pNext = l_pTarget->m_pNext;
        else
            m_pHead = l_pTarget->m_pNext;

        if( l_pTarget->m_pNext )
            l_pTarget->m_pNext->m_pPrev = l_pTarget->m_pPrev;
        else
            m_pTail = l_pTarget->m_pPrev;

        // Destroy target
        delete l_pTarget;
        m_Count--;
    }
}

//==============================================================================
template <class T> void List<T>::DelFromToTail( void * pPtr )
{
    while( pPtr )
        DelAtAndNext( pPtr );
}

//==============================================================================
template <class T> void List<T>::DelElts( const T & p_T )
{
    ListElt<T> * iP = m_pHead;
    while( iP )
    {
        if( iP->m_T == p_T )
        {
            // Del at and next
            ListElt<T> * lpNext = iP->m_pNext;
            ListElt<T> * lpPrev = iP->m_pPrev;

            if( lpPrev )
                lpPrev->m_pNext = lpNext;
            else
                m_pHead = lpNext;

            if( lpNext )
                lpNext->m_pPrev = lpPrev;
            else
                m_pTail = lpPrev;

            delete iP;

            m_Count--;

            iP = lpNext;
        }
        else
            iP = iP->m_pNext;
    }
}

//==============================================================================
template <class T> void List<T>::DelAtAndNext( void * & Ptr )
{
    if( Ptr )
    {
        ListElt<T> * l_pTarget = ( ListElt<T>* )Ptr;

        ListElt<T> * lpNext = l_pTarget->m_pNext;
        ListElt<T> * lpPrev = l_pTarget->m_pPrev;

        if( lpPrev )
            lpPrev->m_pNext = lpNext;
        else
            m_pHead = lpNext;

        if( lpNext )
            lpNext->m_pPrev = lpPrev;
        else
            m_pTail = lpPrev;

        delete l_pTarget;

        m_Count--;

        Ptr = lpNext;
    }
    else
        Ptr = nullptr;
}

//==============================================================================
template <class T> void List<T>::FindAndDel( const T & pT, bool pExcIfNotFound )
{
    void * lpPos = FindElt( pT );

    if( lpPos )
        DelAt( lpPos );
    else
    if( pExcIfNotFound )
        LzLogException("", "Cannot find and delete element! Element not found in list (count= "<<Count()<<")")
}
#pragma endregion


#pragma region "Debug"
//==============================================================================
template <class T> std::string List<T>::ListToString() const
{
	std::stringstream lStrStr;
	lStrStr << "H >> ";
	BrowseList( i, *this )
		lStrStr << GetAt( i ) << " >> ";
	lStrStr << "T";

	return lStrStr.str();
}

//==============================================================================
template <class T> std::ostream & operator<<( std::ostream & pOut, const List<T> & pList )
{
	pOut << pList.ListToString();
	return pOut;
}

//==============================================================================
template <class T> std::string List<T>::ListOfListToString() const
{
	std::stringstream lStrStr;
	lStrStr << "H >> ";
	BrowseList( i, *this )
		lStrStr << "{ " << GetAt( i ).ListToString() << " } >> ";
	lStrStr << "T";

	return lStrStr.str();
}
#pragma endregion
}
