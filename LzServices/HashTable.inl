#include <typeinfo>

#include <LzServices/LzLog.h>
//#include <_LzServices/LzLogSelector.h>


namespace LzServices
{
#pragma region "Construction & destruction"
//================================================================================
template <class T> HashTable<T>::HashTable( unsigned int pSize/*=1000*/ )
    : mSize( pSize )
{
    mpLists = new List<T>[mSize];
}

//================================================================================
template <class T> HashTable<T>::~HashTable()
{
    Free();
}

//================================================================================
template <class T> void HashTable<T>::SetSize( unsigned int pSize )
{
    Free();

    mSize = pSize;
    mpLists = new List<T>[mSize];
}

//================================================================================
template <class T> void HashTable<T>::Free()
{
    delete [] mpLists;
    mpLists = nullptr;
    mSize = 0;
}
#pragma endregion


#pragma region "Key generators"
//================================================================================
template <class T> unsigned int HashTable<T>::HashKey( int pFromInteger )
{
    unsigned int lCode;

    // Select floating point container
    if( sizeof(unsigned int) == sizeof(float) )
    {
        float lCos = static_cast<float>( std::cos( pFromInteger ) );
        lCode = *( reinterpret_cast<unsigned int *>( &lCos ) );
    }
    else
/*  if( sizeof(unsigned int) == sizeof(double) )
    {
        double lCos = std::cos( pFromInteger );
        lpCode = reinterpret_cast<unsigned int *>( &lCos );
    }
    else
*/
        LzLogException("", "Cannot use default hash function due to format mismatch! sizeof(double)= "<<sizeof(double)<<", sizeof(float)= "<<sizeof(float)<<", sizeof(unsigned int)="<<sizeof(unsigned int)<<".")

    return lCode;

// Crappy hash functions ==> return max_uint very quickly
//    return static_cast<unsigned int>( pFromInteger * std::cos( pFromInteger ) );
//    return static_cast<unsigned int>( pFromInteger * std::cos( 3.1415926 * pFromInteger / 180.0 ) );
}

//================================================================================
template <class T> unsigned int HashTable<T>::HashKey( const std::string & pFromString )
{
    // Sum all char values in string
    int lTotal = 0;
    for( unsigned int i = 0 ; i < pFromString.length() ; i++ )
        lTotal += ( int )pFromString[i];

    return HashKey( lTotal );
}
#pragma endregion


#pragma region "Usage"
//================================================================================
template <class T> void HashTable<T>::FlushToVector( vector<T> & pToVector ) const
{
    pToVector.clear();

    for( unsigned int n=0 ; n<mSize ; n++ )
    {
        List<T> & lList = mpLists[n];
        BrowseList( iT, lList )
            pToVector.push_back( lList.GetAt( iT ) );
    }
}

//================================================================================
template <class T> void HashTable<T>::DelAll()
{
    for( unsigned int n = 0 ; n < mSize ; n++ )
        mpLists[n].DelAll();
}
#pragma endregion


#pragma region "Browsing"
//================================================================================
template <class T> HashPos HashTable<T>::HeadPos() const
{
    // Look for first position in table
    HashPos lHPos;
    for( unsigned int i = 0 ; i < mSize ; i++ )
    {
        if( mpLists[i].Count() )
        {
            lHPos.mIdx = i;
            lHPos.mPos = mpLists[i].HeadPos();

            return lHPos;
        }
    }

    // Tabula rasa: returne empty position
    lHPos.mIdx = 0;
    lHPos.mPos = nullptr;
    return lHPos;
}

//================================================================================
template <class T> void HashTable<T>::Next( HashPos & pPos ) const
{
    if( pPos.mPos )
    {
        // Next position in current list is ok?
        mpLists[pPos.mIdx].Next( pPos.mPos );
        if( pPos.mPos )
            return;

        // Nope...

        // Look for next position in table
        for( unsigned int i = pPos.mIdx + 1 ; i < mSize ; i++ )
        {
            if( mpLists[i].Count() )
            {
                pPos.mIdx = i;
                pPos.mPos = mpLists[i].HeadPos();

                return ;
            }
        }

        // Still nope... return with pPos.mPos==nullptr
    }
}

//================================================================================
template <class T> T & HashTable<T>::GetAt( const HashPos & pPos )
{
    return mpLists[ pPos.mIdx ].GetAt( pPos.mPos );
}

//================================================================================
template <class T> const T & HashTable<T>::GetAt( const HashPos & pPos ) const
{
    return mpLists[ pPos.mIdx ].GetAt( pPos.mPos );
}
#pragma endregion


#pragma region "Debug"
//================================================================================
template <class T> void HashTable<T>::Log( const std::string & pTag/*="HashTable"*/ ) const
{
    int lEmptyLists = 0;
    int lMaxCount = 0;

    // Statistics on hash table
    {
//        LzLogN("", "List counts")

        for( unsigned int n=0 ; n<mSize ; n++ )
        {
            int lCount = mpLists[n].Count();

            if( !lCount )
                lEmptyLists++;

            if( lMaxCount < lCount )
                lMaxCount = lCount;


//            LzLogM("", "List "<<n<<" = "<<lCount)
        }
    }

	// Log
    LzLogMsg( "", "Statistics for '"<<pTag<<"<"<<typeid(T).name()<<","<< mSize<<">': "<<lEmptyLists<<" empty list(s), max list count="<<lMaxCount<<"." )
}
#pragma endregion
}
