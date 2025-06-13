#include "LzLog.h"
#include "List.h"
#ifdef USING_QT
    #include <QCoreApplication>
    #include <QDir>
#endif
//
#include <chrono>
#include <time.h>
//
#include <string>
#include <iostream>
#include <iomanip> // std::precision
#include <algorithm> // std::transform
#include <atomic> // std::atomic_bool
#include <sys/stat.h> // FileExists
//#include <cctype>     // testing caracters, toupper, ... include needed?
//#if defined(OS_LINUX)
//	#include <stdio.h>
//#endif
#ifdef USING_Cpp_17
    #include <filesystem>
#endif
#include <mutex>


namespace LzServices
{
using std::cout;
using std::endl;
using std::flush;







#pragma region "Data subsampling"
//================================================================================
size_t GetSubSampSmallStep( size_t pOriginalSize, size_t pSampleSize )
{
    // Check
    if( pOriginalSize==0 || pSampleSize==0 )
        LzLogException("", "Null original or sample size!")

    // Check
    if( pSampleSize >= pOriginalSize )
        return 1;

    // Return smaller step >= 1, since here pOriginalSize > pSampleSize
    return std::floor( pOriginalSize / static_cast<float>(pSampleSize) );
}

//================================================================================
size_t GetActualSamplesCount( size_t pOriginalSize, size_t pSampleStep )
{
    // Check
    if( pOriginalSize==0 || pSampleStep==0 )
        LzLogException("", "Null original size or sample step!")

    // Get largest sample index: M, such as M*Step < Size, i.e. M*Step <= Size - 1
    // Using floor(x) <= x
    //
    const size_t M = std::floor( (pOriginalSize - 1) / static_cast<float>(pSampleStep) );

    // Number of samples is M+1 since indices start at 0
    return M + 1;
}
/* TESTING DIFFERENCE BETWEEN APPROX SAMPLES COUNT AND ACTUAL SAMPLES COUNT
 *
    for( size_t vs=1999 ; vs<=50999 ; vs+=500 )
    {
        const size_t lStep = LzServices::GetSubSampSmallStep( vs, 1000 );
        const size_t lActual = LzServices::GetActualSamplesCount( vs, lStep );
        const size_t lLastIdx = (lActual-1) * lStep;

        LzLogM("", "VS= "<<vs<<"\tStep= "<<lStep<<"\tActual samples= "<<lActual<<"\tLast index= "<<lLastIdx
                   <<"\tUnsampled elements= "<<(vs-1-lLastIdx))
    }
 *
 * OUTPUT
 *
 * VS= 1999    Step= 1    Actual samples= 1999    Last index= 1998     Unsampled elements= 0
 * VS= 2499    Step= 2    Actual samples= 1250    Last index= 2498     Unsampled elements= 0
 * VS= 2999    Step= 2    Actual samples= 1500    Last index= 2998     Unsampled elements= 0
 * VS= 3499    Step= 3    Actual samples= 1167    Last index= 3498     Unsampled elements= 0
 * VS= 3999    Step= 3    Actual samples= 1333    Last index= 3996     Unsampled elements= 2
 * VS= 4499    Step= 4    Actual samples= 1125    Last index= 4496     Unsampled elements= 2
 * VS= 4999    Step= 4    Actual samples= 1250    Last index= 4996     Unsampled elements= 2
 * VS= 5499    Step= 5    Actual samples= 1100    Last index= 5495     Unsampled elements= 3
 * VS= 5999    Step= 5    Actual samples= 1200    Last index= 5995     Unsampled elements= 3
 * VS= 6499    Step= 6    Actual samples= 1084    Last index= 6498     Unsampled elements= 0
 * VS= 6999    Step= 6    Actual samples= 1167    Last index= 6996     Unsampled elements= 2
 * VS= 7499    Step= 7    Actual samples= 1072    Last index= 7497     Unsampled elements= 1
 * VS= 7999    Step= 7    Actual samples= 1143    Last index= 7994     Unsampled elements= 4
 * VS= 8499    Step= 8    Actual samples= 1063    Last index= 8496     Unsampled elements= 2
 * VS= 8999    Step= 8    Actual samples= 1125    Last index= 8992     Unsampled elements= 6
 * VS= 9499    Step= 9    Actual samples= 1056    Last index= 9495     Unsampled elements= 3
 * VS= 9999    Step= 9    Actual samples= 1111    Last index= 9990     Unsampled elements= 8
 * VS= 10499   Step= 10   Actual samples= 1050    Last index= 10490    Unsampled elements= 8
 * VS= 10999   Step= 10   Actual samples= 1100    Last index= 10990    Unsampled elements= 8
 * VS= 11499   Step= 11   Actual samples= 1046    Last index= 11495    Unsampled elements= 3
 * VS= 11999   Step= 11   Actual samples= 1091    Last index= 11990    Unsampled elements= 8
 * VS= 12499   Step= 12   Actual samples= 1042    Last index= 12492    Unsampled elements= 6
 * VS= 12999   Step= 12   Actual samples= 1084    Last index= 12996    Unsampled elements= 2
 * VS= 13499   Step= 13   Actual samples= 1039    Last index= 13494    Unsampled elements= 4
 * VS= 13999   Step= 13   Actual samples= 1077    Last index= 13988    Unsampled elements= 10
 * VS= 14499   Step= 14   Actual samples= 1036    Last index= 14490    Unsampled elements= 8
 * VS= 14999   Step= 14   Actual samples= 1072    Last index= 14994    Unsampled elements= 4
 * VS= 15499   Step= 15   Actual samples= 1034    Last index= 15495    Unsampled elements= 3
 * VS= 15999   Step= 15   Actual samples= 1067    Last index= 15990    Unsampled elements= 8
 * VS= 16499   Step= 16   Actual samples= 1032    Last index= 16496    Unsampled elements= 2
 * VS= 16999   Step= 16   Actual samples= 1063    Last index= 16992    Unsampled elements= 6
 * VS= 17499   Step= 17   Actual samples= 1030    Last index= 17493    Unsampled elements= 5
 * VS= 17999   Step= 17   Actual samples= 1059    Last index= 17986    Unsampled elements= 12
 * VS= 18499   Step= 18   Actual samples= 1028    Last index= 18486    Unsampled elements= 12
 * VS= 18999   Step= 18   Actual samples= 1056    Last index= 18990    Unsampled elements= 8
 * VS= 19499   Step= 19   Actual samples= 1027    Last index= 19494    Unsampled elements= 4
 * VS= 19999   Step= 19   Actual samples= 1053    Last index= 19988    Unsampled elements= 10
 * VS= 20499   Step= 20   Actual samples= 1025    Last index= 20480    Unsampled elements= 18
 * [...]
 *
 */

//================================================================================
/*
 * The sampling is done in two batches. The first one with a sampling step that
 * ensures that all entries are covered. The second batch distributes any remaining
 * samples among the intervals between the samples of the first batch.
 *
 * VS : vector size = number of entries to sample from
 *  N : number of samples
 *
 *
 * Example 1: VS=19, N=7
 * 19 = 2*7+5 = QN+R
 * Sampling period = Q+1 = 2+1 = 3
 *                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
 * First batch (7 elements)   X        X        X        X        X        X        X
 * Second batch (0 elements)
 *
 *
 * Example 2: VS=15, N=7
 * 15 = 2*7+1 = QN+R
 * Sampling period = Q+1 = 2+1 = 3
 *                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
 * First batch (5 elements)   X        X        X        X        X
 * Second batch (2 elements)        X                                   X
 *
 *
 * Example 3: VS=23, N=11
 * 23 = 2*11+1 = QN+R
 * Sampling period = Q+1 = 2+1 = 3
 *                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
 * First batch (8 elements)   X        X        X        X        X        X        X        X
 * Second batch (3 elements)        X                          X                          X
 *
 */
void GetSubSampIdx( size_t pNbEntries, size_t pNbSamples, vector<size_t> & pTo )
{
    // Check
    if( pNbEntries == 0 )
        LzLogException("", "Cannot subsample from empty source!")

    // Check
    if( pNbSamples == 0 )
        LzLogException("", "Cannot subsample to empty set!")

    // Check
    if( pNbSamples > pNbEntries )
        LzLogException("", "More samples than entries! Impossible subsampling.")

    // Check
    if( pNbSamples == pNbEntries )
        LzLogException("", "As manyu samples as entries! Pointless subsampling.")

    // General case: more data than samples
    //
    // N: number of samples
    // VS: vector size, i.e. number of entries
    //
    // N>=1, VS>N -> VS>=2
    const size_t VS = pNbEntries;
    const size_t  N = pNbSamples;
    //
    pTo.resize( N );

    if( VS % N == 0 )
    {
        // N | VS ==> VS = QN
        // ==> Step = Q
        const size_t Q = VS / N;
        for( size_t i=0 ; i<N ; i++ )
            pTo[ i ] = i*Q;
    }
    else
    {
        // not( N | VS ) ==> VS = QN + R; Q>=1, R>=1, N>1 (N cannot be 1 since 1 | VS)
        // ==> Step = Q+1 >= 2
        const size_t Q = VS / N;
#if 0
        const size_t R = VS % N;
#endif
        //
        const size_t M = (VS - 1) / (Q + 1); // Integer division, equivalent to floor( ... )
        // M: last index in first batch of samples, covering the full span of the source vector
        // floor( X ) <= X
        // Thus, M(Q + 1) <= VS - 1; thus M(Q + 1) < VS
        // M: largest int such as M(Q + 1) < VS, thus (M + 1)(Q + 1) >= VS
#if 0
        LzLogN("", VS<<" = "<<Q<<" * "<<N<<" + "<<R)
        LzLogM("", "Step= Q+1= "<<(Q+1))
        LzLogM("", "M= "<<M)
        LzLogM("", "First: "<<(M+1))
        LzLogM("", "Last:  "<<(N-M-1))
#endif
        // Sample first batch
        for( size_t m=0 ; m<=M ; m++ )
           pTo[m] = m*(Q + 1);

        // No second batch?
        if( M == N-1 ) // i.e. N-M-1 == 0
            return;

        // Here, N-1>M, i.e. N-M-1>0, i.e. N-M-1>=1

        // Offset for second batch, at half sample step
        const size_t F = std::ceil( (Q + 1) / 2.0 ); // F>=1 since Q+1>=2 (since Q>=1)
        // Other admissible value: const int F = std::floor( (Q + 1) / 2.0 ) + 1;

        // Constraint on max value of F
        // F + (Q + 1)*(N - M - 2) < VS; last accessed index must be within the range of entries
        //
        // F + QN - QM - 2Q + N - M - 2 < QN + R
        // F      - QM - 2Q + N - M - 2 <      R
        // F < R + QM + 2Q + M + 2 - N
        // Does this leave room for at least F = 1 ?
        //
        // F < R + M(Q + 1) + (Q + 1) + (Q + 1) - N
        // F < R + (M + 1)(Q + 1) + Q - N
        // F < R + Q + (M + 1)(Q + 1) - N
        //
        // We know that (M + 1)(Q + 1) >= VS and that VS > N, hence (M + 1)(Q + 1) > N
        // and (M + 1)(Q + 1) - N > 0, i.e. (M + 1)(Q + 1) - N >= 1
        // If F is such as: F < R + Q + 1, then it will comply with the constraint above
        // since R + Q + 1 <= R + Q + (M + 1)(Q + 1) - N
        // Furthermore R>=1 and Q>=1, thus the minimal constraint for F is
        // F < 3. Conclusion: F = 1 is an admissible value.
        //
        // Finally, F = ceil( (Q + 1) / 2.0 ) is also an admissible value.
        // A sufficient condition for admissibility of F is: F < 1 + (Q+1) = Q+2, since R>=1
        // And floor( (Q+1)/2 ) <= (Q+1)/2, and (Q+1)/2 < Q+2, since Q+1 < 2Q+4, since 0<Q+3
        //
        // Same proof to show that F = ceil( (Q + 1) / 2.0 ) + 1 is ALSO an admissible value.

        // Spread samples across entries if more than one sample in second batch
        if( N-M-1 >= 2 ) // i.e. N - M - 2 >= 1 i.e. N - M - 2 != 0 so the division below makes sense
        {
            // Boost jump length so that the whole spectrum of entry values is covered again
            // The maximal BOOST value is such that
            // F + BOOST*(Q + 1)*(N - M - 2) < VS, i.e. ... <= VS - 1
            const int BOOST = std::floor( (VS - 1 - F) / static_cast<float>( (Q+1)*(N - M - 2) ) );
            // We know from the above that BOOST = 1 is admissible, thus BOOST >= 1

            // Sample second batch
            for( size_t m=0 ; m<N-M-1 ; m++ )
               pTo[M + 1 + m] = F + m*BOOST*(Q + 1);
        }
        // N-M-1 == 1 ?
        else
        {
            // One single entry in second batch
            pTo[N - 1] = F;
        }
    }
}

//================================================================================
void TEST_GetSubSampIdx()
{
    LzLogN("", "TEST_GetSubSampIdx")

    // Test substrate and param
    vector<bool> lSampled;
    const size_t MAX = 2000;

    // Run test
#if 1
    for( int VS=1 ; VS<=MAX ; VS++ ) // avoid VS==0
    for( int  N=1 ; N < VS  ;  N++ ) // avoid N==0 or N==VS
#else
    int VS = 23;
    int N = 11;
    for( int i=0 ; i<1 ; i++ )
#endif
    {
        // Init entries
        lSampled.assign( VS, false );

        // Get subsampled indices
        vector<size_t> lSubIdx;
        GetSubSampIdx( VS, N, lSubIdx );

        // Check
        for( size_t i=0 ; i<lSubIdx.size() ; i++ )
        {
            // Get subsampled index
            const size_t lIdx = lSubIdx[ i ];

            // Check range
            if( lIdx >= VS )
                LzLogException("", "Index "<<lIdx<<" is out of range!")

            // Check unicity
            if( lSampled[ lIdx ] )
                LzLogException("", "Index "<<lIdx<<" was already sampled!")

            // Mark as sampled
            lSampled[ lIdx ] = true;
        }
#if 0
        LzLogM("", "Sampled indices")
        for( size_t i=0 ; i<lSubIdx.size() ; i++ )
            LzLogM("", lSubIdx[i])
#endif
    }

    LzLogM("", "Done. Max nb. of entries= "<<MAX)
}

//================================================================================
void TEST_GetSubSampIdx_2()
{
    LzLogN("", "TEST_GetSubSampIdx_2")

    // Use detonator to disable and enable log
    Detonator lDet( [] { LzServices::DisableLog(); },
                    [] { LzServices::EnableLog(); } );

const int MAX = 1000;
#if 1
    for( int VS=1 ; VS<=MAX ; VS++ )
    for( int  N=1 ; N < VS  ;  N++ )
#else
    int VS = 19;
    int N = 7;
    for( int i=0 ; i<1 ; i++ )
#endif
    {
        // N < VS

        //**** Ajouter check de unicite des prelevements

        if( VS % N == 0 )
        {
            // N | VS ==> VS = QN
            // ==> Step = Q

        }
        else
        {
            // N /|/ VS ==> VS = QN + R; Q>=1, R>=1, N>1
            // ==> Step = Q+1 >= 2

            const int Q = VS / N;
            const int R = VS % N;

            // floor( X ) <= X
            // Thus, M(Q + 1) <= VS - 1; thus M(Q + 1) < VS
            int M = std::floor( (VS - 1) / (double)(Q + 1) );
            int M2 = (VS - 1) / (Q + 1);

            LzLogN("", VS<<" = "<<Q<<" * "<<N<<" + "<<R)
            LzLogM("", "Step= Q+1= "<<(Q+1))
            LzLogM("", "M= "<<M)
            LzLogM("", "First: "<<(M+1))
            LzLogM("", "Last:  "<<(N-M-1))

            // Check
            if( M != M2 )
                LzLogException("", "Error -2")

            // Check
            if( M*(Q + 1) >= VS )
                LzLogException("", "Error -1")
            //
            // Check
            if( (M + 1)*(Q + 1) < VS )
                LzLogException("", "Error 0")

            // No second part?
            if( M+1 == N )
                continue;

            // Here, N-M-1 >= 1

            const int MAX_F = R + (Q+1)*(M+1) + (Q+1) - N; // + 2;

            // Check
            {
                const int lF_OK = (MAX_F - 1) + (Q+1)*(N - M - 2);
                if( lF_OK >= VS )
                    LzLogException("", "Error 1")

                const int F_KO = (MAX_F + 0) + (Q+1)*(N - M - 2);
                if( F_KO < VS )
                    LzLogException("", "Error 2")

                LzLogM("", lF_OK<<" should be "<<(VS-1))
            }

            const int F = std::ceil( (Q+1) / 2.0 ); // F>=1 since Q+1>=2
            // Other admissible value: const int F = std::floor( (Q + 1) / 2.0 ) + 1;
            LzLogM("", "F= "<<F)

            // Check
            {
                const int LAST_IDX = F + (Q+1)*(N - M - 2);

                if( LAST_IDX >= VS )
                    LzLogException("", "Error 3")

                LzLogM("", "Last index= "<<LAST_IDX<<", vs "<<(VS-1))
            }

            // More than 2 samples in second batch?
            if( N - M - 1 >= 2 ) // i.e. N - M - 2 >= 1 i.e. N - M - 2 != 0
            {
                // Alpha jump
                const int ALPHA = std::floor( (VS - 1 - F) / (double)( (Q+1)*(N - M - 2) ) );
                const int LAST_ALPHA_IDX = F + ALPHA*(Q+1)*(N - M - 2);

                LzLogM("", "ALPHA= "<<ALPHA)
                LzLogM("", "Last ALPHA index= "<<LAST_ALPHA_IDX)

                // Check
                {
                    if( ALPHA < 1 )
                        LzLogException("", "Error 4")

                    if( LAST_ALPHA_IDX >= VS )
                        LzLogException("", "Error 5")
                }
            }
        }
    }
}
#pragma endregion


#pragma region "Read Write using ifstream"
//================================================================================
double readDouble64( ifstream & pFrom, bool pSwap/*=false*/ )
{
    /* http://www.cplusplus.com/reference/ios/skipws/
     * For standard streams, the skipws flag is set on initialization.
     * This flag can be unset with the noskipws manipulator, forcing extraction operations to consider
     * leading whitepaces as part of the content to be extracted.
     */
    unsigned char lDouble[8];
    for( int byte=0 ; byte<8 ; byte++ )
    {
        if( pSwap )
            pFrom >> std::noskipws >> lDouble[7 - byte];
        else
            pFrom >> std::noskipws >> lDouble[byte];
    }

    return *( (double*)lDouble );
}

//================================================================================
float readFloat32( ifstream & pFrom, bool pSwap/*=false*/ )
{
    /* http://www.cplusplus.com/reference/ios/skipws/
     * For standard streams, the skipws flag is set on initialization.
     * This flag can be unset with the noskipws manipulator, forcing extraction operations to consider
     * leading whitepaces as part of the content to be extracted.
     */
    unsigned char lFloat[4];
    for( int byte=0 ; byte<4 ; byte++ )
    {
        if( pSwap )
            pFrom >> std::noskipws >> lFloat[3 - byte];
        else
            pFrom >> std::noskipws >> lFloat[byte];
    }

    return *( (float*)lFloat );
}

//================================================================================
unsigned int readUnsignedInt32( ifstream & pFrom, bool pSwap/*=false*/ )
{
    /* http://www.cplusplus.com/reference/ios/skipws/
     * For standard streams, the skipws flag is set on initialization.
     * This flag can be unset with the noskipws manipulator, forcing extraction operations to consider
     * leading whitepaces as part of the content to be extracted.
     */
    unsigned char lUInt[4];
    for( int byte=0 ; byte<4 ; byte++ )
    {
        if( pSwap )
            pFrom >> std::noskipws >> lUInt[3 - byte];
        else
            pFrom >> std::noskipws >> lUInt[byte];
    }

    return *reinterpret_cast<unsigned int*>( lUInt );
}

//================================================================================
unsigned char * readBinFile( const string & pFileName, unsigned int & pSize )
{
    // Open & check
    ifstream lFile( pFileName, ios::in | ios::binary | ios::ate );
    if( !lFile.is_open() )
        LzLogException("", "Could not read binary file from '"<<pFileName<<"'!")

    // Read size
    lFile.seekg(0, ios::end);
    pSize = lFile.tellg();

    // Read content
    char * lpChars = new char [pSize];
    lFile.seekg (0, ios::beg);
    lFile.read( lpChars, pSize );
    lFile.close();

    // Return as uchars
    return (unsigned char *)lpChars;
}

//================================================================================
void ReadFileToString( const string & pFile, string & pOut )
{
    // Open file
    std::ifstream lFile;
    lFile.open( pFile );

    // Check
    if( !lFile.is_open() )
        LzLogException("", "Impossible to open file: "<<pFile<<"!")

    lFile.seekg(0, std::ios::end);
    pOut.reserve( (unsigned int)lFile.tellg());
    lFile.seekg(0, std::ios::beg);

    pOut.assign((std::istreambuf_iterator<char>(lFile)), std::istreambuf_iterator<char>());
}

//================================================================================
double readDouble64( const unsigned char * pFrom, unsigned int & pIdx, bool pSwap/*=false*/ )
{
    unsigned char lDouble[8];
    for( int byte=0 ; byte<8 ; byte++ )
    {
        if( pSwap )
            lDouble[7 - byte] = pFrom[ pIdx++ ];
        else
            lDouble[byte] = pFrom[ pIdx++ ];
    }

    return *( (double*)lDouble );
}

//================================================================================
float readFloat32( const unsigned char * pFrom, unsigned int & pIdx, bool pSwap/*=false*/ )
{
    unsigned char lFloat[4];
    for( int byte=0 ; byte<4 ; byte++ )
    {
        if( pSwap )
            lFloat[3 - byte] = pFrom[ pIdx++ ];
        else
            lFloat[byte] = pFrom[ pIdx++ ];
    }

    return *( (float*)lFloat );
}

//================================================================================
unsigned int readUnsignedInt32( const unsigned char * pFrom, unsigned int & pIdx, bool pSwap/*=false*/ )
{
    unsigned char lUInt[4];
    for( int byte=0 ; byte<4 ; byte++ )
    {
        if( pSwap )
            lUInt[3 - byte] = pFrom[ pIdx++ ];
        else
            lUInt[byte] = pFrom[ pIdx++ ];
    }

    return *( (unsigned int*)lUInt );
}

#if 0
//================================================================================
// The result of the read is placed in here
vector< vector<double> > lData;

// Read one line at a time into the variable line:
string lLine;
while( std::getline(lFile, lLine) )
{
    vector<double> lLineData;
    std::stringstream lineStream(lLine);

    // Read one value at a time from the line
    double lValue;
    while(lineStream >> lValue)
    {
        // Add the integers from a line to a 1D array (vector)
        lLineData.push_back(lValue);
    }
    // When all the integers have been read, add the 1D array
    // into a 2D array (as one line in the 2D array)
    lData.push_back(lLineData);
}
#endif
#pragma endregion


#pragma region "Logical cores"
unsigned int CountLogicalCores()
{
    return std::thread::hardware_concurrency();
}
#pragma endregion


//================================================================================
//static recursive_mutex sLogMex;
static recursive_mutex & GetLogMex()
{
    static recursive_mutex * spLogMex = nullptr;
    if( !spLogMex )
        spLogMex = new recursive_mutex;

    return *spLogMex;
}
#define LOCK_API std::unique_lock<recursive_mutex> lLock(GetLogMex());


#pragma region "Entropy"
//================================================================================
thread_local std::mt19937 sRandomEngine;
std::mt19937 & RandomEngine()
{
    return sRandomEngine;
}

//================================================================================
thread_local class RandomEnginInit
{
public:
	RandomEnginInit()
	{
        // Init random engine
        auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        RandomEngine().seed( static_cast<unsigned int>(seed) );
    }
}
sRandomEnginInit;
#pragma endregion


#pragma region "Dir / file operations"
//================================================================================
string StartUpPath( const string & pPath/*=""*/ )
{
#if defined(USING_QT)
    // Check
    if( pPath.length() )
        LzLogException("", "LzServices::StartUpPath("<<pPath<<") should be called without any argument! Using QT to determine start-up path.")

    return QCoreApplication::applicationDirPath().toStdString();
#else
////***** OS specific: on windows
//    DWORD v;
//    size_t size = 1;
//    do {
//        size += MAX_PATH;
//        space.resize(int(size));
//        v = GetModuleFileName(NULL, space.data(), DWORD(space.size()));
//    } while (Q_UNLIKELY(v >= size));
////***** OS specific: on windows

    // Check
    if( !pPath.length() )
        LzLogException("", "LzServices::StartUpPath() requires an argument to determine start-up path!")

    // Find slash or back-slash separator
    std::size_t lLastSlash = pPath.rfind("/");
    std::size_t lLastBackS = pPath.rfind("\\");

    // Found any path marker?
    if( lLastSlash==std::string::npos && lLastBackS==std::string::npos )
    {
//        LzLogException("", "Could not find neither a / nor a \\ in path '"<<pPath<<"'!")

        // Return local folder
        return ".";
    }
    else
    {
        // Crop string
        if( lLastSlash != std::string::npos )
            return pPath.substr(0, lLastSlash);
        else
            return pPath.substr(0, lLastBackS);
    }
#endif
}

//================================================================================
string WorkingPath()
{
#if defined(USING_QT)
    return QDir::currentPath().toStdString();
#else
//    // Check
//    if( !pPath.length() )
//        LzLogException("", "LzServices::WorkingPath() requires an argument to determine working directory!")
    /*
     * https://www.tutorialspoint.com/find-out-the-current-working-directory-in-c-cplusplus
    #ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
    #else
    #include <unistd.h>
    #define GetCurrentDir getcwd
    #endif

    #include<iostream>
    using namespace std;

    std::string get_current_dir() {
       char buff[FILENAME_MAX]; //create string buffer to hold path
       GetCurrentDir( buff, FILENAME_MAX );
       string current_working_dir(buff);
       return current_working_dir;
    }

    main() { cout << get_current_dir() << endl; }
     *
     */

    LzLogException("", "LzServices::WorkingPath() is not available!")
#endif
}

//================================================================================
void CreateOrPurgeDir( const string & pDir )
{
#ifdef USING_Cpp_17
    LOCK_API

    // Check existence
    if( std::filesystem::exists(pDir) )
    {
        // Log
        LzLogM("", "Recursively purging directory '"<<pDir<<"'")

        // Dir exists: purge
        std::filesystem::remove_all(pDir);
    }

    // Log
    LzLogM("", "Creating directory '"<<pDir<<"'")

    // Create (or re-create) dir
    try
    {
        if( !std::filesystem::create_directory(pDir) )
            LzLogErr("", "Could not create directory '"<<pDir<<"'!")
    }
    catch( const std::exception & pExc )
    {
        // Log underlying exception and rethrow with higher semantics
        LzLogErr("", "Exception: "<<pExc.what())
        LzLogException("", "Could not create directory '"<<pDir<<"'!")
    }
#elif defined(USING_QT)
    QString lDir( pDir.c_str() );

    // Check existence
    if( QDir(lDir).exists() )
    {
        // Log
        LzLogM/*sg*/("", "Purging recursively directory: "<<lDir.toStdString())

        // Dir exists: purge
        QDir(lDir).removeRecursively();
    }

    // Log
    LzLogM/*sg*/("", "Creating directory: "<<lDir.toStdString())

    // Create (or re-create) dir
    if( !QDir().mkpath(lDir) )
        LzLogErr("", "QDir::mkpath failed!")
#else
    LzLogException("", "LzServices::CreateOrPurgeDir("<<pDir<<") is not available!")
#endif
}

//================================================================================
bool FileOrDirExists( const string & pFile )
{
    struct stat buf;
    if( stat(pFile.c_str(), &buf) != -1 )
        return true;

    return false;
}

//================================================================================
void SplitFile( const string & pFile, string & pPath, string & pFilename )
{
    // Clear previous
    pPath = "";
    pFilename = "";

    // Split string into items using "."
    vector<string> lItems;
    SplitString( pFile, { "/", "\\" }, lItems ); // All items are non-empty iff pFile is non-empty

    // No items iff pFile is empty
    if( lItems.size() == 0 )
        return;

    // File ends with a slash or back-slash?
    if( pFile.back()=='/' || pFile.back()=='\\' )
    {
        // The file contains just a path
        pPath = pFile;
        // pFilename remains ""
    }
    else
    {
        // HERE: Last item is the filename

        // Filename is the last item
        pFilename = lItems.back();
        // pPath is still "" here

        // Do we have a path?
        if( lItems.size() > 1 )
        {
            // Path is everything before filename
            pPath = pFile.substr( 0, pFile.length()-pFilename.length() );
        }
    }
}

//================================================================================
/**
 * @brief These functions are not published (only used internally)
 *
 * SplitFilename has a similar name to SplitFile and the exact same prototype
 * ==> high risk of confusion
 */
static unsigned int GetExtensionLevel( const string & pStr )
{
    // Integers are accepted as extensions if extension level is set to 2
    if( IsInt(pStr) )
        return 2;

    // Integers followed by nun num caracters are accepted as extensions if
    // extension level is set to 1
    if( IsIntAlpha(pStr) )
        return 1;

    // Any other strings are accepted as extensions since extension level is
    // always >= 0
    return 0;
}
//
static unsigned int sExtensionToleranceLevel = 0;   // Alpha-int ("+123 meters^2") (and hence int ("+123")) extensions are rejected
//static unsigned int sExtensionToleranceLevel = 1; // Int ("+123") extensions are rejected; alpha-int extensions are accepted
//static unsigned int sExtensionToleranceLevel = 2; // Any string is accepted as an extension
//
static void SplitFilename( const string & pFilename, string & pNameNoExt, string & pExt )
{
    // Check if we have a filename
    if( pFilename.length() == 0 )
        LzLogException("", "An empty filename cannot be split!")

    // HERE: lFilename is non-empty

    // Split string into items using "."
    vector<string> lItems;
    SplitString( pFilename, { "." }, lItems );

    // All items are non-empty, and at least one item since lFilename is non-empty

    // No extension at all?
    if( lItems.size() == 1 )
    {
        pNameNoExt = lItems[0];
        pExt = "";
    }
    else
    {
        // HERE: at least 2 items, the filename without extension, and the (possibly many) extensions

        // Set filename without extension
        pNameNoExt = lItems.front();

        // Collect extensions
        pExt = "";
        for( size_t i=1 ; i<lItems.size() ; i++ )
        {
            // Is it an acceptable extension?
            // Extension level must not be higher than the tolerance level
            if( GetExtensionLevel(lItems[i]) > sExtensionToleranceLevel )
            {
                // Reject extension!

                // Any other extensions so far?
                if( pExt.length() )
                {
                    // All extensions so far go to filename without ext
                    pNameNoExt += "." + pExt;

                    // Purge extensions
                    pExt = "";
                }

                // Add fake extension
                pNameNoExt += "." + lItems[i];
            }
            else
            {
                // Need a '.'?
                if( pExt.length() )
                    pExt += ".";

                // Add extension
                pExt += lItems[i];
            }
        }
    }
}

//================================================================================
string GetPath( const string & pFile )
{
    // Split file
    string lPath;
    {
        string lFilename_unused;
        SplitFile( pFile, lPath, lFilename_unused );
    }

    return lPath;
}

//================================================================================
void SetExtensionLevel( unsigned int pLevel )
{
    // Check
    if( pLevel > 2 )
        LzLogException("", "Unexpected extension tolerance level! Should be 0, 1 or 2.")

    // Log
    LzLogM("", "Setting extension tolerance level to "<<pLevel<<".")

    // Set
    sExtensionToleranceLevel = pLevel;
}

//================================================================================
void ResetExtensionLevel()
{
    // Log
    LzLogM("", "Resetting extension tolerance level to 0.")

    // Reset
    sExtensionToleranceLevel = 0;
}

//================================================================================
string GetExtension( const string & pFile )
{
    // Split file
    string lFilename;
    {
        string lPath_unused;
        SplitFile( pFile, lPath_unused, lFilename );
    }

    // Split filename
    string lNameNoExt;
    string lExt;
    SplitFilename( lFilename, lNameNoExt, lExt );

    return lExt;
}

//================================================================================
string GetFilenameWithoutExt( const string & pFile )
{
    // Split file
    string lFilename;
    {
        string lPath_unused;
        SplitFile( pFile, lPath_unused, lFilename );
    }

    // Split filename
    string lNameNoExt;
    string lExt;
    SplitFilename( lFilename, lNameNoExt, lExt );

    return lNameNoExt;
}

//================================================================================
bool VerifyExtension( const string & pFile, const string & pExtension )
{
    // Get extension from file
    // Might be an empty string
    // Will throw exception if does not end with a filename
    string lFileExt = GetExtension( pFile );

    // Process trivial cases

    // Impossible fit
    if( lFileExt.length() < pExtension.length() )
        return false;

    // No extension
    if( lFileExt=="" || pExtension=="" )
        return lFileExt == pExtension;

    // Convert to uppercase
    std::transform( lFileExt.begin(), lFileExt.end(), lFileExt.begin(), ::toupper );
    //
    string lUpExt = pExtension;
    std::transform( lUpExt.begin(), lUpExt.end(), lUpExt.begin(), ::toupper );

    // HERE: lFileExt and lUpExt are both non-empty,
    //       and lFileExt.length() >= lUpExt.length()

    // FileExt must end with UpExt
    string::size_type lPos = lFileExt.rfind( lUpExt ); // makes sense since both non-empty, can be npos
    string::size_type lExpectedPos = lFileExt.length() - lUpExt.length(); // >= 0
    //
    return lPos == lExpectedPos;
}

//================================================================================
string ChangeExtension( const string & pFile, const string & pNewExtension )
{
    // Split file
    string lPath;
    string lFilename;
    SplitFile( pFile, lPath, lFilename );

    // Split filename
    string lNameNoExt;
    string lExt;
    SplitFilename( lFilename, lNameNoExt, lExt );

    // Return file with the new extension
    if( pNewExtension == "" )
    {
        // Remove extension
        return lPath + lNameNoExt;
    }
    else
    {
        // Check if a dot is already provided
        if( pNewExtension[0] == '.' )
            return lPath + lNameNoExt + pNewExtension;
        else
            return lPath + lNameNoExt + "." + pNewExtension;
    }
}
#pragma endregion


#pragma region "String ops"
//================================================================================
bool StartsWidth( const string & pStr, const string & pHeader )
{
    return pStr.find( pHeader ) == 0;
}

//================================================================================
int StrToInt( const string & pStr )
{
    size_t lIdx;
    int lRes;

    try { lRes = std::stoi( pStr, &lIdx ); }
    catch( ... ) { LzLogException("", "Syntax error in '"<<pStr<<"'!")  }

    // Check
    if( lIdx != pStr.length() )
        LzLogException("", "Syntax error in '"<<pStr<<"'!!")

    return lRes;
}

//================================================================================
bool IsInt( const string & pStr )
{
    // Disabling and enabling log
    Detonator lDet([]{ DisableLog(); }, []{ EnableLog(); });

    try
    {
        // Try to parse string EXACTLY as an int
        StrToInt( pStr );

        // Success
        return true;
    }
    catch(...)
    {
        // Failure
        return false;
    }
}

//================================================================================
bool IsIntAlpha( const string & pStr )
{
    // Disabling and enabling log
    Detonator lDet([]{ DisableLog(); }, []{ EnableLog(); });

    try
    {
        // Try to parse string as STARTING WITH an int
        std::stoi( pStr );

        // Success
        return true;
    }
    catch(...)
    {
        // Failure
        return false;
    }

}

//================================================================================
size_t StrToUnsInt( const string & pStr )
{
#if 1
    std::size_t lIdx;
//    unsigned long lRes;
    size_t lRes;

    try { lRes = std::stoul( pStr, &lIdx ); }
    catch( ... ) { LzLogException("", "Syntax error in '"<<pStr<<"'!")  }

    // Check
    if( lIdx != pStr.length() )
        LzLogException("", "Syntax error in '"<<pStr<<"'!!")

//    return static_cast<unsigned int>( lRes );
    return lRes;
#else
    // Check
    const char * lpChar = pStr.c_str();
    while( *lpChar )
        if( *(lpChar++) == '-' )
            LzLogException("", "Found a '-' in string! Cannot be a proper unsigned.")

    // Parse
    std::size_t lIdx;
    unsigned long lLongRes = std::stoul( pStr, &lIdx );

    // Check
    if( lIdx != pStr.length() )
        LzLogException("", "Syntax error!")

    // Check
    if( lLongRes > std::numeric_limits<unsigned int>::max() )
        LzLogException("", "Out of range!")

    return (unsigned int)lLongRes;
#endif
}

//================================================================================
double StrToDouble( const string & pStr )
{
    std::size_t lIdx;
    double lRes;

    try { lRes = std::stod( pStr, &lIdx ); }
    catch( ... ) { LzLogException("", "Syntax error in '"<<pStr<<"'!")  }

    // Check
    if( lIdx != pStr.length() )
        LzLogException("", "Syntax error in '"<<pStr<<"'!!")

    return lRes;
}

//================================================================================
size_t FindNth( const string & haystack, size_t pos, const string & needle, size_t nth )
{
    size_t found_pos = haystack.find( needle, pos );
    if (1 == nth || string::npos == found_pos)
        return found_pos;

    return FindNth(haystack, found_pos + 1, needle, nth - 1);
}

//================================================================================
string FormatNumber( double pX, size_t pDecimals )
{
    stringstream ss;
    ss << std::fixed << std::setprecision(pDecimals) << pX;
    return ss.str();
}

//================================================================================
double RoundNumber( double pX, size_t pDecimals )
{
    const double multiplier = std::pow(10.0, pDecimals);
    return std::nearbyint(pX * multiplier) / multiplier;
}

//================================================================================
void LTrim( string & s )
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

//================================================================================
void RTrim( string & s )
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base()/*base converts a reverse iterator into the corresponding forward iterator*/, s.end());
}

//================================================================================
void Trim( std::string & s )
{
    LTrim(s);
    RTrim(s);
}

//================================================================================
string GetTrimmed( const string & s )
{
    string lRes = s;
    Trim( lRes );
    return lRes;
}

//================================================================================
std::string ReplaceAll( const string & pPattern, const string & pBy, const string & pIn )
{
    // Look for pattern
    size_t lFound = pIn.find( pPattern );

    // Nothing found
    if( lFound == string::npos )
        return pIn;

    // Found pattern: split string, replace, and stitch together
    string lLeft  = pIn.substr( 0, lFound );
    string lRight = pIn.substr( lFound+pPattern.length(), string::npos );
    return lLeft + pBy + ReplaceAll( pPattern, pBy, lRight );
}

//================================================================================
void SplitString( const string & pStr, const vector<string> & pSeps, vector<string> & pItems, bool pKeepEmptyParts/*=false*/ )
{
    // Clear previous
    pItems.clear();

    // If string is empty, my work here is done
    if( pStr == "" )
        return;

    // Check if any separators
    if( pSeps.size() == 0 )
        LzLogException("", "Cannot split a string without any separator!")

    // Check that all separators have non 0 length
    vector<size_t> lSepLens( pSeps.size() );
    for( unsigned int s=0 ; s<pSeps.size() ; s++ )
    {
        // Compute length and check
        lSepLens[s] = pSeps[s].length();
        if( lSepLens[s] == 0 )
            LzLogException("", "Cannot split string '"<<pStr<<"' using empty separator!")
    }

    // Split
    size_t lCurrPos = 0;
    while( true )
    {
        // Find first separator
        size_t lSepPos = string::npos;
        unsigned int lSepIdx = 0; // init to 0 to silence a warning
        for( size_t s=0 ; s<pSeps.size() ; s++ )
        {
            size_t lPos = pStr.find(pSeps[s], lCurrPos);
            if( lPos != string::npos )
            {
                // Update
                if( lSepPos==string::npos || lPos<lSepPos )
                {
                    lSepPos = lPos;
                    lSepIdx = s;
                }
            }
        }

        // Found any?
        if( lSepPos != string::npos )
        {
            // Item length
            const size_t lLen = lSepPos - lCurrPos;

            // Extract only non empty items?
            if( pKeepEmptyParts )
            {
                // Keep whether empty or not
                pItems.push_back( pStr.substr(lCurrPos, lLen) );
            }
            else
            {
                // Keep non-empty only
                if( lLen > 0 )
                    pItems.push_back( pStr.substr(lCurrPos, lLen) );
            }

            // Skip separator by moving cursor
            lCurrPos = lSepPos + lSepLens[lSepIdx];
        }
        else
        {
            // Stash remainder, if not empty
            if( lCurrPos < pStr.length() )
                pItems.push_back( pStr.substr(lCurrPos) );

            // Finished
            break;
        }
    }
}

//================================================================================
string FormatTimestamp( std::chrono::time_point<std::chrono::system_clock> p_TS )
{
    LOCK_API

    // Why is there no C++11 threadsafe alternative to std::localtime and std::gmtime?
    // https://howardhinnant.github.io/date/date.html
    // https://stackoverflow.com/questions/25618702/why-is-there-no-c11-threadsafe-alternative-to-stdlocaltime-and-stdgmtime
    // https://stackoverflow.com/questions/10894504/is-there-any-stdchrono-thread-safety-guarantee-even-with-multicore-context

    const time_t l_TimeT = std::chrono::system_clock::to_time_t(p_TS);
    const auto l_MSecs = std::chrono::duration_cast<std::chrono::milliseconds>(p_TS.time_since_epoch()) % 1000;

    // https://en.cppreference.com/w/c/chrono/tm
    struct tm l_Tm;
    // localtime_r(&l_Tm, &l_TimeT);

    std::stringstream l_SS;
    l_SS << std::setfill('0')
         << std::setw(4) << (1900 + l_Tm.tm_year) << "-"
         << std::setw(2) << (l_Tm.tm_mon + 1)     << "-"
         << std::setw(2) << l_Tm.tm_mday          << " - "
         << std::setw(2) << l_Tm.tm_hour << "-"
         << std::setw(2) << l_Tm.tm_min  << "-"
         << std::setw(2) << l_Tm.tm_sec  << "-"
         << std::setw(3) << l_MSecs.count();

//*** Alternative formatting
//    stringstream lTimestamp;
//    lTimestamp << std::put_time(std::localtime(&lTimeT), "%H:%M:%S")
//               << '.' << std::setfill('0') << std::setw(3) << lMSecs.count() << " ";
//*** Alternative formatting

    return l_SS.str();
//#endif
}

//================================================================================
string FormatTimestampNow()
{
    // Get current time
    // Yes, calls to some_clock::now() from different threads should be thread safe.
    // https://stackoverflow.com/questions/10894504/is-there-any-stdchrono-thread-safety-guarantee-even-with-multicore-context
    // Fingers crossed... not a sensitive item though (log only)
    std::chrono::time_point<std::chrono::system_clock> lNow = std::chrono::system_clock::now();

    // Format timestamp
    return FormatTimestamp( lNow );
}
#pragma endregion


#pragma region "Scope terminator"
//================================================================================
Detonator::Detonator( Type pType/*=Type::Defusable*/ )
    : mType(pType), mDest(nullptr)
{}

//================================================================================
Detonator::Detonator( std::function<void()> pDest, Type pType/*=Type::Defusable*/ )
    : mType(pType), mDest(pDest)
{}

//================================================================================
Detonator::Detonator( std::function<void()> pConst, std::function<void()> pDest, Type pType/*=Type::Defusable*/ )
    : mType(pType), mDest(pDest)
{
    // Exceptions are allowed in ctor
    pConst();
}

//================================================================================
Detonator::~Detonator()
{
    try
    {
        // Exceptions are NOT allowed in dtor
        if( mDetonate && mDest )
            mDest();
    }
    catch( const std::exception & pExc )
    {
        LzLogNode("", "*** Caught an exception in Detonator::~Detonator!")
        LzLogE("", "*** Exception: "<<pExc.what())
        LzLogE("", "*** This exception is silenced as no exceptions are allowed in destructors...")
    }
    catch( ... )
    {
        LzLogNode("", "*** Caught an unknown exception in Detonator::~Detonator!")
        LzLogE("", "*** This exception is silenced as no exceptions are allowed in destructors...")
    }

    // On exceptions from Destructor and stack unwinding
    //
    // https://www.cs.technion.ac.il/users/yechiel/c++-faq/dtors-shouldnt-throw.html
    // https://scc.ustc.edu.cn/zlsc/sugon/intel/ssadiag_docs/pt_reference/references/sc_cpp_ex_from_delete.htm
    // https://stackoverflow.com/questions/130117/throwing-exceptions-out-of-a-destructor
}

//================================================================================
void Detonator::set( std::function<void()> pDest )
{
    mDest = pDest;
}

//================================================================================
void Detonator::defuse()
{
    // Check
    if( mType == Type::NonDefusable )
        LzLogException("", "This detonator cannot be defused! Run for cover!")

    // Defuse
    mDetonate = false;
}
#pragma endregion


#pragma region "Timer"
//================================================================================
static mutex sChronoMex;

//================================================================================
std::chrono::time_point<std::chrono::steady_clock> Chrono_Now()
{
    std::unique_lock<mutex> lLock(sChronoMex);
    return std::chrono::steady_clock::now();
}

//================================================================================
double Chrono_MSecs( std::chrono::time_point<std::chrono::steady_clock> pStart )
{
    std::unique_lock<mutex> lLock(sChronoMex);
    return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now() - pStart).count();
}
#pragma endregion


//================================================================================
void StdoutToFile( const string & pFile )
{
#if defined(OS_WINDOWS)
    FILE * lLogFile = nullptr;
    auto lErr = freopen_s(&lLogFile, (LzServices::StartUpPath()+"/"+pFile).c_str(), "w+", stdout);
  
    if( lErr != 0 )
        LzLogException("", "Could not redirect stdout to file '"<<pFile<<"'!")
#elif defined(OS_LINUX)
    FILE * lLogFile = nullptr;
    FILE * lpNew = freopen((LzServices::StartUpPath()+"/"+pFile).c_str(), "w+", stdout);
  
    if( lpNew == nullptr )
	    LzLogException("", "Could not redirect stdout to file '"<<pFile<<"'!")
#else
	*** StdoutToFile not defined! ***
#endif
}

//================================================================================
//static std::ofstream sLogFile; // Protected by sLogMex / LOCK_API
static std::ofstream & GetLogFile()
{
    static std::ofstream * spLogFile = nullptr;
    if( !spLogFile )
        spLogFile = new std::ofstream;

    return *spLogFile;
}
#define sLogFile GetLogFile()

//================================================================================
static std::thread::id sMainThreadId; // Protected by sLogMex / LOCK_API
const std::thread::id & MainThreadId() { return sMainThreadId; }

//================================================================================
enum class LogMode { TXT, YAML };
LogMode sLogMode; // Protected by sLogMex / LOCK_API

//================================================================================
void OpenLogFile( const string & pFile )
{
    LOCK_API

    // Check
    if( sLogFile.is_open() )
        LzLogException("", "Cannot open log file '"<<pFile<<"'! Log file already opened.")

    // Open log file
    sLogFile.open( pFile.c_str(), std::ofstream::out | std::ofstream::trunc );

    // Check
    if( !sLogFile.is_open() )
        LzLogException("", "Could not open log file '"<<pFile<<"'!")

    // Store main thread id
    sMainThreadId = std::this_thread::get_id();

    // Check extension
    if( VerifyExtension(pFile, "txt") )
        sLogMode = LogMode::TXT;
    else
    if( VerifyExtension(pFile, "yaml") )
        sLogMode = LogMode::YAML;
    else
        LzLogException("", "Unexpected log file extension '"<<GetExtension(pFile)<<"'!")

    // Log
    LzLogM("", "<<< Log file opened by LzServices::OpenLogFile >>>")
    LzLogM("", "<<< Log file= '"<<pFile<<"' >>>")
}

//================================================================================
void CloseLogFile()
{
    LOCK_API

    // Log
    LzLogM("", "<<< Log file closed >>>")

    // Close
    sLogFile.close();
}

//================================================================================
static atomic_bool sLogEnabled( true );
//
static List<string> sLogTail;
static unsigned int sLogTailLength = 10;

//================================================================================
void EnableLog()
{
    sLogEnabled = true;
}

//================================================================================
void DisableLog()
{
    sLogEnabled = false;
}

//================================================================================
List<string> LogTail()
{
    LOCK_API

    return sLogTail;
}

//================================================================================
void SetLogTailLength( unsigned int pLen )
{
    LOCK_API

    // Set max length
    sLogTailLength = pLen;

    // Adjust content
    while( sLogTail.Count() > sLogTailLength )
        sLogTail.DelHead();
}

//================================================================================
void LogLine( const string & pLine )
{
#if 1
    // NOP if log disabled
    if( !sLogEnabled )
        return;

    LOCK_API

    // Format timestamp
    const string lTimestamp = FormatTimestampNow();

    // Line
    string lLogLine;

    // Mode = TXT
    if( sLogMode == LogMode::TXT )
    {
        // Init with timestamp
        lLogLine = lTimestamp;

        // Add tab string
        lLogLine += LzLogNodeObject::TabString();

        // Actual line
        lLogLine += pLine;
    }
    else
    // Mode = YAML
    if( sLogMode == LogMode::YAML )
    {
        // Init with space string
        lLogLine = LzLogNodeObject::SpaceString();

        // Add timestamp
        lLogLine += lTimestamp;

        // Actual line
        lLogLine += pLine;
    }
    else
        LzLogException("", "Unexpected log mode!")

    // Log line
    if( sLogFile.is_open() )
        sLogFile << lLogLine << endl;

    // To console
    cout << lLogLine << endl << flush;

    // To tail
    sLogTail.AddTail( lLogLine );

    // Prune head
    if( sLogTail.Count() > sLogTailLength )
        sLogTail.DelHead();

////*** DEBUG
//std::cout << ".... TAIL LEN= " << sLogTail.Count() << std::endl;
//BrowseList(iL, sLogTail )
//    std::cout << "---- " << sLogTail.GetAt(iL) << std::endl;
////*** DEBUG
#else
    // NOP if log disabled
    if( !sLogEnabled )
        return;

    LOCK_API

    // Format timestamp
    const string lTimestamp = FormatTimestampNow();

    // Mode = TXT
    if( sLogMode == LogMode::TXT )
    {
//        LzServices::LogFile() << lTimestamp;

//        // Add tab string
//        LzServices::LogFile() << LzLogNodeObject::TabString();

        // Log line
        LzServices::LogFile() << lTimestamp /*<< " "*/ << pLine << endl;

        // To console
        cout << lTimestamp /*<< " "*/ << pLine << endl << flush;

        // To tail
        sLogTail.AddTail( pLine );
    }
    else
    // Mode = YAML
    if( sLogMode == LogMode::YAML )
    {


//        zzee
    }
    else
        LzLogException("", "Unexpected log mode!")

    // Crop tail
    if( sLogTail.Count() > sLogTailLength )
        sLogTail.DelHead();

//*** DEBUG
std::cout << ".... TAIL LEN= " << sLogTail.Count() << std::endl;
BrowseList(iL, sLogTail )
    std::cout << "---- " << sLogTail.GetAt(iL) << std::endl;
//*** DEBUG
#endif
}

/*================================================================================
 *void StdoutToLogFile( const string & pDir )
 *{
 *	//std::time_t lTime = std::time(nullptr);
 *	////std::chrono::steady_clock
 *	//cout << "Time= " << lTime << endl;
 *	//cout << "Time= " << std::asctime(std::localtime(&lTime)) << endl;
 *
 *	//auto lNow = std::chrono::system_clock::now();
 *	//std::time_t lST = std::chrono::system_clock::to_time_t( lNow );
 * //	cout << "SysClock= " << lST << endl;
 *
 *	//std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
 * //   std::chrono::system_clock::duration tp = now.time_since_epoch();
 *
 * //   tp -= std::chrono::duration_cast<std::chrono::seconds>(tp);
 *
 * //   std::time_t tt = std::chrono::system_clock::to_time_t(now);
 *
 * //   print_time(*std::gmtime(&tt), tp);
 *	////for( int i=0 ; i<20 ; i++ )
 *	//{
 *	//	std::this_thread::sleep_for(std::chrono::seconds(2));
 *	//    print_time(*std::localtime(&tt), tp);
 *	//}
 *    string lFilePath = "./_tmp_TEST_LOG_FILE.txt";
 *}
 */

//================================================================================
template <typename Duration> void print_time( tm t, Duration fraction )
{
//https://stackoverflow.com/questions/27136854/c11-actual-system-time-with-milliseconds
	
	using namespace std::chrono;
    std::printf("[%04u-%02u-%02u %02u:%02u:%02u.%03u]\n", t.tm_year + 1900,
                t.tm_mon + 1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec,
                static_cast<int>(fraction / milliseconds(1)));

// VS2013's library has a bug which may require you to replace
// "fraction / milliseconds(1)" with
// "duration_cast<milliseconds>(fraction).count()"
}

//================================================================================
void GetTime( int & /*pYear*/, int & /*pMonth*/, int & /*pDay*/, int & /*pHour*/, int & /*pMinute*/, int & /*pSecond*/, int & /*pMillisecond*/ )
{
    //*** TODO
}


#if defined(USING_LzLog) && defined(USING_QT) && !defined(NO_GUI)
//================================================================================
unsigned int LzWaitCursor::sNested = 0;

//================================================================================
LzWaitCursor::LzWaitCursor()
{
    if( sNested == 0 )
        QApplication::setOverrideCursor(Qt::WaitCursor);

    sNested++;
}

//================================================================================
LzWaitCursor::~LzWaitCursor()
{
    sNested--;

    if( sNested == 0 )
        QApplication::restoreOverrideCursor();
}
#endif


#ifdef USING_LzLog
/* TEST LzLogNode
 *
    bool lEnable1 = true;
    bool lEnable2 = true;
    {
        LzOptLogN(true, "", "Coucou 1")
        LzLogM("", "Hello 1")
    }
    LzLogMsg("", "Hey")
    {
        LzOptLogN(lEnable1, "", "Coucou 2")
        if( lEnable1 )
            LzLogM("", "Hello 2")

        {
            LzOptLogN(lEnable2, "", "Inner Coucou 2")
            if( lEnable2 )
                LzLogM("", "Inner Hello 2")
        }
    }
    LzLogMsg("", "Heyyy")
    {
        LzOptLogN(true, "", "Coucou 3")
        LzLogM("", "Hello 3")
    }
    return 0;
 *
 */
//================================================================================
#if defined(EXPORT_LzServices) && defined(__MSVC_LIB__)
/*thread_local not allowed in MSVC libs*/ unsigned int LzLogNodeObject::sTabs = 0;
#else
thread_local DLL_EXPORT_LzServices unsigned int LzLogNodeObject::sTabs = 0;
#endif

//================================================================================
LzLogNodeObject::LzLogNodeObject( bool pEnabled, bool pTimer, const string & pMsg/*=""*/, const string & pThreadInfo/*=""*/ )
 : mEnabled(pEnabled)
{
    if( mEnabled )
    {
        sTabs++;

        mMsg = pMsg;

        mTimer = pTimer;
        if( mTimer )
        {
            mThreadInfo = pThreadInfo;
            mStart = LzServices::Chrono_Now();
        }
    }
}

//================================================================================
LzLogNodeObject::~LzLogNodeObject()
{
    if( mEnabled )
    {
        sTabs--;

        if( mTimer )
            LINKZ_LOG_TXT("\\-- "<<mThreadInfo<<" "<<mMsg<<" completed in "<<LzServices::Chrono_MSecs(mStart)<<" msec")
        else
            LINKZ_LOG_TXT("\\-- "<<mThreadInfo<<" "<<mMsg)
    }
}

//================================================================================
string LzLogNodeObject::TabString()
{
    string lRet = "";
    for( unsigned int t=0 ; t<sTabs ; t++ )
        lRet += "| ";

	return lRet;
}

//================================================================================
string LzLogNodeObject::SpaceString()
{
    string lRet = "";
    for( unsigned int t=0 ; t<sTabs ; t++ )
        lRet += "  ";

    return lRet;
}

#endif
}
