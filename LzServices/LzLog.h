#pragma once

#include "LzLib_Services.h"
#include <fstream>
#include <thread>
#if defined(USING_QT) && !defined(NO_GUI)
    #include <QApplication>
    #include <QMessageBox>
#endif
#include <sstream>
#include <random>
#include <functional>


// DLL exports
namespace LzServices
{
    class DLL_EXPORT_LzServices LzLogNodeObject;
}


namespace LzServices
{
template <class T> class List;

using std::string;
using std::vector;
using std::ifstream;


// Min and max
template<class T> T Min( T pA, T pB ) { return pA < pB ? pA : pB ; }
template<class T> T Max( T pA, T pB ) { return pA > pB ? pA : pB ; }


#pragma region "Data subsampling"
// Extracting an approximate number of samples
//
// Returns a smaller step, possibly yielding more samples than desired
// Increment step by 1 to get less samples than desired
// If pSampleSize | pOriginalSize, returns a step that yields exactly
// the desired number of samples
size_t GetSubSampSmallStep( size_t pOriginalSize, size_t pSampleSize );

// Check the actual number of samples for a given sample step
size_t GetActualSamplesCount( size_t pOriginalSize, size_t pSampleStep );

// Extracting an exact number of samples
//
// Returns the set of indices of the desired size, which samples uniformely the data
void GetSubSampIdx( size_t pOriginalSize, size_t pSampleSize, vector<size_t> & pTo );
template<class T> void SubSampleTo( const vector<T> & pFrom, vector<T> & pTo, size_t pSampleSize )
{
    // Get sampled indices
    vector<size_t> lIdx;
    GetSubSampIdx( pFrom.size(), pSampleSize, lIdx );

    // Copy sampled elements
    pTo.resize( pSampleSize );
    for( size_t i=0 ; i<lIdx.size() ; i++ )
        pTo[i] = pFrom[ lIdx[i] ];
}
//
void TEST_GetSubSampIdx();
void TEST_GetSubSampIdx_2();
#pragma endregion


#pragma region "Read Write using ifstream"
// Reading Double, Float and UInt values from binary ifstream
//
// Use pSwap to swap byte ordering
//
double readDouble64( ifstream & pFrom, bool pSwap=false );
float readFloat32( ifstream & pFrom, bool pSwap=false );
unsigned int readUnsignedInt32( ifstream & pFrom, bool pSwap=false );

// Read file content to memory
unsigned char * readBinFile( const string & pFileName, unsigned int & pSize );

// Read text file to string
void ReadFileToString( const string & pFile, string & pOut );

// Reading values from an array
double readDouble64( const unsigned char * pFrom, unsigned int & pIdx, bool pSwap=false );
float readFloat32( const unsigned char * pFrom, unsigned int & pIdx, bool pSwap=false );
unsigned int readUnsignedInt32( const unsigned char * pFrom, unsigned int & pIdx, bool pSwap=false );
#pragma endregion


#pragma region "Local IP"
//*** TODO
//https://stackoverflow.com/questions/122208/get-the-ip-address-of-local-computer
//*** TODO
#pragma endregion


#pragma region "Base64 encoding"
//*** TODO
// LzServices::base64Encode(...) --> https://stackoverflow.com/questions/6385319/how-to-convert-a-binary-string-into-base64-encoded-data
//*** TODO
#pragma endregion


#pragma region "Logical cores"
    /* Count available logical processors
     *
     * This function will return the number of processor cores + cores that behave as two or more logical cores due to hyperthreading.
     * May return 0 when not able to detect.
     *
     * You should use this value from your target platform to plan the max number of Runnable threads your program should concurrently use.
     * You can also designate a core for all your waitable threads and use the remaining number of cores for runnable threads. For example,
     * on a quad-core system, use one core for ALL waitable threads and use three runnable threads for the remaining three cores. Depending
     * on your thread schedulers efficiency, a few of your runnable threads might get context switched out (due to page faults etc.) leaving
     * the core idle for some amount of time. If you observe this situation during profiling, you should create a few more runnable threads
     * than the number of your cores and tune it for your system.
     */
    DLL_EXPORT_LzServices unsigned int CountLogicalCores();
#pragma endregion


#pragma region "Entropy"
	// Random generator
    DLL_EXPORT_LzServices std::mt19937 & RandomEngine();
#pragma endregion


#pragma region "Dir / file operations"
    // Current exe directory (not including last "/")
    // If not using Qt, caller must provide the path to the current exe (main's argv[0])
    DLL_EXPORT_LzServices string StartUpPath( const string & pPath="" );

    //Working directory (not including last "/")
    DLL_EXPORT_LzServices string WorkingPath();

    // Create or purge directory
    DLL_EXPORT_LzServices void CreateOrPurgeDir( const string & pDir );

    /**
     * Definition: File = Path + Filename
     * ===========
     *
     * Here using '/'
     *
     *      File = 'D:/Users/Marek/Desktop/Big_Lebowsky.mp4'
     *      Path = 'D:/Users/Marek/Desktop/' or 'D:/Users/Marek/Desktop' if final '/' is implicit
     * NB: all functions below return paths that are either empty (no path) or ending with '/'
     *
     *  Filename = 'Big_Lebowsky.mp4'
     * Extension = 'mp4'
     *
     * Filename without extension = 'Big_Lebowsky'
     *
     * NB: an extension cannot be of a purely numerical type e.g., '.+123'. A file thus cannot
     * === have this type of extension i.e., in 'file.-123' there is no extension and the
     *     filename without extension is 'file.-123'.
     *
     */

    // Checks existance of file
    DLL_EXPORT_LzServices bool FileOrDirExists( const string & pFile );

    // Splits full path into path to containing directory, and filename
    DLL_EXPORT_LzServices void SplitFile( const string & pFile, string & pPath, string & pFilename );

    // Get path from full path
    DLL_EXPORT_LzServices string GetPath( const string & pFile );

    // Extension tolerance level
    DLL_EXPORT_LzServices void SetExtensionLevel( unsigned int pLevel );
    DLL_EXPORT_LzServices void ResetExtensionLevel();

    // Get extension from full path
    DLL_EXPORT_LzServices string GetExtension( const string & pFile );

    // Get filename without extension from full path
    DLL_EXPORT_LzServices string GetFilenameWithoutExt( const string & pFile );

    // Extension manipulation
    DLL_EXPORT_LzServices bool VerifyExtension( const string & pFile, const string & pExtension );
    DLL_EXPORT_LzServices string ChangeExtension( const string & pFile, const string & pNewExtension );
#pragma endregion


#pragma region "String ops"
    // Computes char length of an integer or double. "-123" ==> 4, "0" ==> 1, "(+)123" ==> 3
    template <typename T> size_t StringLength( T pT )
    {
        return std::to_string(pT).length();
    }

    // Starts width
    bool StartsWidth( const string & pStr, const string & pHeader );

    // Parsing numerical values from strings
    DLL_EXPORT_LzServices int StrToInt( const string & pStr );
    DLL_EXPORT_LzServices bool IsInt( const string & pStr );
    DLL_EXPORT_LzServices bool IsIntAlpha( const string & pStr );
    DLL_EXPORT_LzServices size_t StrToUnsInt( const string & pStr );
    DLL_EXPORT_LzServices double StrToDouble( const string & pStr );
    DLL_EXPORT_LzServices size_t FindNth(const string & haystack, size_t pos, const string & needle, size_t nth);
    template <typename T> string ToString( const T & t )
    {
        std::locale::global( std::locale("C") );
        string str{std::to_string (t)};
        size_t offset{1};
        if (str.find_last_not_of('0') == str.find('.'))
        {
            offset = 0;
        }
        str.erase(str.find_last_not_of('0') + /*(unsigned int)*/offset, string::npos);
        return str;
    }

    // Format number
    string FormatNumber( double pX, size_t pDecimals );

    // Round number
    double RoundNumber( double pX, size_t pDecimals );

    // Trimming: Left, Right, and both sides
    DLL_EXPORT_LzServices void LTrim( string & s );
    DLL_EXPORT_LzServices void RTrim( string & s );
    DLL_EXPORT_LzServices void Trim( string & s );
    DLL_EXPORT_LzServices string GetTrimmed( const string & s );

    // String replacement
    DLL_EXPORT_LzServices std::string ReplaceAll( const string & pPattern, const string & pBy, const string & pIn );

    // String splitting
    //
    // If pKeepEmptyParts==true, then only those parts surrounded on *both sides* by separators will be kept.
    //
    DLL_EXPORT_LzServices void SplitString( const string & pStr, const vector<string> & pSeps, vector<string> & pItems, bool pKeepEmptyParts=false );

    // Timestamp formatting
    DLL_EXPORT_LzServices string FormatTimestamp( std::chrono::time_point<std::chrono::system_clock> p_TS );
    DLL_EXPORT_LzServices string FormatTimestampNow();
#pragma endregion


#pragma region "Scope terminator"
    // Capture: http://en.cppreference.com/w/cpp/language/lambda
    class DLL_EXPORT_LzServices Detonator
    {
    public:
        enum class Type { Defusable, NonDefusable };

        Detonator( Type pType=Type::Defusable );
        Detonator( std::function<void()> pDest, Type pType=Type::Defusable );
        Detonator( std::function<void()> pConst, std::function<void()> pDest, Type pType=Type::Defusable );
        ~Detonator();

        void set( std::function<void()> pDest );
        void defuse();

    protected:
        Type mType;
        bool mDetonate{ true };
        std::function<void()> mDest;
    };
#pragma endregion


    // Thread safe perf meter
    DLL_EXPORT_LzServices std::chrono::time_point<std::chrono::steady_clock> Chrono_Now();
    DLL_EXPORT_LzServices double Chrono_MSecs( std::chrono::time_point<std::chrono::steady_clock> pStart );


//---------------------------------------------
// BEGIN: LinkZ Log
//---------------------------------------------
#ifdef USING_LzLog

    // Redirect stdout to a file
    DLL_EXPORT_LzServices void StdoutToFile( const string & pFile );

    // Log file manipulation
    DLL_EXPORT_LzServices const std::thread::id & MainThreadId();
    //
    DLL_EXPORT_LzServices void OpenLogFile( const string & pFile );
    DLL_EXPORT_LzServices void CloseLogFile();
    DLL_EXPORT_LzServices void EnableLog();
    DLL_EXPORT_LzServices void DisableLog();

    // Log file tail
    DLL_EXPORT_LzServices List<string> LogTail();
    DLL_EXPORT_LzServices void SetLogTailLength( unsigned int pLen );

    // Log line
    DLL_EXPORT_LzServices void LogLine( const string & pLine );

#if 1
#define LINKZ_LOG_TXT(TXT)\
{\
    std::stringstream lSStr;                           \
    lSStr << TXT;                                      \
    LzServices::LogLine(lSStr.str());                  \
}
#else
//    A SUPPRIMER APRES REFCTO LOG YAML
//    A SUPPRIMER APRES REFCTO LOG YAML
//    A SUPPRIMER APRES REFCTO LOG YAML
//    A SUPPRIMER APRES REFCTO LOG YAML
//    A SUPPRIMER APRES REFCTO LOG YAML
//    A SUPPRIMER APRES REFCTO LOG YAML
//
#define LINKZ_LOG_TXT(TXT)\
{\
    std::stringstream lSStr;                           \
    lSStr << LzServices::LzLogNodeObject::TabString(); \
    lSStr << TXT;                                      \
    LzServices::LogLine(lSStr.str());                  \
}
#endif

    // Log node: log message with indentation
    class DLL_EXPORT_LzServices LzLogNodeObject
    {
    public:
        LzLogNodeObject( bool pEnabled, bool pTimer, const string & pMsg="", const string & pThreadInfo="" );
        ~LzLogNodeObject();

        static string TabString();
        static string SpaceString();

    private:
#if defined(EXPORT_LzServices) && defined(__MSVC_LIB__)
        static /*thread_local not allowed in MSVC libs*/ unsigned int sTabs;
#else
        static thread_local unsigned int sTabs;
#endif
        bool mEnabled;

        string mMsg;

        bool mTimer;
        string mThreadInfo;
        std::chrono::time_point<std::chrono::steady_clock> mStart;
    };

// Thread info: name + id. Name is automatically set to main if main thread is detected
#define THREAD_ID std::this_thread::get_id()
#define THREAD_INFO(THREAD) "["<<(THREAD_ID == LzServices::MainThreadId() ? "main" : (THREAD))<<" "<<THREAD_ID<<"]"

// Node: function + file + line number
#define LzLogNode(THREAD,MSG)\
    LINKZ_LOG_TXT("/-- "<<THREAD_INFO(THREAD)<<" "<<MSG<<" ["<<__FILE__<<"("<<__LINE__<<"), "<<__FUNCTION__<<"]")\
    LzServices::LzLogNodeObject lLzLogNodObj(true, false);

// Node: function + file + line number + timer
#define LzLogTimeNode(THREAD,MSG)\
    std::stringstream lLzLogNodeObjStrStr; \
    lLzLogNodeObjStrStr << THREAD_INFO(THREAD); \
    const std::string lLzLogNodeObjStr = lLzLogNodeObjStrStr.str(); \
    \
    std::stringstream lLzLogNodeObj_Msg_StrStr; \
    lLzLogNodeObj_Msg_StrStr << MSG; \
    const std::string lLzLogNodeObj_Msg_Str = lLzLogNodeObj_Msg_StrStr.str(); \
    \
    LINKZ_LOG_TXT("/-- "<<lLzLogNodeObjStr<<" "<<lLzLogNodeObj_Msg_Str<<" ["<<__FILE__<<"("<<__LINE__<<"), "<<__FUNCTION__<<"]")\
    LzServices::LzLogNodeObject lLzLogNodObj(true, true, lLzLogNodeObj_Msg_Str, lLzLogNodeObjStr);

// Node: simple
#define LzLogN(THREAD,MSG)\
    LINKZ_LOG_TXT("/-- "<<THREAD_INFO(THREAD)<<" "<<MSG)\
    LzServices::LzLogNodeObject lLzLogNodObj(true, false);

// Node: simple + timer
#define LzLogTimeN(THREAD,MSG)\
    std::stringstream lLzLogNodeObj_Thread_StrStr; \
    lLzLogNodeObj_Thread_StrStr << THREAD_INFO(THREAD); \
    const std::string lLzLogNodeObj_Thread_Str = lLzLogNodeObj_Thread_StrStr.str(); \
    \
    std::stringstream lLzLogNodeObj_Msg_StrStr; \
    lLzLogNodeObj_Msg_StrStr << MSG; \
    const std::string lLzLogNodeObj_Msg_Str = lLzLogNodeObj_Msg_StrStr.str(); \
    \
    LINKZ_LOG_TXT("/-- "<<lLzLogNodeObj_Thread_Str<<" "<<lLzLogNodeObj_Msg_Str)\
    LzServices::LzLogNodeObject lLzLogNodObj(true, true, lLzLogNodeObj_Msg_Str, lLzLogNodeObj_Thread_Str);

// Optional Node: function + file + line number
#define LzOptLogNode(ENABLED,THREAD,MSG)\
    if( ENABLED ) LINKZ_LOG_TXT("/-- "<<THREAD_INFO(THREAD)<<" "<<MSG<<" ["<<__FILE__<<"("<<__LINE__<<"), "<<__FUNCTION__<<"]")\
    LzServices::LzLogNodeObject lLzLogNodObj(ENABLED, false);

// Optional Node: simple
#define LzOptLogN(ENABLED,THREAD,MSG)\
    if( ENABLED ) LINKZ_LOG_TXT("/-- "<<THREAD_INFO(THREAD)<<" "<<MSG)\
    LzServices::LzLogNodeObject lLzLogNodObj(ENABLED, false);

// Optional Node: simple + timer
#define LzOptLogTimeN(ENABLED,THREAD,MSG)\
    if( ENABLED ) LINKZ_LOG_TXT("/-- "<<THREAD_INFO(THREAD)<<" "<<MSG)\
    LzServices::LzLogNodeObject lLzLogNodObj(ENABLED, true);

// Message: function + file + line number
#define LzLogMsg(THREAD,MSG)\
    LINKZ_LOG_TXT("   "<<THREAD_INFO(THREAD)<<" "<<MSG<<" ["<<__FILE__<<"("<<__LINE__<<"), "<<__FUNCTION__<<"]")

// Message: simple
#define LzLogM(THREAD,MSG)\
    LINKZ_LOG_TXT("   "<<THREAD_INFO(THREAD)<<" "<<MSG)

// Error: function + file + line number
#define LzLogErr(THREAD,MSG)\
    LINKZ_LOG_TXT(" # "<<THREAD_INFO(THREAD)<<" "<<MSG<<" ["<<__FILE__<<"("<<__LINE__<<"), "<<__FUNCTION__<<"]")

// Error: simple
#define LzLogE(THREAD,MSG)\
    LINKZ_LOG_TXT(" # "<<THREAD_INFO(THREAD)<<" "<<MSG)

// Exception
#define LzLogException(THREAD, MSG)\
{\
    std::stringstream lStrStream;                           \
    lStrStream <<THREAD_INFO(THREAD)<<": "<<MSG<<" ["<< __FILE__ << "(" << __LINE__ << "), " << __FUNCTION__ << "]";\
    LzLogE(THREAD, "About to throw: " << lStrStream.str()); \
    throw std::runtime_error( lStrStream.str().c_str() );   \
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
#ifdef NO_GUI

// Wait cursor
#define _LzWaitCursor {/** nothing **/}

#define LzLogErrorMessage(TXT, CAP)   {/** nothing **/}
#define LzLogWarningMessage(TXT, CAP) {/** nothing **/}
#define LzLogInfoMessage(TXT, CAP)    {/** nothing **/}

#define LzLogErrorMessageCSS(TXT, CAP, CSS)   {/** nothing **/}
#define LzLogWarningMessageCSS(TXT, CAP, CSS) {/** nothing **/}
#define LzLogInfoMessageCSS(TXT, CAP, CSS)    {/** nothing **/}

#elif defined(USING_QT)

// Wait cursor
class DLL_EXPORT_LzServices LzWaitCursor
{
public:
    LzWaitCursor();
    ~LzWaitCursor();
protected:
    static unsigned int sNested;
};
#define _LzWaitCursor LzServices::LzWaitCursor lLzWaitCursor;


// See QMessageBox width adjustment topics:
// https://stackoverflow.com/questions/37668820/how-can-i-resize-qmessagebox


// Pop-ups without style
#define LzLogErrorMessage(TXT, CAP)\
{\
    std::stringstream lStrStream_TXT, lStrStream_CAP;                  \
    lStrStream_TXT << TXT;                                             \
    lStrStream_CAP << CAP << " ";                                      \
    LzLogE("main", "*** Error message to user: "<<TXT<<" ("<<CAP<<")") \
    QMessageBox::critical(nullptr, lStrStream_CAP.str().c_str(), lStrStream_TXT.str().c_str(), QMessageBox::Ok, QMessageBox::Ok);\
}
#define LzLogWarningMessage(TXT, CAP)\
{\
    std::stringstream lStrStream_TXT, lStrStream_CAP;                    \
    lStrStream_TXT << TXT;                                               \
    lStrStream_CAP << CAP << " ";                                        \
    LzLogE("main", "--- Warning message to user: "<<TXT<<" ("<<CAP<<")") \
    QMessageBox::warning(nullptr, lStrStream_CAP.str().c_str(), lStrStream_TXT.str().c_str(), QMessageBox::Ok, QMessageBox::Ok);\
}
#define LzLogInfoMessage(TXT, CAP)\
{\
    std::stringstream lStrStream_TXT, lStrStream_CAP;                 \
    lStrStream_TXT << TXT;                                            \
    lStrStream_CAP << CAP << " ";                                     \
    LzLogM("main", "... Info message to user: "<<TXT<<" ("<<CAP<<")") \
    QMessageBox::information(nullptr, lStrStream_CAP.str().c_str(), lStrStream_TXT.str().c_str(), QMessageBox::Ok, QMessageBox::Ok);\
}

// Pop-ups WITH style
#define LzLogErrorMessageCSS(TXT, CAP, CSS)\
{\
    std::stringstream lStrStream_TXT, lStrStream_CAP;                  \
    lStrStream_TXT << TXT;                                             \
    lStrStream_CAP << CAP << " ";                                      \
    LzLogE("main", "*** Error message to user: "<<TXT<<" ("<<CAP<<")") \
    QMessageBox lBox;                                                  \
    lBox.setStyleSheet(CSS);                                           \
    lBox.setIcon(QMessageBox::Critical);                               \
    lBox.setWindowTitle(lStrStream_CAP.str().c_str());                 \
    lBox.setText(lStrStream_TXT.str().c_str());                        \
    lBox.exec();                                                       \
}
#define LzLogWarningMessageCSS(TXT, CAP, CSS)\
{\
    std::stringstream lStrStream_TXT, lStrStream_CAP;                    \
    lStrStream_TXT << TXT;                                               \
    lStrStream_CAP << CAP << " ";                                        \
    LzLogE("main", "--- Warning message to user: "<<TXT<<" ("<<CAP<<")") \
    QMessageBox lBox;                                                    \
    lBox.setStyleSheet(CSS);                                             \
    lBox.setIcon(QMessageBox::Warning);                                  \
    lBox.setWindowTitle(lStrStream_CAP.str().c_str());                   \
    lBox.setText(lStrStream_TXT.str().c_str());                          \
    lBox.exec();                                                         \
}
#define LzLogInfoMessageCSS(TXT, CAP, CSS)\
{\
    std::stringstream lStrStream_TXT, lStrStream_CAP;                 \
    lStrStream_TXT << TXT;                                            \
    lStrStream_CAP << CAP << " ";                                     \
    LzLogM("main", "... Info message to user: "<<TXT<<" ("<<CAP<<")") \
    QMessageBox lBox;                                                 \
    lBox.setStyleSheet(CSS);                                          \
    lBox.setIcon(QMessageBox::Information);                           \
    lBox.setWindowTitle(lStrStream_CAP.str().c_str());                \
    lBox.setText(lStrStream_TXT.str().c_str());                       \
    lBox.exec();                                                      \
}

#else
    Some GUI must be defined or use NO_GUI!
#endif
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// Try catch
#define _LzLogTry try {

#define _LzLogFinalCatch(MSG, CAP) }\
    catch( const std::exception & pExc )\
    {\
        LzLogErr("main", "Final catch. "<<MSG<<" ("<<CAP<<")")\
        LzLogErr("main", "Final catch. Exception: "<<pExc.what())\
        LzLogErrorMessage(MSG, CAP)\
    }\
    catch( ... )\
    {\
        LzLogErr("main", "Final catch. "<<MSG<<" ("<<CAP<<")")\
        LzLogErr("main", "Final catch. Unknown exception!")\
        LzLogErrorMessage(MSG, CAP)\
    }

#define _LzLogFinalCatchCSS(MSG, CAP, CSS) }\
    catch( const std::exception & pExc )\
    {\
        LzLogErr("main", "Final catch. "<<MSG<<" ("<<CAP<<")")\
        LzLogErr("main", "Final catch. Exception: "<<pExc.what())\
        LzLogErrorMessageCSS(MSG, CAP, CSS)\
    }\
    catch( ... )\
    {\
        LzLogErr("main", "Final catch. "<<MSG<<" ("<<CAP<<")")\
        LzLogErr("main", "Final catch. Unknown exception!")\
        LzLogErrorMessageCSS(MSG, CAP, CSS)\
    }

#define _LzLogCatchAndThrow(TXT) }\
    catch( const std::exception & pExc )\
    {\
        LzLogErr("", "Catch and throw. Caught exception: "<<pExc.what())\
        /*LzLogErr("", "Stack trace:\r\n"+pExc.);*/\
        LzLogException("", TXT)\
    }\
    catch( ... )\
    {\
        LzLogErr("", "Catch and throw. Caught unknown exception.")\
        LzLogException("", TXT)\
    }
#endif
//---------------------------------------------
// END: LinkZ Log
//---------------------------------------------
}


//---------------------------------------------
// Third party Log
//---------------------------------------------
#ifndef USING_LzLog

    //// Node: function + file + line number
    //#define LzLogNode(THREAD,MSG) xxxx
	//
    //// Node: simple
    //#define LzLogN(THREAD,MSG) xxxx
	//
    //// Optional Node: function + file + line number
    //#define LzOptLogNode(ENABLED, THREAD,MSG) xxxx
    //
    //// Optional Node: simple
    //#define LzOptLogN(ENABLED, THREAD,MSG) xxxx
    //
    //// Message: function + file + line number
    //#define LzLogMsg(THREAD,MSG) xxxx
	//
    //// Message: simple
    //#define LzLogM(THREAD,MSG) xxxx
	//
    //// Error: function + file + line number
    //#define LzLogErr(THREAD,MSG) xxxx
	//
    //// Error: simple
    //#define LzLogE(THREAD,MSG) xxxx
	//
	//// Exception
    //#define LzLogException(THREAD, MSG) xxxx

// Try catch
//    #define _LzLogTry try {

//    #define _LzLogFinalCatch(MSG, CAP) }\
//        catch( const std::exception & pExc )\
//        {\
//            LzLogErr("main", "Final catch. "<<MSG<<" ("<<CAP<<")")\
//            LzLogErr("main", "Final catch. Exception: "<<pExc.what())\
//            LzLogErrorMessage(MSG, CAP)\
//        }

//    #define _LzLogFinalCatchCSS(MSG, CAP, CSS) }\
//        catch( const std::exception & pExc )\
//        {\
//            LzLogErr("main", "Final catch. "<<MSG<<" ("<<CAP<<")")\
//            LzLogErr("main", "Final catch. Exception: "<<pExc.what())\
//            LzLogErrorMessageCSS(MSG, CAP, CSS)\
//        }



#include <ZglurbServices/Log.h>

// Node: function + file + line number
#define LzLogNode(THREAD,MSG) ZglurbLogNode("[" << THREAD << " " << THREAD_ID << "] " << MSG)

// Node: simple
#define LzLogN(THREAD,MSG) ZglurbLogN("[" << THREAD << " " << THREAD_ID << "] " << MSG)

// Message: function + file + line number
#define LzLogMsg(THREAD,MSG) ZglurbLogMsg("[" << THREAD << " " << THREAD_ID << "] " << MSG)

// Message: simple
#define LzLogM(THREAD,MSG) ZglurbLogM("[" << THREAD << " " << THREAD_ID << "] " << MSG)

// Error: function + file + line number
#define LzLogErr(THREAD,MSG) ZglurbLogErr("[" << THREAD << " " << THREAD_ID << "] " << MSG)

// Error: simple
#define LzLogE(THREAD,MSG) ZglurbLogE("[" << THREAD << " " << THREAD_ID << "] " << MSG)

// Exception
#define LzLogException(THREAD, MSG) ZglurbLogException("[" << THREAD << " " << THREAD_ID << "] " << MSG)

#endif
