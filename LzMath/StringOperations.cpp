#include "StringOperations.h"
#include <LzServices/LzLog.h>

namespace LzMath
{
namespace ToolBox
{
#ifdef USE_CLI
#pragma region "String ops"
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^StrToInt32( String ^pStr, System::Int32 & pInt, System::Int32 pClampMin/*=+1*/, System::Int32 pClampMax/*=-1*/ )
{
    try
    {
        if( pStr == "+" || pStr == "-" )
            pInt = 0;
        else
            pInt = System::Int32::Parse( pStr );

        // Clamp?
        if( pClampMin <= pClampMax )
        {
            if( pInt < pClampMin )
            {
                pInt = pClampMin;
                return "" + pInt;
            }

            if( pInt > pClampMax )
            {
                pInt = pClampMax;
                return "" + pInt;
            }
        }

        return pStr; // to allow +/-
    }
    catch( ... )
    {
        TxLogErr( "String is not in a numerical format: " << TxServices::StdStrFromCLI(pStr) << "!" );

        return StrToInt32( "0", pInt, pClampMin, pClampMax );
    }
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^StrToUInt32( String ^pStr, System::UInt32 & pInt, System::UInt32 pClampMin/*=+1*/, System::UInt32 pClampMax/*=0*/ )
{
    try
    {
        pInt = System::UInt32::Parse( pStr );

        // Clamp?
        if( pClampMin <= pClampMax )
        {
            if( pInt < pClampMin )
                pInt = pClampMin;

            if( pInt > pClampMax )
                pInt = pClampMax;
        }

        return "" + pInt;
    }
    catch( ... )
    {
        TxLogErr( "String is not in a numerical format: " << TxServices::StdStrFromCLI(pStr) << "!" );

        return StrToUInt32( "0", pInt, pClampMin, pClampMax );
    }
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^StrToDouble( String ^pStr, System::Double & pDouble, System::Double pClampMin/*=+1*/, System::Double pClampMax/*=-1*/ )
{
    try
    {
        // Do not interrupt user when starts typing...
        if( pStr == "+" || pStr == "-" || pStr == "+." || pStr == "-." || pStr == "." )
            pDouble = 0;
        else
            pDouble = System::Double::Parse( pStr, System::Globalization::NumberStyles::AllowLeadingSign
                                             | System::Globalization::NumberStyles::AllowDecimalPoint
                                             | System::Globalization::NumberStyles::AllowExponent
                                             | System::Globalization::NumberStyles::AllowLeadingWhite
                                             | System::Globalization::NumberStyles::AllowTrailingWhite );
        // Clamp?
        if( pClampMin <= pClampMax )
        {
            if( pDouble < pClampMin )
            {
                pDouble = pClampMin;
                return "" + pDouble;
            }

            if( pDouble > pClampMax )
            {
                pDouble = pClampMax;
                return "" + pDouble;
            }
        }

        return pStr;
    }
    catch( ... )
    {
        TxLogErr( "String is not in a numerical format: " << TxServices::StdStrFromCLI(pStr) << "!" );

        return StrToDouble( "0", pDouble, pClampMin, pClampMax );
    }
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^StrToByte( String ^pStr, System::Byte & pByte, System::Byte pClampMin/*=0xFF*/, System::Byte pClampMax/*=0x00*/ )
{
    try
    {
        pByte = System::Byte::Parse( pStr );

        // Clamp?
        if( pClampMin <= pClampMax )
        {
            if( pByte < pClampMin )
                pByte = pClampMin;

            if( pByte > pClampMax )
                pByte = pClampMax;
        }

        return "" + pByte;
    }
    catch( ... )
    {
        TxLogErr( "String is not in a numerical format: " << TxServices::StdStrFromCLI(pStr) << "!" );

        return StrToByte( "0", pByte, pClampMin, pClampMax );
    }
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//std::string StrToStdStr( String ^s )
//{
//    std::string lStr;
//    const char* chars = ( const char* )( System::Runtime::InteropServices::Marshal::StringToHGlobalAnsi( s ) ).ToPointer();
//    lStr = chars;
//    System::Runtime::InteropServices::Marshal::FreeHGlobal( System::IntPtr( ( void* )chars ) );
//    return lStr;
//}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
std::wstring StrToStdUnicodeStr( String ^ws )
{
    std::wstring os;
    const wchar_t* chars = ( const wchar_t* )( System::Runtime::InteropServices::Marshal::StringToHGlobalUni( ws ) ).ToPointer();
    os = chars;
    System::Runtime::InteropServices::Marshal::FreeHGlobal( System::IntPtr( ( void* )chars ) );
    return os;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^StdStrToStr( const std::string & s )
{
    return gcnew String( s.c_str() );
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^StdUnicodeStrToStr( const wstring & ws )
{
    return gcnew String( ws.c_str() );
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
String ^UInt32ToBinaryStr( unsigned int pI )
{
    String ^lRes = "";

    do
    {
        if( pI & 1 )
            lRes = "1" + lRes;
        else
            lRes = "0" + lRes;

        pI >>= 1;
    }
    while( pI );

    return lRes;
}
#pragma endregion

#endif
}
}
