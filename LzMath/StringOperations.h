#pragma once

//#include <vector>
#include <string>


namespace LzMath {
namespace ToolBox
{
//using std::vector;
using std::string;
//using std::wstring;


#pragma region "String ops"
	// Computes char length of an integer or double. "-123" ==> 4, "0" ==> 1, "(+)123" ==> 3
    template <typename T> size_t StringLength( T pT )
	{
        return std::to_string(pT).length();
	}

#ifdef USE_CLI
	// Convert System::String to numerical value
	//***** find a way to convert "live" (pb with "-" and ".")
	String ^StrToInt32( String ^pStr, System::Int32 & pInt, System::Int32 pClampMin = +1, System::Int32 pClampMax = -1 );
	String ^StrToUInt32( String ^pStr, System::UInt32 & pInt, System::UInt32 pClampMin = +1, System::UInt32 pClampMax = 0 );
	String ^StrToDouble( String ^pStr, System::Double & pDouble, System::Double pClampMin = +1, System::Double pClampMax = -1 );
	String ^StrToByte( String ^pStr, System::Byte & pByte, System::Byte pClampMin = 0xFF, System::Byte pClampMax = 0x00 );

	// Convert System::String to std::string and std::wstring
	// More info: http://www.codeproject.com/Articles/10400/StringConvertor-A-convertor-class-for-managed-unma
//string StrToStdStr( String ^ws );
	wstring StrToStdUnicodeStr( String ^s );
	String ^StdStrToStr( const string & s );
	String ^StdUnicodeStrToStr( const wstring & ws );

	// Convert an integer to binary string
	String ^UInt32ToBinaryStr( unsigned int pI );
#endif
#pragma endregion
}
}
