#pragma once

#include "Polynom.h"
#include <LzServices/List.h>


namespace LzMath
{
using LzServices::List;
class Matrix;


class PolyMatrix
{
public:
	// Construction / destruction
    PolyMatrix();
    ~PolyMatrix();

	// Sizing
    void SetDims( unsigned int p_H/*ROWS*/, unsigned int p_W/*COLS*/ );
    unsigned int Rows() const { return m_H; }
    unsigned int Cols() const { return m_W; }
	void Free();

	// Access
    Polynom & Elt( unsigned int p_I/*ROW*/, unsigned int p_J/*COL*/ );
    Polynom Elt( unsigned int p_I/*ROW*/, unsigned int p_J/*COL*/ ) const;
    void SubMatrixTo( unsigned int p_FirstRow, unsigned int p_LastRow,
                      unsigned int p_FirstCol, unsigned int p_LastCol, PolyMatrix & p_To );

	// Operations
    void Mult( const Matrix & p_A, PolyMatrix & p_Res );
    void Mult( const PolyMatrix & p_A, PolyMatrix & p_Res );
    void TransposeTo( PolyMatrix & p_Res );

	// Evaluation
    void EvalTo( double X, double Y, Matrix & p_To ) const;
    void EvalTo( double X, double Y, double Z, Matrix & p_To ) const;

	// Integration
    void IntegralTo( unsigned int p_MaxDegX, unsigned int p_MaxDegY, unsigned int p_MaxDegZ,
                     double * p_pSums, Matrix & p_Res );

	// Debug
    void Log( string p_Prefix="" );

	// Co-Facteurs
    Polynom CF_Det( bool p_EmptyMat=false );
    void CF_InverseTo( PolyMatrix & p_Cof, Polynom & p_Det );


protected:
	// Data
    unsigned int m_H;
    unsigned int m_W;
    Polynom ** m_pMat;

	// Co-Facteurs
	int MaxZeroLocalRow();
    Polynom CF_RecDet( bool p_EmptyMat=false );
	void InitLocalRowsCols();
    static thread_local List<int> m_LocalRows;
    static thread_local List<int> m_LocalCols;
    unsigned int LocalCols();
	int RemoveLocalRow( int p_LocalI );
	int RemoveLocalCol( int p_LocalJ );
	void AddGlobalRow( int p_GlobalI );
	void AddGlobalCol( int p_GlobalJ );
};

}
