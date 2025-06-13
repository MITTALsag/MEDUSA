#include "PolyMatrix.h"
#include "Matrix.h"


namespace LzMath
{
//==============================================================================
PolyMatrix::PolyMatrix()
{
	m_H = m_W = 0;
    m_pMat = nullptr;
}

//==============================================================================
PolyMatrix::~PolyMatrix()
{
	Free();
}

//==============================================================================
void PolyMatrix::Free()
{
    // Non empty matrix?
    if( m_pMat )
	{
		// Free
		for( int i=0 ; i<m_H ; i++ )
			delete [] m_pMat[i];

		delete [] m_pMat;
        m_pMat = nullptr;
	}

	// Mark as freed
	m_H = m_W = 0;
}

//==============================================================================
void PolyMatrix::SetDims( unsigned int p_H, unsigned int p_W )
{
	// Same dims ?
	if( p_H==m_H && p_W==m_W )
        return;

	// Free previous
	Free();

    // Non empty matrix?
    if( p_H*p_W != 0 )
    {
        // Alloc
        m_pMat = new Polynom * [p_H];
        for( int i=0 ; i<p_H ; i++ )
            m_pMat[i] = new Polynom [p_W];

        // Store new dimensions
        m_H = p_H;
        m_W = p_W;
    }
}

//==============================================================================
Polynom & PolyMatrix::Elt( unsigned int p_I, unsigned int p_J )
{
    // Check
    if( p_I>=m_H || p_J>=m_W )
        LzLogException("", "[PolyMatrix::Elt] : Accessing invalid element (i="<<p_I<<",j="<<p_J<<")!")

	return m_pMat[p_I][p_J];
}

//==============================================================================
Polynom PolyMatrix::Elt( unsigned int p_I, unsigned int p_J ) const
{
    // Check
    if( p_I>=m_H || p_J>=m_W )
        LzLogException("", "[PolyMatrix::Elt] : Accessing invalid element (i="<<p_I<<",j="<<p_J<<") !")

	return m_pMat[p_I][p_J];
}

//==============================================================================
void PolyMatrix::SubMatrixTo( unsigned int p_FirstRow, unsigned int p_LastRow,
                              unsigned int p_FirstCol, unsigned int p_LastCol, PolyMatrix & p_To )
{
    // Check
    if( p_FirstRow>p_LastRow || p_FirstCol>p_LastCol
     || p_LastRow>=m_H || p_LastCol>=m_W )
	{
        LzLogException("", "[PolyMatrix::SubMatrixTo] : Invalid range!")
	}

	// Fill new matrix
	p_To.SetDims( p_LastRow-p_FirstRow+1, p_LastCol-p_FirstCol+1 );
    for( unsigned int i=0 ; i<p_To.Rows() ; i++ )
    for( unsigned int j=0 ; j<p_To.Cols() ; j++ )
	{
		p_To.Elt(i,j) = Elt( p_FirstRow+i, p_FirstCol+j );
	}
}

////////////////////////////////////////////////////////////////////////////////
//
// Operations
//
//==============================================================================
void PolyMatrix::Mult( const Matrix & p_A, PolyMatrix & p_Res )
{
    // Check
    if( !Rows() || !Cols() || p_A.Rows()!=Cols() || !p_A.Cols() )
        LzLogException("", "[PolyMatrix::Mult] : Unable to multiplicate!")

    // Set dimensions
    p_Res.SetDims( Rows(), p_A.Cols() );

    // Mult
    for( unsigned int i=0 ; i<p_Res.Rows() ; i++ )
    for( unsigned int j=0 ; j<p_Res.Cols() ; j++ )
	{
        for( unsigned int k=0 ; k<Cols() ; k++ )
			p_Res.Elt(i,j) = p_Res.Elt(i,j) + Elt(i,k) * p_A.Elt(k,j);
	}
}
	
//==============================================================================
void PolyMatrix::Mult( const PolyMatrix & p_A, PolyMatrix & p_Res )
{
    // Check
	if( !Rows() || !Cols() || p_A.Rows()!=Cols() || !p_A.Cols() )
        LzLogException("", "[PolyMatrix::Mult] : Unable to multiplicate!")

    // Set dimensions
    p_Res.SetDims( Rows(), p_A.Cols() );

    // Mult
    for( unsigned int i=0 ; i<p_Res.Rows() ; i++ )
    for( unsigned int j=0 ; j<p_Res.Cols() ; j++ )
	{
        for( unsigned int k=0 ; k<Cols() ; k++ )
			p_Res.Elt(i,j) = p_Res.Elt(i,j) + Elt(i,k) * p_A.Elt(k,j);
	}
}

//==============================================================================
void PolyMatrix::TransposeTo( PolyMatrix & p_Res )
{
    // Set dimensions
    p_Res.SetDims(Cols(), Rows());

    // Transpose
    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
	{
		p_Res.Elt(j,i) = Elt(i,j);
	}
}



////////////////////////////////////////////////////////////////////////////////
//
// Evaluation
//
//==============================================================================
void PolyMatrix::EvalTo( double X, double Y, Matrix & p_To ) const
{
    // Set dimensions
    p_To.SetDims(m_H, m_W);

    // Eval to
    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
        p_To.Elt(i,j) = Elt(i,j).EvalAt(X,Y);
}

//==============================================================================
void PolyMatrix::EvalTo( double X, double Y, double Z, Matrix & p_To ) const
{
    // Set dimensions
    p_To.SetDims(Rows(), Cols());

    // Eval to
    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
        p_To.Elt(i,j) = Elt(i,j).EvalAt(X,Y,Z);
}


////////////////////////////////////////////////////////////////////////////////
//
// Integration
//
//==============================================================================
void PolyMatrix::IntegralTo( unsigned int p_MaxDegX, unsigned int p_MaxDegY, unsigned int p_MaxDegZ, double * p_pSums, Matrix & p_Res )
{
    // Set dimensions
    p_Res.SetDims(m_H, m_W);

	// Compute integration matrix
    for( unsigned int i=0 ; i<m_H ; i++ )
    for( unsigned int j=0 ; j<m_W ; j++ )
	{
        // Integrate polynom
		double SumPol;
        Elt(i,j).IntegralTo(p_MaxDegX, p_MaxDegY, p_MaxDegZ, p_pSums, SumPol);

        // Commit
        p_Res.Elt(i,j) = SumPol;
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// Debug
//
//==============================================================================
void PolyMatrix::Log( string p_Prefix/*=""*/ )
{
    // Log
    LzLogN("", "PolyMatrix "<<p_Prefix)

	if( !m_H )
	{
        LzLogM("", "Empty matrix")
	}
	else
	{
        for( unsigned int i=0 ; i<m_H ; i++ )
		{
            string l_Row = p_Prefix+" ["+std::to_string(i)+"] = [\t";

            for( unsigned int j=0 ; j<m_W ; j++ )
			{
				l_Row += "\t\t" + Elt(i,j).ToString();
			}
			l_Row += "\t\t]";

            LzLogM("", l_Row)
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// Co-Facteurs
//
//==============================================================================
thread_local List<int> PolyMatrix::m_LocalRows;
thread_local List<int> PolyMatrix::m_LocalCols;

//==============================================================================
Polynom PolyMatrix::CF_Det( bool p_EmptyMat/*=false*/ )
{
	// Ok ?
	if( !m_H || m_H!=m_W )
        LzLogException("", "[Matrix::CF_Det] : Unable to compute determinant!")

	// Init local cols and rows
	InitLocalRowsCols();

	return CF_RecDet(p_EmptyMat);
}

//==============================================================================
int PolyMatrix::MaxZeroLocalRow()
{
	int l_MaxZeroCount = -1;
	int l_LocalRow = 0;
	int l_BestRow;

    BrowseList( i_RowPos, m_LocalRows )
	{
		int l_GlobalRow = m_LocalRows.GetAt(i_RowPos);
		int l_ZeroCount = 0;

        BrowseList( i_ColPos, m_LocalCols )
		{
			int l_GlobalCol = m_LocalCols.GetAt(i_ColPos);

			if( Elt(l_GlobalRow,l_GlobalCol).IsZero() )
				l_ZeroCount++;
		}

		if( l_ZeroCount > l_MaxZeroCount )
		{
			l_MaxZeroCount = l_ZeroCount;
			l_BestRow = l_LocalRow;
		}

		l_LocalRow++;
	}

	return l_BestRow;
}

//==============================================================================
Polynom PolyMatrix::CF_RecDet( bool p_EmptyMat/*=false*/ )
{
	// 0x0 ?
	if( LocalCols() == 0 )
	{
		return 1;
	}
	// 1x1 ?
	else if( LocalCols() == 1 )
	{
		return Elt(m_LocalRows.GetHead(),m_LocalCols.GetHead());
	}
	// 2x2 ?
	else if( LocalCols() == 2 )
	{
        Polynom l_00 = Elt(m_LocalRows.GetHead(),m_LocalCols.GetHead());
        Polynom l_11 = Elt(m_LocalRows.GetTail(),m_LocalCols.GetTail());
        Polynom l_01 = Elt(m_LocalRows.GetHead(),m_LocalCols.GetTail());
        Polynom l_10 = Elt(m_LocalRows.GetTail(),m_LocalCols.GetHead());

		return l_00*l_11 - l_01*l_10;
	}

	//--- Decrease dim and recurse

	// Remove first local row
	int l_LocalI;
	if( p_EmptyMat )
		l_LocalI = MaxZeroLocalRow();
	else
		l_LocalI = 0;

	int l_GlobalRow = RemoveLocalRow(l_LocalI);

	// Develop along first row
    Polynom l_Det = 0;
	for( int i_LocalJ=0 ; i_LocalJ<LocalCols() ; i_LocalJ++ )
	{
		// Remove local col
		int l_GlobalCol = RemoveLocalCol(i_LocalJ);

		// Compute sub det only if element non null
        Polynom l_Elt = Elt(l_GlobalRow,l_GlobalCol);
		if( !l_Elt.IsZero() )
			l_Det = l_Det + ((l_LocalI+i_LocalJ)%2 ? -1 : +1) * l_Elt * CF_RecDet(p_EmptyMat);

		// Add global col for further computing
		AddGlobalCol( l_GlobalCol );
	}

	// Add global row further computing
	AddGlobalRow( l_GlobalRow );

	return l_Det;
}

//==============================================================================
void PolyMatrix::InitLocalRowsCols()
{
	m_LocalRows.DelAll();
	m_LocalCols.DelAll();
	for( int i=0 ; i<m_H ; i++ )
	{
		m_LocalRows.AddTail(i);
		m_LocalCols.AddTail(i);
	}
}
unsigned int PolyMatrix::LocalCols()
{
    return m_LocalCols.Count();
}
int PolyMatrix::RemoveLocalRow( int p_LocalI )
{
	if( p_LocalI<0 || p_LocalI>=m_LocalRows.Count() )
        LzLogException("", "[PolyMatrix::RemoveLocalRow] : Could not remove "<<p_LocalI<<"!")

    void * l_Pos = m_LocalRows.FindIdx(p_LocalI);
	int l_GlobalRow = m_LocalRows.GetAt( l_Pos );
	m_LocalRows.DelAt( l_Pos );
	return l_GlobalRow;
}
int PolyMatrix::RemoveLocalCol( int p_LocalJ )
{
	if( p_LocalJ<0 || p_LocalJ>=m_LocalCols.Count() )
        LzLogException("", "[PolyMatrix::RemoveLocalCol] : Could not remove "<<p_LocalJ<<"!")

    void * l_Pos = m_LocalCols.FindIdx(p_LocalJ);
	int l_GlobalCol = m_LocalCols.GetAt( l_Pos );
	m_LocalCols.DelAt( l_Pos );
	return l_GlobalCol;
}
void PolyMatrix::AddGlobalRow( int p_GlobalI )
{
	void * i_Pos = m_LocalRows.HeadPos();
	while( i_Pos && m_LocalRows.GetAt(i_Pos)<p_GlobalI )
		m_LocalRows.Next(i_Pos);
	if( i_Pos )
		m_LocalRows.AddBefore( i_Pos, p_GlobalI );
	else
		m_LocalRows.AddTail(p_GlobalI);
}
void PolyMatrix::AddGlobalCol( int p_GlobalJ )
{
	void * i_Pos = m_LocalCols.HeadPos();
	while( i_Pos && m_LocalCols.GetAt(i_Pos)<p_GlobalJ )
		m_LocalCols.Next(i_Pos);
	if( i_Pos )
		m_LocalCols.AddBefore( i_Pos, p_GlobalJ );
	else
		m_LocalCols.AddTail(p_GlobalJ);
}

//==============================================================================
void PolyMatrix::CF_InverseTo( PolyMatrix & p_Cof, Polynom & p_Det )
{
	// Ok ?
	if( !Rows() || Rows()!=Cols() )
        LzLogException("", "[PolyMatrix::CF_InverseTo] : Unable to inverse!")

	// Inversible ?
	p_Det = CF_Det();
	if( p_Det.IsZero() )
        LzLogException("", "[PolyMatrix::CF_InverseTo] : Matrix not inversible!")

	// Dims
	p_Cof.SetDims(Rows(),Cols());

	// Fill with transpose(co-factors)
    for( unsigned int i=0 ; i<Rows() ; i++ )
    for( unsigned int j=0 ; j<Cols() ; j++ )
	{
		// Remove
		RemoveLocalRow(i);
		RemoveLocalCol(j);

		// Compute co-factor and transpose
		p_Cof.Elt(j,i) = ((i+j)%2 ? -1 : +1) * CF_RecDet();

		// Add
		AddGlobalRow(i);
		AddGlobalCol(j);
	}
}

}
