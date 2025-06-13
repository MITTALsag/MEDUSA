#include "Polynom.h"
#include <LzGeom/Point3D.h>
#include <LzServices/LzLog.h>
//#include <float.h>
//#include <math.h>
//#include <Math/Vector.h>


namespace LzMath
{
//==============================================================================
static inline int MAX( int A, int B ) { return A > B ? A : B ; }

// What is zero and what is not ....
const double Polynom::TOLERANCE = 1e-15;

//==============================================================================
const Polynom Polynom::XYZ( 1, 1, 1, 1 );
const Polynom Polynom::XY ( 1, 1, 0, 1 );
const Polynom Polynom::XZ ( 1, 0, 1, 1 );
const Polynom Polynom::YZ ( 0, 1, 1, 1 );
const Polynom Polynom::X  ( 1, 0, 0, 1 );
const Polynom Polynom::Y  ( 0, 1, 0, 1 );
const Polynom Polynom::Z  ( 0, 0, 1, 1 );

//==============================================================================
Polynom::Polynom()
{
	// Init polynom
    m_pCoefs = nullptr;

	// Default polynom is P(x,y,z) = 0
	SetDegs(0,0,0);
}

//==============================================================================
Polynom::Polynom( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ, double p_Coef )
{
	// Init polynom
    m_pCoefs = nullptr;
	SetDegs( p_DegX, p_DegY, p_DegZ );
	Monom( p_DegX, p_DegY, p_DegZ) = p_Coef;
}

//==============================================================================
Polynom::Polynom( const Polynom & p_Poly )
{
	// Init polynom
    m_pCoefs = nullptr;
	*this = p_Poly;
}

//==============================================================================
Polynom::Polynom( double p_Monom0 )
{
	// Init polynom
    m_pCoefs = nullptr;
	*this = p_Monom0;
}

//==============================================================================
Polynom::~Polynom()
{
	Free();
}

//==============================================================================
void Polynom::Free()
{
	if( m_pCoefs )
	{
		delete m_pCoefs;
        m_pCoefs = nullptr;
		m_Monoms = 0;
	}
}

//==============================================================================
void Polynom::SetDegs( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ )
{
	Free();

	// Init polynom
	m_DegX = p_DegX;
	m_DegY = p_DegY;
	m_DegZ = p_DegZ;

	m_Monoms = (m_DegX+1)*(m_DegY+1)*(m_DegZ+1);
	m_pCoefs = new double [ m_Monoms ];

	// Reset polynom
	for( int i=0 ; i<m_Monoms ; i++ )
	{
        m_pCoefs[i] = 0.0;
	}
}

//==============================================================================
double & Polynom::Monom( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ )
{
    // Check
    if( p_DegX>m_DegX || p_DegY>m_DegY || p_DegZ>m_DegZ )
        LzLogException("", "[Polynom::Monom] : Degree overflow ! My degs=("<<m_DegX<<","<<m_DegY<<","<<m_DegZ<<")")

	// Get monom
	return m_pCoefs[ p_DegX + (m_DegX+1)*(p_DegY+(m_DegY+1)*p_DegZ) ];
}

//==============================================================================
const double & Polynom::Monom( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ ) const
{
    return ((Polynom*)this)->Monom(p_DegX,p_DegY,p_DegZ);
}


////////////////////////////////////////////////////////////////////////////////
//
// Operators
//
//==============================================================================
Polynom Polynom::operator*( double p_Scal ) const
{
    Polynom l_Res;
	
	l_Res.SetDegs( m_DegX, m_DegY, m_DegZ );

	for( int i=0 ; i<m_Monoms ; i++ )
	{
		l_Res.m_pCoefs[i] = p_Scal * m_pCoefs[i];
	}

	return l_Res;
}

//==============================================================================
Polynom Polynom::operator/( double p_Scal ) const
{
	// Division by 0 ?
    if( std::abs(p_Scal) < TOLERANCE )
        LzLogException("", "[Polynom::operator/] : Division by 0!")
	
    Polynom l_Res;

    l_Res.SetDegs( m_DegX, m_DegY, m_DegZ );

	for( int i=0 ; i<m_Monoms ; i++ )
		l_Res.m_pCoefs[i] = m_pCoefs[i] / p_Scal;

	return l_Res;
}

//==============================================================================
Polynom Polynom::operator*( const Polynom & p_Poly ) const
{
    Polynom l_Res;

	l_Res.SetDegs( m_DegX+p_Poly.m_DegX, m_DegY+p_Poly.m_DegY, m_DegZ+p_Poly.m_DegZ );

	for( int x1=0 ; x1<=m_DegX ; x1++ )
	for( int y1=0 ; y1<=m_DegY ; y1++ )
	for( int z1=0 ; z1<=m_DegZ ; z1++ )
	{
		for( int x2=0 ; x2<=p_Poly.m_DegX ; x2++ )
		for( int y2=0 ; y2<=p_Poly.m_DegY ; y2++ )
		for( int z2=0 ; z2<=p_Poly.m_DegZ ; z2++ )
		{
			l_Res.Monom(x1+x2,y1+y2,z1+z2) += Monom(x1,y1,z1)*p_Poly.Monom(x2,y2,z2);
		}
	}

	return l_Res;
}

//==============================================================================
Polynom Polynom::operator+( double p_Scal ) const
{
    Polynom l_Res;

	l_Res = *this;
	l_Res.m_pCoefs[0] += p_Scal;

	return l_Res;
}

//==============================================================================
Polynom Polynom::operator-( double p_Scal ) const
{
    Polynom l_Res;

	l_Res = *this;
	l_Res.m_pCoefs[0] -= p_Scal;

	return l_Res;
}

//==============================================================================
Polynom Polynom::operator+( const Polynom & p_Poly ) const
{
    Polynom l_Res;

	l_Res.SetDegs( MAX(m_DegX,p_Poly.m_DegX), MAX(m_DegY,p_Poly.m_DegY), MAX(m_DegZ,p_Poly.m_DegZ) );
		
	for( int x=0 ; x<l_Res.m_DegX+1 ; x++ )
	for( int y=0 ; y<l_Res.m_DegY+1 ; y++ )
	for( int z=0 ; z<l_Res.m_DegZ+1 ; z++ )
	{
		if( x<m_DegX+1 && y<m_DegY+1 && z<m_DegZ+1 )
			l_Res.Monom(x,y,z) += Monom(x,y,z);

		if( x<p_Poly.m_DegX+1 && y<p_Poly.m_DegY+1 && z<p_Poly.m_DegZ+1 )
			l_Res.Monom(x,y,z) += p_Poly.Monom(x,y,z);
	}

	return l_Res;
}

//==============================================================================
Polynom Polynom::operator-( const Polynom & p_Poly ) const
{
    Polynom l_Res;

	l_Res.SetDegs( MAX(m_DegX,p_Poly.m_DegX), MAX(m_DegY,p_Poly.m_DegY), MAX(m_DegZ,p_Poly.m_DegZ) );

    for( unsigned int x=0 ; x<l_Res.m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<l_Res.m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<l_Res.m_DegZ+1 ; z++ )
	{
		if( x<m_DegX+1 && y<m_DegY+1 && z<m_DegZ+1 )
			l_Res.Monom(x,y,z) += Monom(x,y,z);

		if( x<p_Poly.m_DegX+1 && y<p_Poly.m_DegY+1 && z<p_Poly.m_DegZ+1 )
			l_Res.Monom(x,y,z) -= p_Poly.Monom(x,y,z);
	}

	return l_Res;
}

//==============================================================================
//Polynom & Polynom::operator+=( const Polynom & p_Poly )
//{
//	*this = *this + p_Poly;
//	return *this;
//}

//==============================================================================
Polynom & Polynom::operator=( const Polynom & p_Poly )
{
	SetDegs( p_Poly.m_DegX, p_Poly.m_DegY, p_Poly.m_DegZ );

    for( unsigned int i=0 ; i<m_Monoms ; i++ )
	{
		m_pCoefs[i] = p_Poly.m_pCoefs[i];
	}

	return *this;
}

//==============================================================================
Polynom & Polynom::operator=( double p_Monom0 )
{
	SetDegs( 0, 0, 0 );
	m_pCoefs[0] = p_Monom0;
	return *this;
}

//==============================================================================
static const Polynom One = 1;
Polynom Polynom::PolPow( int Pow ) const
{
	if( Pow )
		return (*this) * PolPow(Pow-1);
	else
		return One;
}


////////////////////////////////////////////////////////////////////////////////
//
// Evaluation
//
//==============================================================================
// pow( double x, double y ) as defined in <math.h> is much slower for integer powers
//static double IntPow( double X, int N )
//{
//	double lRes = 1;
//	while( N )
//	{
//		lRes *= X;
//		N--;
//	}

//	return lRes;
//}

//==============================================================================
double Polynom::EvalAt( double X ) const
{
	double Res = 0;

    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
	{
        Res += Monom(x,0,0) * std::pow(X, x); //IntPow(X,x);
	}

	return Res;
}

//==============================================================================
double Polynom::EvalAt( double X, double Y ) const
{
	double Res = 0;

    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
	{
        Res += Monom(x,y,0) * std::pow(X, x) * std::pow(Y, y); //IntPow(X,x) * IntPow(Y,y);
	}

	return Res;
}

//==============================================================================
double Polynom::EvalAt( double X, double Y, double Z ) const
{
	double Res = 0;

    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
	{
        Res += Monom(x,y,z) * std::pow(X, x) * std::pow(Y, y) * std::pow(Z, z); //IntPow(X,x) * IntPow(Y,y) * IntPow(Z,z);
	}

	return Res;
}

//==============================================================================
double Polynom::EvalAt( const Point3D & pPt3D ) const
{
    return EvalAt( pPt3D.X(), pPt3D.Y(), pPt3D.Z() );
}

//==============================================================================
bool Polynom::IsZero( double tol/*=TOLERANCE*/ ) const
{
	for( int x=0 ; x<m_DegX+1 ; x++ )
	for( int y=0 ; y<m_DegY+1 ; y++ )
	for( int z=0 ; z<m_DegZ+1 ; z++ )
	{
        if( std::abs(Monom(x,y,z)) >= tol )
            return false;
	}

    return true;
}


////////////////////////////////////////////////////////////////////////////////
//
// Composition
//
//==============================================================================
Polynom Polynom::Compose( const Polynom & PolX, const Polynom & PolY, const Polynom & PolZ ) const
{
    Polynom Res;

	// Sum for all monoms
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
	{
		// If monom non null
		double Coef = Monom(x,y,z);

		if( fabs(Coef) > TOLERANCE )
		{
			Res = Res + Coef * PolX.PolPow(x) * PolY.PolPow(y) * PolZ.PolPow(z);
		}
	}

	return Res;
}

//==============================================================================
Polynom Polynom::Compose( const Polynom & PolX, const Polynom & PolY ) const
{
    Polynom Res;

	// Sum for all monoms
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
	{
		// If monom non null
		double Coef = Monom(x,y,0);

        if( std::abs(Coef) > TOLERANCE )
		{
			Res = Res + Coef * PolX.PolPow(x) * PolY.PolPow(y);
		}
	}

	return Res;
}


////////////////////////////////////////////////////////////////////////////////
//
// Derivation
//
//==============================================================================
Polynom Polynom::DerX() const
{
	// Adjust degs
	int MaxDX = 0;
	int MaxDY = 0;
	int MaxDZ = 0;
    for( unsigned int x=1 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
	{
		if( fabs(Monom(x,y,z)) > TOLERANCE )
		{
			MaxDX = x;
			MaxDY = y>MaxDY ? y : MaxDY ;
			MaxDZ = z>MaxDZ ? z : MaxDZ ;
		}
	}

	// Check if X dependent
	if( !MaxDX )
		return 0;

	// Result
    Polynom Res;
	Res.SetDegs( MaxDX-1, MaxDY, MaxDZ );

	// Init degs
    for( unsigned int x=0 ; x<MaxDX   ; x++ )
    for( unsigned int y=0 ; y<MaxDY+1 ; y++ )
    for( unsigned int z=0 ; z<MaxDZ+1 ; z++ )
	{
		Res.Monom(x,y,z) = (x+1) * Monom(x+1,y,z);
	}

	return Res;
}

//==============================================================================
Polynom Polynom::DerY() const
{
	// Adjust degs
	int MaxDX = 0;
	int MaxDY = 0;
	int MaxDZ = 0;
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=1 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
	{
		if( fabs(Monom(x,y,z)) > TOLERANCE )
		{
			MaxDX = x;
			MaxDY = y>MaxDY ? y : MaxDY ;
			MaxDZ = z>MaxDZ ? z : MaxDZ ;
		}
	}

	// Check if Y dependent
	if( !MaxDY )
		return 0;

	// Result
    Polynom Res;
	Res.SetDegs( MaxDX, MaxDY-1, MaxDZ );

	// Init degs
    for( unsigned int x=0 ; x<MaxDX+1 ; x++ )
    for( unsigned int y=0 ; y<MaxDY   ; y++ )
    for( unsigned int z=0 ; z<MaxDZ+1 ; z++ )
	{
		Res.Monom(x,y,z) = (y+1) * Monom(x,y+1,z);
	}

	return Res;
}

//==============================================================================
Polynom Polynom::DerZ() const
{
	// Adjust degs
	int MaxDX = 0;
	int MaxDY = 0;
	int MaxDZ = 0;
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=1 ; z<m_DegZ+1 ; z++ )
	{
		if( fabs(Monom(x,y,z)) > TOLERANCE )
		{
			MaxDX = x;
			MaxDY = y>MaxDY ? y : MaxDY ;
			MaxDZ = z>MaxDZ ? z : MaxDZ ;
		}
	}

	// Check if Z dependent
	if( !MaxDZ )
		return 0;

	// Result
    Polynom Res;
	Res.SetDegs( MaxDX, MaxDY, MaxDZ-1 );

	// Init degs
    for( unsigned int x=0 ; x<MaxDX+1 ; x++ )
    for( unsigned int y=0 ; y<MaxDY+1 ; y++ )
    for( unsigned int z=0 ; z<MaxDZ   ; z++ )
	{
		Res.Monom(x,y,z) = (z+1) * Monom(x,y,z+1);
	}

	return Res;
}

//==============================================================================
Polynom Polynom::IntX() const
{
    // Result
    Polynom Res;
    Res.SetDegs( m_DegX+1, m_DegY, m_DegZ );

    // Init degs
    for( unsigned int x=1 ; x<m_DegX+2 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
    {
        Res.Monom(x,y,z) = (1.0/x) * Monom(x-1,y,z);
    }

    return Res;
}

//==============================================================================
Polynom Polynom::IntY() const
{
    // Result
    Polynom Res;
    Res.SetDegs( m_DegX, m_DegY+1, m_DegZ );

    // Init degs
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=1 ; y<m_DegY+2 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
    {
        Res.Monom(x,y,z) = (1.0/y) * Monom(x,y-1,z);
    }

    return Res;
}

//==============================================================================
Polynom Polynom::IntZ() const
{
    // Result
    Polynom Res;
    Res.SetDegs( m_DegX, m_DegY, m_DegZ+1 );

    // Init degs
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=1 ; z<m_DegZ+2 ; z++ )
    {
        Res.Monom(x,y,z) = (1.0/z) * Monom(x,y,z-1);
    }

    return Res;
}

/**
 * TESTING integration and derivation
 *
    LzMath::Polynom TEST = 1.0 + 2.0*_X_ + 3.0*_Y_ + 4.0*_XY_ + 5.0*_XZ_;
    TEST.Log("TEST ");

    TEST = TEST.IntZ();
    TEST.Log("INT 1");

    TEST = TEST.IntZ();
    TEST.Log("INT 2");

    TEST = TEST.DerZ();
    TEST.Log("DER 1");

    TEST = TEST.DerZ();
    TEST.Log("DER 2"); // Should produce the initial TEST polynom
 *
 */


////////////////////////////////////////////////////////////////////////////////
//
// Integration
//
//==============================================================================
void Polynom::IntegralTo( unsigned int p_MaxDegX, unsigned int p_MaxDegY, unsigned int p_MaxDegZ, double * p_pSums, double & p_Res ) const
{
	// Check available info
	if( m_DegX>p_MaxDegX || m_DegY>p_MaxDegY || m_DegZ>p_MaxDegZ )
        LzLogException("", "[Polynom::IntegralTo] : Integration info not available for all my monoms ! My degs=("<<m_DegX<<","<<m_DegY<<","<<m_DegZ<<")")

	// Integration of all monoms
	p_Res = 0;
    for( unsigned int x=0 ; x<m_DegX+1 ; x++ )
    for( unsigned int y=0 ; y<m_DegY+1 ; y++ )
    for( unsigned int z=0 ; z<m_DegZ+1 ; z++ )
	{
		// If monom non null
		double Coef = Monom(x,y,z);

		if( fabs(Coef) > TOLERANCE )
		{
			p_Res = p_Res + Coef * p_pSums[ x + (p_MaxDegX+1)*(y + (p_MaxDegY+1)*z) ];
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// Debug
//
//==============================================================================
string Polynom::ToString() const
{
    string Pol = "";

	for( int x=m_DegX ; x>=0 ; x-- )
	for( int y=m_DegY ; y>=0 ; y-- )
	for( int z=m_DegZ ; z>=0 ; z-- )
	{
        double lCoefXYZ = Monom(x,y,z);

        if( lCoefXYZ )
		{
            string MonomX;
			if( x )
			{
				if( x == 1)
					MonomX = "X";
				else
                    MonomX = "X^" + std::to_string(x);
			}

            string MonomY;
			if( y )
			{
				if( y == 1)
					MonomY = "Y";
				else
                    MonomY = "Y^" + std::to_string(y);
			}

            string MonomZ;
			if( z )
			{
				if( z == 1)
					MonomZ = "Z";
				else
                    MonomZ = "Z^" + std::to_string(z);
			}

            if( Pol.length() )
				Pol += "+";

            string MonomXYZ;
            MonomXYZ = "("+std::to_string(lCoefXYZ)+")"+MonomX+MonomY+MonomZ;
			Pol += MonomXYZ;
		}
	}

    if( !Pol.length() )
		Pol = "(0)";

	return Pol;
}

//==============================================================================
void Polynom::Log( string p_Prefix/*=""*/ ) const
{
    LzLogM("", p_Prefix<<" : (DX="<<m_DegX<<",DY="<<m_DegY<<",DZ="<<m_DegZ<<") = "<<ToString())
}

}
