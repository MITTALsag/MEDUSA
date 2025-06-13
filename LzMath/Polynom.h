#pragma once

#include <string>

namespace LzGeom { class Point3D; }

namespace LzMath
{
using LzGeom::Point3D;
using std::string;


//------------------------------------------------------------------------------
// POLYNOMS
//------------------------------------------------------------------------------
class Polynom
{
public:
	// Construction / destruction
    Polynom();
    Polynom( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ, double p_Coef );
    Polynom( const Polynom & p_Poly );
    Polynom( double p_Monom0 );
    virtual ~Polynom();

	// Set
    void SetDegs( unsigned int p_DegX, unsigned int p_DegY,unsigned  int p_DegZ );
    double & Monom( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ );
    const double & Monom( unsigned int p_DegX, unsigned int p_DegY, unsigned int p_DegZ ) const;

	// Operators
    Polynom operator*( double p_Scal ) const;
    Polynom operator/( double p_Scal ) const;
    friend Polynom operator*( double p_Scal, const Polynom & p_Poly ) { return p_Poly*p_Scal; }
    friend Polynom operator+( double p_Scal, const Polynom & p_Poly ) { return p_Poly+p_Scal; }
    friend Polynom operator-( double p_Scal, const Polynom & p_Poly ) { return p_Poly*(-1)+p_Scal; }
    friend Polynom operator-( const Polynom & p_Poly ) { return p_Poly*(-1); }
    Polynom operator*( const Polynom & p_Poly ) const;
    Polynom operator+( double p_Scal ) const;
    Polynom operator-( double p_Scal ) const;
    Polynom operator+( const Polynom & p_Poly ) const;
    Polynom operator-( const Polynom & p_Poly ) const;
//	Polynom & operator+=( const Polynom & p_Poly );
    Polynom & operator=( const Polynom & p_Poly );
    Polynom & operator=( double p_Monom0 );
    Polynom PolPow( int Pow ) const;

	// Evaluation
	double EvalAt( double X ) const;
	double EvalAt( double X, double Y ) const;
	double EvalAt( double X, double Y, double Z ) const;
    double EvalAt( const Point3D & pPt3D ) const;
    bool IsZero( double tol=TOLERANCE ) const;

	// Composition
    Polynom Compose( const Polynom & PolX, const Polynom & PolY ) const;
    Polynom Compose( const Polynom & PolX, const Polynom & PolY, const Polynom & PolZ ) const;

	// Derivation
    Polynom DerX() const;
    Polynom DerY() const;
    Polynom DerZ() const;

    // Integration
    Polynom IntX() const;
    Polynom IntY() const;
    Polynom IntZ() const;

	// Integration
    void IntegralTo( unsigned int p_MaxDegX, unsigned int p_MaxDegY, unsigned int p_MaxDegZ, double * p_pSums, double & p_Res ) const;

	// Monoms
    static const Polynom XYZ;
    static const Polynom XY;
    static const Polynom XZ;
    static const Polynom YZ;
    static const Polynom X;
    static const Polynom Y;
    static const Polynom Z;

	// Debug
    string ToString() const;
    void Log( string p_Prefix="" ) const;

	// What is zero and what is not ....
	static const double TOLERANCE;


protected:
	void Free();
	
    unsigned int m_DegX;
    unsigned int m_DegY;
    unsigned int m_DegZ;

    unsigned int m_Monoms;

    double * m_pCoefs; /// utiliser un std vector ici
};


//------------------------------------------------------------------------------
#define _XYZ_ (LzMath::Polynom::XYZ)
#define _XY_  (LzMath::Polynom::XY)
#define _XZ_  (LzMath::Polynom::XZ)
#define _YZ_  (LzMath::Polynom::YZ)
#define _X_   (LzMath::Polynom::X)
#define _Y_   (LzMath::Polynom::Y)
#define _Z_   (LzMath::Polynom::Z)
//------------------------------------------------------------------------------

}
