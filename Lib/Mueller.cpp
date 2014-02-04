//==============================================================================
#include "Mueller.hpp"

//==============================================================================
// deriving Mueller matrix from a given Jones matrix
matrix Mueller(const matrixC& in)
{
	matrix M(4,4);
    const double a11 = norm(in[0][0]), a12 = norm(in[0][1]),
				 a21 = norm(in[1][0]), a22 = norm(in[1][1]);
	//----------------------------------------------------------------------------
	// first & second lines
	double A1 = a11+a21, 	A2 = a12+a22;
	M[0][0] = (A1+A2)/2.0; 	M[0][1] = (A1-A2)/2.0;
	A1 = a11-a21; 			A2 = a12-a22;
	M[1][0] = (A1+A2)/2.0; 	M[1][1] = (A1-A2)/2.0;
	complex C1 = in[0][0]*conj(in[0][1]), 	C2 = in[1][1]*conj(in[1][0]);
	M[0][2] = -real(C1)-real(C2); 			M[0][3] = imag(C2)-imag(C1);
	M[1][2] = real(C2)-real(C1); 			M[1][3] = -imag(C1)-imag(C2);
	//----------------------------------------------------------------------------
	// forth & third lines
	C1 = in[0][0]*conj(in[1][0]); 			C2 = in[1][1]*conj(in[0][1]);
	M[2][0] = -real(C1)-real(C2); 			M[2][1] = real(C2)-real(C1);
	M[3][0] = imag(C1)-imag(C2); 			M[3][1] = imag(C2)+imag(C1);
	C1 = in[0][0]*conj(in[1][1]); 			C2 = in[0][1]*conj(in[1][0]);
	M[2][2] = real(C1)+real(C2); 			M[2][3] = imag(C1)-imag(C2);
	M[3][2] = -imag(C1)-imag(C2); 			M[3][3] = real(C1)-real(C2);
	return M;
}

//==============================================================================

void RightRotateMueller(matrix& m, double cs, double sn)
{
	// do not change first and last columns
	const double A01 = m[0][1], A11 = m[1][1], A21 = m[2][1], A31 = m[3][1];
	// second column
	m[0][1] = A01*cs - m[0][2]*sn;
	m[1][1] = A11*cs - m[1][2]*sn;
	m[2][1] = A21*cs - m[2][2]*sn;
	m[3][1] = A31*cs - m[3][2]*sn;
	// third column
	m[0][2] = A01*sn + m[0][2]*cs;
	m[1][2] = A11*sn + m[1][2]*cs;
	m[2][2] = A21*sn + m[2][2]*cs;
	m[3][2] = A31*sn + m[3][2]*cs;
}

//==============================================================================
//
void LeftRotateMueller(matrix& m, double cs, double sn)
{
	// do not change first and last rows
	const double A10 = m[1][0], A11 = m[1][1], A12 = m[1][2], A13 = m[1][3];
	// second row
	m[1][0] = A10*cs + m[2][0]*sn;
	m[1][1] = A11*cs + m[2][1]*sn;
	m[1][2] = A12*cs + m[2][2]*sn;
	m[1][3] = A13*cs + m[2][3]*sn;
	// third row
	m[2][0] = m[2][0]*cs-A10*sn;
	m[2][1] = m[2][1]*cs-A11*sn;
	m[2][2] = m[2][2]*cs-A12*sn;
	m[2][3] = m[2][3]*cs-A13*sn;
}

//==============================================================================
void
RotateMueller(matrix& m, double _cs, double _sn, double cs, double sn)
{
	// corner elements does not change!

	// first column
	double A1 = m[1][0], A2 = m[2][0];
	m[1][0] = A1*_cs + A2*_sn;
	m[2][0] = A2*_cs - A1*_sn;
	// second column
	double B01 = m[0][1], B11 = m[1][1], B21 = m[2][1], B31 = m[3][1];
	A1 = B11*cs - m[1][2]*sn;
	A2 = B21*cs - m[2][2]*sn;
	m[0][1] = B01*cs - m[0][2]*sn;
	m[1][1] = A1*_cs + A2*_sn;
	m[2][1] = A2*_cs - A1*_sn;
	m[3][1] = B31*cs - m[3][2]*sn;
	// third column
	A1 = B11*sn + m[1][2]*cs;
	A2 = B21*sn + m[2][2]*cs;
	m[0][2] = B01*sn + m[0][2]*cs;
	m[1][2] = A1*_cs + A2*_sn;
	m[2][2] = A2*_cs - A1*_sn;
	m[3][2] = B31*sn + m[3][2]*cs;
	// forth column
	A1 = m[1][3]; A2 = m[2][3];
	m[1][3] = A1*_cs + A2*_sn;
	m[2][3] = A2*_cs - A1*_sn;
} //*/
//==============================================================================
void ForwardScattering(matrix& m)
{
	double A = (m[1][1]+m[2][2])/2.0;
	m[1][1] = m[2][2] = A;

	m[0][1] = m[0][2] = 0;
	m[1][0] = m[1][2] = m[1][3] = 0;
	m[2][0] = m[2][1] = m[2][3] = 0;
	m[3][1] = m[3][2] = 0;
}
//==============================================================================
void BackwardScattering(matrix& m)
{
	double A = (m[1][1]-m[2][2])/2.0;
	m[1][1] = A; m[2][2] = -A;
	
	m[0][1] = m[0][2] = 0;
	m[1][0] = m[1][3] = 0;
	m[2][0] = m[2][3] = 0;
	m[3][1] = m[3][2] = 0;
}
//==============================================================================

