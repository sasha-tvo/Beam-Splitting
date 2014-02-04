#ifndef COMPL_HPP
#define COMPL_HPP

#ifndef __MATH_H
#include <math.h>
#endif

#ifndef SERVICE_HPP
#include "service.hpp"
#endif

/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */

/** @addtogroup AxFunc Axillary functions
 * @{
 */


class complex;
complex  exp(const complex& z);
complex  exp_im(double x);
void  SinCos(long double, long double&, long double &);
//-----------------------------------------------------------------------------
/**
 @brief This class provides a complex numbers and operation with them
*/
class complex {
	double re, im;
public:
	complex(double r = 0, double i = 0) : re(r), im(i) {}
	// members
	complex&  operator=(const complex& b)
		{ this->re = b.re; this->im = b.im; return *this; }
	complex  operator+()         const { return *this; }
	complex  operator+(double x) const
		{ return complex(this->re+x, this->im); }
	complex  operator+(const complex& z) const
		{ return complex(this->re+z.re, this->im+z.im); }
	complex  operator+=(const complex& z)
		{ return *this = complex(this->re+z.re, this->im+z.im); }
	complex  operator+=(double x)
		{ return *this = complex(this->re+x, this->im); }
	complex  operator-() const
		{ return complex(-this->re, -this->im); }
	complex  operator-(double x) const
		{ return complex(this->re-x, this->im); }
	complex  operator-(const complex& z) const
		{ return complex(this->re-z.re, this->im-z.im); }
	complex  operator-=(const complex& z)
		{ return *this = complex(this->re-z.re, this->im-z.im); }
	complex  operator-=(double x)
		{ return *this = complex(this->re-x, this->im); }
	complex  operator*(double x) const
		{ return complex(x*this->re, x*this->im); }
	complex  operator*(const complex& z) const
		{ 
			return complex(this->re*z.re - this->im*z.im,
			this->im*z.re + this->re*z.im); 
		}
	complex  operator*=(const complex& z)
		{ 
			return *this = complex(this->re*z.re - this->im*z.im,
			this->im*z.re + this->re*z.im); 
		}
	complex  operator*=(double x)
		{ return *this = complex(x*this->re, x*this->im); }
	complex  operator/(double x) const
		{ return complex(this->re/x, this->im/x); }
	complex  operator/(const complex &z) const
		{
			double r, den;
			return (fabs(z.re) >= fabs(z.im)) ?
			(r=z.im/z.re, den=z.re+r*z.im,
			complex(this->re+r*this->im, this->im-r*this->re)/den):
			(r=z.re/z.im, den=z.im+r*z.re,
			complex(r*this->re+this->im, r*this->im-this->re)/den); 
		}
	complex  operator/=(double x)
		{ return *this = complex(this->re/x, this->im/x); };
	complex  operator/=(const complex &z) { return *this = (*this)/z; }
	bool  operator==(const complex& z) const
		{ return (z.re==this->re) && (z.im==this->im); }
	bool  operator!=(const complex& z) const
		{ return (z.re!=this->re) || (z.im!=this->im); }
	bool  operator==(double x) const
		{ return (this->re==x) && (this->im==0); }
	bool  operator!=(double x) const
		{ return (this->re!=x) || (this->im!=0); }
	// friends
	friend complex  operator+(double x, const complex& z)
		{ return complex(x+z.re, z.im); }
	friend complex  operator-(double x, const complex& z)
		{ return complex(x-z.re,-z.im); }
	friend complex  operator/(double x, const complex &z)
		{ 
			double r, den;
			complex c = (fabs(z.re) >= fabs(z.im)) ?
			(r=z.im/z.re, den=z.re+r*z.im, complex(x, -r*x)):
			(r=z.re/z.im, den=z.im+r*z.re, complex(r*x, -x));
			return c/den; 
		}
	friend double&  real(complex& z)       { return z.re; }
	friend double&  imag(complex& z)       { return z.im; }
	friend double   real(const complex& z) { return z.re; }
	friend double   imag(const complex& z) { return z.im; }
	friend double   norm(const complex& z)
		{ return SQR(z.re) + SQR(z.im); }
	friend complex  conj(const complex& z)
		{ return complex(z.re, -z.im); }
	friend complex  exp(const complex& z)
		{ return exp(z.re)*exp_im(z.im); }
	friend complex  operator*(double x, const complex& z)
		{ return complex(x*z.re, x*z.im); }
	friend double  arg(const complex &z)
		{ return z!=0 ? atan2(z.re,z.im) : 0; }
	friend double   abs(const complex &z) { return hypot(z.re, z.im); }
	friend complex  sqrt(const complex &);
 };
//==============================================================================

/** @}
*/

/** @}
*/

#endif
