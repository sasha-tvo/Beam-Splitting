#include "compl.hpp"

//==============================================================================

// returns sine and cosine of an angle < 28.07.2002 >
// Copyright (C) 1996,99 Inprise Corporation
// SinCos is 2x faster than calling Sin and Cos separately for the same angle
void  SinCos(long double x, long double &sn, long double &cs)
{
	sn = sin(x); cs = cos(x);
}
//------------------------------------------------------------------------------

// squre root of the complex number z < 28.07.2002 >
// (c) Press W.H., Teukolsky S.A., Flannery B.P.
//     NUMERICAL RECIPES IN C: the art of scientific computing,
//     - 2nd edition, p. 507
complex  sqrt(const complex &z)
{
	if(z == 0.0) return 0;
	double	r = fabs(real(z)), i = fabs(imag(z)), t,
			w = (r>=i) ? (t=i/r, sqrt(r)*sqrt(0.5*(1.0+sqrt(1.0+t*t)))):
				(t=r/i, sqrt(i)*sqrt(0.5*(1.0+sqrt(1.0+t*t))));

	return real(z)>0 ? complex(w,imag(z)/(2.0*w)) :
			((t = imag(z)>=0.0 ? w : -w), complex(imag(z)/(2.0*t),t));
}
//------------------------------------------------------------------------------

complex  exp_im(double x)
{ 
	return complex(cos(x), sin(x)); 
};
//------------------------------------------------------------------------------
