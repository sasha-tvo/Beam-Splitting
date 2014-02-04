#ifndef SERVICE_HPP
#define SERVICE_HPP


//=============================================================================

template<class T> T     SQR(const T& x)  { return x*x; }
//-----------------------------------------------------------------------------

template<class T> T     CUB(const T& x)  { return x*x*x; }
//-----------------------------------------------------------------------------

template<class T> T     ABS(const T& x)  { return x>0 ? x : -x; }
//-----------------------------------------------------------------------------

template<class T> int   SGN(const T& x)
 { return x>0 ? 1 : (x<0?-1:0); }
//-----------------------------------------------------------------------------

template<class T> void  SWAP(T& x, T& y)
 { const T tmp = y; y = x; x = tmp; }
//-----------------------------------------------------------------------------

template<class T> int   SCAL(const T& x)
 { return x<0 ? -1 : 1; } 
//-----------------------------------------------------------------------------

template<class T> T  MIN(const T& x, const T& y) { return x>y?y:x; }
//-----------------------------------------------------------------------------

template<class T> T  MAX(const T& x, const T& y) { return x<y?y:x; }
//-----------------------------------------------------------------------------

template<class T> T  MAX(const T& x, const T& y, const T& z)
 { return x>y ? (x>z?x:z) : (y>z?y:z); }
//-----------------------------------------------------------------------------

template<class T> T  MIN(const T& x, const T& y, const T& z)
 { return x<y ? (x<z?x:z) : (y<z?y:z); }
//-----------------------------------------------------------------------------

template<class T> int  nMX(const T& x, const T& y) { return x>y?0:1; }
//-----------------------------------------------------------------------------

template<class T> int  nMN(const T& x, const T& y) { return x<y?0:1; }
//-----------------------------------------------------------------------------

template<class T> int  nMX(const T& x, const T& y, const T& z)
 { return x>y ? (x>z?0:2) : (y>z?1:2); }
//-----------------------------------------------------------------------------

template<class T> int  nMN(const T& x, const T& y, const T& z)
 { return x<y ? (x<z?0:2) : (y<z?1:2); }
//-----------------------------------------------------------------------------

#endif
