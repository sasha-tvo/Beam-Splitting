#ifndef PARTICLE_H
#define PARTICLE_H

#ifndef CRYSTAL_H
#include "Crystal.hpp"
#endif

//----------------------------------------------------------------------
/** @addtogroup Tracer Beam splitting algorithm
* @{
*/

/** @addtogroup particles Crystal particles
* @{
*/

/**
@brief This class defines hexagonal prism.

The prism is defined by its radius and half height
*/
class Prism : public Crystal {
	double w, h;       ///< sizes of prism: w - radius and h - half height
	double tt, ps, fi; // orientation of prism   
	double  Third(void) { throw " Prism::Third(): Error! "; }
	double  Forth(void) { throw " Prism::Forth(): Error! "; }
public:
	// constructor and destructor
    Prism(const complex& r, double W, double H, unsigned int Itr, Point3D _k, Point3D _Ey) :
        Crystal(r,12,8,6,Itr,_k, _Ey), w(W), h(H), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets(); }
	virtual  ~Prism(void) {}
	// parameters
	double   First(void)  { return this->w; }
	double   Second(void) { return this->h; }
	void     SetVertices(void); // sets vertices of prism
	void     SetFacets(void);   // sets facets of prism
};

///This class defines the pyramid
class Pyramid : public Crystal {
	double w, h;       ///< sizes of Pyramid: w - radius and h - half height
	double tt, ps, fi; // orientation of prism
	double  Third(void) { throw " Pyramid::Third(): Error! "; }
	double  Forth(void) { throw " Pyramid::Forth(): Error! "; }
public:
    // constructor and destructor
    Pyramid(const complex& r, double W, double H, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,7,7,6,Itr,_k, _Ey), w(W), h(H), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets(); this->ChangeFacets(); }
	virtual  ~Pyramid(void) {}
	// parameters
	double  First(void)  { return this->w; }
	double  Second(void) { return this->h; }
	void    SetVertices(void); // sets vertices of prism
	void    SetFacets(void); // sets facets of prism
};
//==============================================================================
///This class defines the hexagonal bullet
class Bullet : public Crystal {
	double w, h, hp;  ///< sizes of bullet: w - radius, h - half height, hp - height of tip
	double tt, ps, fi; // orientation of bullet (angles)
	double  Forth(void) { throw " Bullet::Forth(): Error! "; }
public:
	// constructor and destructor
    Bullet(const complex& r, double W, double H, double HP, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,13,13,6,Itr,_k, _Ey), w(W), h(H), hp(HP), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets(); this->ChangeFacets(); }
	virtual  ~Bullet(void) {}
	// parameters
	double  First(void)  { return this->w; }
	double  Second(void) { return this->h; }
	double  Third(void)  { return this->hp; }
	void    SetVertices(void); // sets vertices of prism
	void    SetFacets(void);   // sets facets of prism
};
//==============================================================================
///This class defines the hexagonal tapered prism
class TaperedPrism : public Crystal {
	double 	w,    ///< w - radius of the prism
			h,    ///< h - semiheight of the crystal
			wp,   ///< wp - radius of the crystal tip
			hp;   ///< hp - height of the crystal tip
	double tt, ps, fi;   // the orientation of prism
public:
	// constructor and destructor
    TaperedPrism (const complex& r, double W, double H, double WP, double HP, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,24,20,6,Itr,_k, _Ey), w(W), h(H), wp(WP), hp(HP), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets(); this->ChangeFacets(); }
	virtual  ~TaperedPrism(void) {}
	// parameters
	double  First(void)  { return this->w; }
	double  Second(void) { return this->h; }
	double  Third(void)  { return this->hp; }
	double  Forth(void)  { return this->wp; }
	void    SetVertices(void); // sets vertices of prism
	void    SetFacets(void);   // sets facets of prism
};
//==============================================================================
///This class defines the hexagonal prism with cavity on bottom
class CavityPrism : public Crystal {
	double w, h, c;    ///< sizes: w - radius, h - half height and c - cavity depth
	double tt, ps, fi; // orientation of prism
	double  Forth(void) { throw " CavityPrism::Forth(): Error! "; }
public:
	// constructor and destructor
    CavityPrism(const complex& r, double W, double H, double C, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,13,13,6,Itr,_k, _Ey), w(W), h(H), c(C), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets();  this->ChangeFacets();}
	virtual  ~CavityPrism(void) {}
	// parameters
	double  First(void)  { return this->w; }
	double  Second(void) { return this->h; }
	double  Third(void)  { return this->c; }
	void    SetVertices(void); // sets vertices of prism
	void    SetFacets(void);   // sets facets of prism
};

//==============================================================================
///This class defines the cupped prism
class Cup : public Crystal {
	double 	wp, hp, ///< sizes: wp - radius prism, hp - half height of prism
			wc, hc; ///< sizes: wc - radius cup, hc - half height of cup
	double tt, ps, fi; // orientation of cup
public:
	// constructor and destructor
    Cup(const complex& r, double WP, double HP, double WC, double HC, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,18,14,6,Itr,_k, _Ey), wp(WP), hp(HP), wc(WC), hc(HC), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets();  this->ChangeFacets();}
	virtual  ~Cup(void) {}
	// parameters
	double  First(void)  { return this->wp; }
	double  Second(void) { return this->hp; }
	double  Third(void)  { return this->wc; }
	double  Forth(void)  { return this->hc; }
	void    SetVertices(void); // sets vertices of prism
	void    SetFacets(void);   // sets facets of prism
};

///This class defines the triangle prism
class TrianglePrism : public Crystal {
	double w, h;       // sizes of TrianglePrism: w - radius and h - half height
	double tt, ps, fi; // orientation of TrianglePrism
	double  Third(void) { throw " TrianglePrism::Third(): Error! "; }
	double  Forth(void) { throw " TrianglePrism::Forth(): Error! "; }
public:
	// constructor and destructor
    TrianglePrism(const complex& r, double W, double H, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,6,5,4,Itr,_k, _Ey), w(W), h(H), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets();  this->ChangeFacets();}
	virtual  ~TrianglePrism(void) {}
	// parameters
	double   First(void)  { return this->w; }
	double   Second(void) { return this->h; }
	void     SetVertices(void); // sets vertices of TrianglePrism
	void     SetFacets(void);   // sets facets of TrianglePrism
};
//==============================================================================
///This class defines the parallelepiped
class Parallelepiped : public Crystal {
	double w, h, wp;   ///< sizes of Parallelepiped
	double tt, ps, fi; // orientation of Parallelepiped
	double  Third(void) { throw " Parallelepiped::Third(): Error! "; }
public:
	// constructor and destructor
    Parallelepiped(const complex& r, double W, double H, double WP, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,8,6,4,Itr,_k, _Ey), w(W), h(H), wp(WP), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets();  this->ChangeFacets();}
	virtual  ~Parallelepiped(void) {}
	// parameters
	double   First(void)  { return this->w; }
	double   Second(void) { return this->h; }
	double   Forth(void)   { return this->wp; }
	void     SetVertices(void); // sets vertices of Parallelepiped
	void     SetFacets(void);   // sets facets of Parallelepiped
};
//==============================================================================
///This class defines the tetrahedron
class Tetrahedron : public Crystal {
	double w;          // size of Tetrahedron: w - side
	double tt, ps, fi; // orientation of Parallelepiped
	double  Second(void) { throw " Tetrahedron::Second(): Error! "; }
	double  Third(void) { throw " Tetrahedron::Third(): Error! "; }
	double  Forth(void) { throw " Tetrahedron::Forth(): Error! "; }
public:
	// constructor and destructor
    Tetrahedron(const complex& r, double W, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,4,4,3,Itr,_k, _Ey), w(W), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets();  this->ChangeFacets();}
	virtual  ~Tetrahedron(void) {}
	// parameters
	double   First(void)  { return this->w; }
	void     SetVertices(void); // sets vertices of Tetrahedron
	void     SetFacets(void);   // sets facets of Tetrahedron
};
//==============================================================================
///This class defines the dodecahedron
class Dodecahedron : public Crystal {
	double w, h;       // sizes of dodecahedron: w - radius and h -half height
	double tt, ps, fi; // orientation of prism
	double  Third(void) { throw " Dodecahedron::Third(): Error! "; }
	double  Forth(void) { throw " Dodecahedron::Forth(): Error! "; }
public:
	// constructor and destructor
    Dodecahedron(const complex& r, double W, double H, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Crystal(r,24,14,12,Itr,_k, _Ey), w(W), h(H), tt(0), ps(0), fi(0)
		{ this->SetVertices(); this->SetFacets();  this->ChangeFacets();}
	virtual  ~Dodecahedron(void) {}
	// parameters
	double   First(void)  { return this->w; }
	double   Second(void) { return this->h; }
	void     SetVertices(void); // sets vertices of prism
	void     SetFacets(void); // sets facets of prism
};

/** @}
*/

/** @}
*/

#endif // PARTICLE_H
