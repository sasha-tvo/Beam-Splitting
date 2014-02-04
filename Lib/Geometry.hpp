#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

//==============================================================================
// Basic geometrical objects as points, lines and so on for different geometrical
// algorithms (Chapter 4)
// (c) M.Laszlo "Computational Geometry and Computer Graphics in C++", 1996
// ISBN 0-13-290842-5 
//==============================================================================

#ifndef  _FLOAT_H___
#include <float.h>
#endif

#ifndef  SERVICE_HPP
#include "service.hpp"
#endif

#ifndef MATRIX_HPP
#include "matrix.hpp"
#endif

/** @addtogroup Tracer Beam splitting algorithm
* @{
*/

/** @addtogroup geomstruct Geometrical data structures
* @{
*/


//==============================================================================
#define m_2pi 6.283185307179586476925286766559  // 2*M_PI
const double E_DBL = 1e7*DBL_EPSILON; // tolerance for precision
//=============================================================================
/// relative position of two edges
enum direction  { NON_PARALLEL, OPPOSITE, ACCORDANT, PERPENDICULAR };
/// direction of rotation (p.84)
enum rotation   { CLOCKWISE, COUNTER_CLOCKWISE};
/// relative state of two edges (p.93)
enum state      { COLLINEAR, PARALLEL, SKEW, SKEW_CROSS, SKEW_NO_CROSS };
/// position of point regarding to vector (p.76)
enum position   { LEFT, RIGHT, BEYOND, BEHIND, BETWEEN, ORIGIN, DESTINATION };
/// mutual position of two 2D polygons (p.158)
enum crossing { UNKNOWN, P_IS_INSIDE, Q_IS_INSIDE};

class Edge;
class SphCrd;

//==============================================================================
/**
@brief The class Point2D contains 2D point in Cartesian coordinates.

See M.Laszlo "Computational Geometry and Computer Graphics in C++",
page 73.
*/
class Point2D {
public:
	double x, y;
	// constructor
	Point2D(double _x=0.0, double _y=0.0) : x(_x), y(_y) {};
	// members
	Point2D   operator-(void) const
		{ return Point2D(-this->x, -this->y); }
	Point2D   operator-(const Point2D& a) const
		{ return Point2D(this->x-a.x, this->y-a.y); }
	Point2D   operator-=(const Point2D& a)
		{ return *this = Point2D(this->x-a.x, this->y-a.y); }
	Point2D   operator+(void) const 
		{ return *this; }
	Point2D   operator+(const Point2D& a) const
		{ return Point2D(this->x+a.x, this->y+a.y); }
	Point2D   operator+=(const Point2D& a)
		{ return *this = Point2D(this->x+a.x, this->y+a.y); }
	double    operator*(const Point2D& a) const///< Dot Product of 2D vectors
		{ return this->x*a.x + this->y*a.y; }
	double   operator%(const Point2D& a) const ///< Cross Product of 2D vectors
		{ return this->x*a.y - this->y*a.x; }
	Point2D   operator*(const double& a)  const
		{ return Point2D(this->x*a, this->y*a); }
	Point2D   operator/(const double& a)  const
		{ return Point2D(this->x/a, this->y/a); }
	Point2D   operator/=(const double& a)
		{ return *this = Point2D(this->x/a, this->y/a); }
	bool      operator==(const Point2D& a) const
		{ return this->x==a.x && this->y==a.y; }
	bool      operator!=(const Point2D& a) const
		{ return this->x!=a.x || this->y!=a.y; }
	position  classify(const Edge &) const;
	position  classify(const Point2D&, const Point2D &) const;
	double    distance(const Edge &) const; ///< signed distance
	double    dist(const Edge &) const;     ///< unsigned distance
	// friends
	double   norm(const Point2D& a) 
		{ return a.x*a.x+a.y*a.y; }
	friend double   length(const Point2D& a) 
		{ return hypot(a.x, a.y); }
	friend double   polarAngle(const Point2D &);
	Point2D  normalTo(const Point2D& a)
		{ return Point2D(a.y, -a.x); }
	friend Point2D  operator*(const double& x, const Point2D& a)
		{ return Point2D(x*a.x, x*a.y); }
	friend double   signedAngle(const Point2D&, const Point2D&, const Point2D &);
};

//==============================================================================

/**
@brief The Edge class represents all forms of lines in 2D.

See M.Laszlo "Computational Geometry and Computer Graphics in C++",
page 88.
*/
class Edge {
public:
	Point2D org, dest;
	Edge(const Point2D &_org, const Point2D &_dest) :
		org(_org), dest(_dest) {};
	Edge(void) : 
		org(Point2D(0,0)), dest(Point2D(1.0,0)) {};
	Edge&    operator*=(double x)
		{ return *this = Edge(x*this->org, x*this->dest); }
	Edge     operator*(double x) const
		{ return Edge(x*this->org, x*this->dest); }
	Edge&    rot90(void) /// rotation of the edge through angle 90 deg.
		{ const Point2D m(this->org+this->dest),
						n(this->dest.y-this->org.y, this->org.x-this->dest.x);
			this->org = (m-n)/2.0;  this->dest = (m+n)/2.0;
			return *this;
		}
	Edge&    flip(void) 
		{ return rot90().rot90(); }
	Point2D  point(double t) const
		{ return Point2D(org+t*(dest-org)); }
	state    intersect(const Edge&, double&) const;
	double   slope(void) const { ///< the slope of the edge
		return (this->org.x != this->dest.x) ?
			(this->dest.y-this->org.y)/(this->dest.x-this->org.x) :
			DBL_MAX;
	}
	double   y(double t) const
		{ return this->slope()*(t-this->org.x)+this->org.y; }
	double   signedAngle(const Point2D& a) const;
	// friends of the class "Edge"
	friend double     length(const Edge& e)
		{ return length(e.dest - e.org); }
	friend direction  parallel(const Edge& e, const Edge& f);
	friend state  crossingPoint(const Edge&, const Edge&, Point2D &);
};

//==============================================================================
// 3D point in Cartesian coordinates

/**
@brief The class Point3D contains 3D point in Cartesian coordinates.

See M.Laszlo "Computational Geometry and Computer Graphics in C++",
page 95.
*/
class Point3D {
public:
	double x, y, z;
	Point3D(double X=0, double Y=0, double Z=0) : x(X), y(Y), z(Z) {};
	Point3D  operator-() const
		{ return Point3D(-this->x, -this->y, -this->z); }
	Point3D  operator-(const Point3D& a) const
		{ return Point3D(this->x-a.x, this->y-a.y, this->z-a.z); }
	Point3D  operator-=(const Point3D& a)
		{ return *this = Point3D(this->x-a.x, this->y-a.y, this->z-a.z); }
	Point3D  operator+(const Point3D& a) const
		{ return Point3D(this->x+a.x, this->y+a.y, this->z+a.z); }
	Point3D  operator+=(const Point3D& a)
		{ return *this = Point3D(this->x+a.x, this->y+a.y, this->z+a.z); }
	Point3D  operator*(double a) const
		{ return Point3D(a*this->x, a*this->y, a*this->z); }
	double   operator*(const Point3D& a) const
		{ return this->x*a.x + this->y*a.y + this->z*a.z; } ///< Dot product
	Point3D  operator/(double a) const
		{ return Point3D(this->x/a, this->y/a, this->z/a); }
	Point3D  operator*=(double a)
		{ return *this = Point3D(this->x*a, this->y*a, this->z*a); }
	Point3D  operator/=(double a)
		{ return *this = Point3D(this->x/a, this->y/a, this->z/a); }
	Point3D  operator%(const Point3D& a) const ///< Cross product
		{ return Point3D(this->y*a.z-this->z*a.y, this->z*a.x-this->x*a.z, this->x*a.y-this->y*a.x);}
	Point3D  rotate(double, const Point3D &) const;
	Point3D  rotateAxisY(double a) const;
	Point3D  rotateAxisZ(double a) const;
	Point3D  RotateEuler(double, double, double);
	bool     operator==(const Point3D& a) const
		{ return (fabs(this->x-a.x)<=1e+2*DBL_EPSILON) && (fabs(this->y-a.y)<=1e+2*DBL_EPSILON) && (fabs(this->z-a.z)<=1e+2*DBL_EPSILON); }
	bool     operator!=(const Point3D& a) const { return !(*this==a); }
	friend Point3D  operator*(const Point3D& a, const matrix &m); 
	friend Point3D  operator*(double x, const Point3D& a)
		{ return a*x; }
	friend double   norm(const Point3D& a) ///< Returns square of length
		{ return a.x*a.x + a.y*a.y + a.z*a.z; }
	friend double   length(const Point3D& a) 
		{ return sqrt(norm(a)); }
	friend Point2D  Proj(const Point3D& a, int i){ 
		return (i==0) ? Point2D(a.y, a.z) :
		( (i==1) ? Point2D(a.x, a.z) :
		Point2D(a.x, a.y) ); 
	}
	friend SphCrd GetSpherical(const Point3D& r);
};
//==============================================================================
// 
/**
@brief Direction in 3D space (spherical coordinates on unit sphere).
*/
class SphCrd {
public:
	double fi, tetta;  // spherical coordinates (radians): fi, tetta
	SphCrd(double _fi=0, double _tt=0) : fi(_fi), tetta(_tt) {};
};

//==============================================================================
/**
@brief A pointer-based implementation of a linked list.

See M.Laszlo "Computational Geometry and Computer Graphics in C++",
page 25.
*/
class Node {
protected:
	Node *_next, // reference to the next node
	*_prev; // reference to the previous node
public:
	Node(void) : _next(this), _prev(this) {}; // null list
	virtual  ~Node(void) {};
	Node*     next(void) const { return this->_next; }
	Node*     prev(void) const { return this->_prev; }
	Node*     insert(Node*);
	Node*     remove(void);
	void      splice(Node*);
};
//==============================================================================
// Vertex of the 2D polygon:
// 2D point + two references to next and previous vertices

/**
@brief The class contains vertexes of 2D polygon.

Vertex of the 2D polygon consisits of 2D point and two references to \b next and \b previous vertexes

See M.Laszlo "Computational Geometry and Computer Graphics in C++",
page 80.
*/

class Vertex: public Node, public Point2D {
public:
	Vertex(double x, double y) : Point2D(x,y) {};
	Vertex(const Point2D &p)   : Point2D(p)   {};
	virtual  ~Vertex(void) {};
	Vertex*  cw(void)    const { return (Vertex*)_next; }
	Vertex*  ccw(void)   const { return (Vertex*)_prev; }
	Point2D  point(void) const { return *((Point2D*)this); }
	Vertex*  insert(Vertex *v) { return (Vertex*)(Node::insert(v)); }
	Vertex*  remove(void)      { return (Vertex*)(Node::remove()); }
	void     splice(Vertex *b) { Node::splice(b); }
	Vertex*  neighbor(rotation r) const
		{ return (r == CLOCKWISE ? this->cw() : this->ccw()); }
	Vertex*  split(Vertex*);
	friend class Polyg;
};
//==============================================================================
/**
@brief The Polygon class represents 2D polygon.

See M.Laszlo "Computational Geometry and Computer Graphics in C++",
page 82.
*/
class Polyg {
private:
	mutable Vertex* _v;
	unsigned int _size;
	void  resize(void);
	void  Copy(const Polyg&);
public:
	// constructors
	Polyg(void): _v(NULL), _size(0) {}; // empty polygon
	Polyg(const Polyg& p)    { this->Copy(p); } // copy constructor
	Polyg(Vertex *v) : _v(v) { this->resize(); } //
	~Polyg(void);
	// members
	Polyg&   operator=(const Polyg& p) // the assignment operation
	{ return (this == &p) ? *this : (this->~Polyg(), this->Copy(p), *this); }
	Vertex*  setV(Vertex *v) const     { return this->_v = v; }
	Vertex*  v(void)     const         { return this->_v; }
	unsigned int  size(void)  const    { return this->_size; }
	Vertex*  cw(void)    const         { return this->_v->cw(); }
	Vertex*  ccw(void)   const         { return this->_v->ccw(); }
	Point2D  point(void) const         { return this->_v->point(); }
	Edge     edge(void)  const
		{ return Edge(this->point(), this->_v->cw()->point()); }
	Vertex*  neighbor(rotation r) const
		{ return this->_v->neighbor(r); }
	Vertex*  advance(rotation r) const // перебирает точки полигона
		{ return setV((this->_v->neighbor(r))); }
	Vertex*  insert(const Point2D& p)
		{ return this->_size++ == 0 ? this->_v=new Vertex(p)
	    : this->_v=this->_v->insert(new Vertex(p)); }
	void     remove(void)
		{ Vertex *v=this->_v; this->_v = (--this->_size==0 ? NULL : _v->ccw());
		delete v->remove(); }
	Polyg*   split(Vertex* b)
		{ Vertex *bp=this->_v->split(b); this->resize();
		return new Polyg(bp); }
	bool     Test(void);
	double  PolygonIntersection (const Polyg&, const Polyg&);
	// friends
	friend void  advance(const Polyg& A, Polyg& R, bool inside);
	friend bool    convexPolygonIntersect(const Polyg&, const Polyg&, Polyg &);
	friend bool    FastPolygonIntersect(const Polyg&, const Polyg&, Polyg &);
	friend double  AreaOfConvexPolygon(const Polyg &);
	friend void    OutPolyg(char*, char*, const Polyg &);
	friend bool    CheckPolygon(const Polyg &);
	friend bool    pointInConvexPolygon(const Point2D&, const Polyg &);
	friend bool    clipPolygon(const Polyg&, const Polyg&, Polyg &);
	friend bool    clipPolygonToEdge(const Polyg&, const Edge&, Polyg &);
};
//==============================================================================
/** @}
*/

/** @}
*/
#endif

