#include "Beam.hpp"

const double Eps = 1e3*DBL_EPSILON; // tolerance precision for methods of the class "Beam"

//==============================================================================
// Rotates matrix by angle between e and E
// Parameters: e - old basis,
// E - new basis,
// v - propagation direction.
// rotation matrix:  |  cos(angle)  sin(angle)  |
//                   | -sin(angle)  cos(angle)  |
void  Rot(const Point3D& e, const Point3D& E, const Point3D& v, matrixC &m)
{
	const double Eps = 1e2*DBL_EPSILON; // acceptable precision
	double cs = E*e;
	if(fabs(1.0-cs) < Eps) return; 
	if(fabs(1.0+cs) < Eps) {       
		m[0][0] = -m[0][0];  m[0][1] = -m[0][1];
		m[1][0] = -m[1][0];  m[1][1] = -m[1][1];
		return;
	}
	#ifdef _DEBUG
		if(fabs(cs) > 1.0+DBL_EPSILON)  throw " Rot: Error!";
	#endif 
	double angl = acos(cs); 
    Point3D k = e%E;
	k /= length(k);
	if(norm(k+v) <= 0.5) angl = -angl;   
	double  sn = sin(angl); // the rotation of matrix "m"
	complex b00 = m[0][0]*cs+m[1][0]*sn, // first row of the result
			b01 = m[0][1]*cs+m[1][1]*sn;
	m[1][0] = m[1][0]*cs-m[0][0]*sn;
	m[1][1] = m[1][1]*cs-m[0][1]*sn;
	m[0][0] = b00;  m[0][1] = b01; //
}
//==============================================================================
// members of the class "Beam"
// projection of the incident beam bm (3D polygon) to the plane of the
// 3D facet i and then to 2D plane k
Polyg  Beam::Projection(const double * const pp, int k) const
{
	Polyg g; // result projection
	const Point3D n = Point3D(pp[0],pp[1],pp[2]); // normal to facet
	const double  ang = n*this->r, D = pp[3];
	std::list<Point3D>::const_iterator p = this->v.begin();
	Point3D pp0 = *p, pp1 = *++p, pcnt = CenterOfBeam(*this);

	Point2D p0 = Proj(pp0-this->r*(n*pp0+D)/ang, k),
			p1 = Proj(pp1-this->r*(n*pp1+D)/ang, k),
			cnt = Proj(pcnt-this->r*(n*pcnt+D)/ang, k);

	if(cnt.classify(p0,p1) == RIGHT)
		for(p = this->v.begin(); p!=this->v.end(); p++) {
			const Point3D r = *p-this->r*(n*(*p)+D)/ang;
			const Point2D pp = Proj(r,k);
			g.insert(pp);
		}
	else {
		p = this->v.end();
		do {
			p--;
			// projection of the beam to facet
			const Point3D r = *p-this->r*(n*(*p)+D)/ang;
			const Point2D pp = Proj(r,k);
			g.insert(pp);
		}
		while(p!=this->v.begin());
	}
	return g;
}
//------------------------------------------------------------------------------

// geometrical center of the input beam bm in 3D space
Point3D  CenterOfBeam(const Beam &bm)
{
	Point3D pt(0,0,0);
	for(std::list<Point3D>::const_iterator p=bm.v.begin(); p!=bm.v.end(); p++)
		pt += *p;
	return pt/bm.v.size();
}
//------------------------------------------------------------------------------

// the area of the cross section of the beam
double  CrossSection(const Beam &bm)
{
	const double Eps = 1e7*DBL_EPSILON;
	std::list<Point3D>::const_iterator p = bm.v.begin();
	const Point3D bgn = *p++;
	Point3D p1 = *p++-bgn, p2 = *p-bgn;
	Point3D v = p2%p1;
	const double e = fabs(v*bm.r);
	if(e < Eps) return 0;
	const double lv = length(v);
	double s = 0;
	Point3D d = p2%p1;
	s += length(d);
	while(p++, p!=bm.v.end()) {
		p1 = p2;
		p2 = *p - bgn;
		d = p2%p1;
		s += length(d);
	}
	return e*s/(2.0*lv);
}
//------------------------------------------------------------------------------

// the area of of the beam
double AreaOfBeam(const Beam &bm)
{
	double s = 0;
	std::list<Point3D>::const_iterator p = bm.v.begin();
	const Point3D& bgn = *p++;
	Point3D p1 = *p-bgn;
	while(p++, p!=bm.v.end()) {
		const Point3D p2 = *p - bgn;
		s += length(p2%p1);
		p1 = p2;
	}
	return s/2.0;
}
//------------------------------------------------------------------------------

// output of the vertexes of the beam
void  OutBeam(char* name, const Beam& bm)
{
	unsigned int i = 0;
	std::ofstream out(name);
	out.precision(15);
	for(std::list<Point3D>::const_iterator p=bm.v.begin(); p!=bm.v.end(); p++)
		out << i++ << "  {" << p->x << ',' << p->y << ',' << p->z << "}\n";
}
//------------------------------------------------------------------------------

Point2D  Proj_Vertex(const Beam &bm , Point3D &pt)
{
	const  Point3D p_pr = pt-bm.r*(bm.r*pt); // 
	return Point2D(p_pr*bm.T, p_pr*bm.F);
}
//------------------------------------------------------------------------------

Beam&  Beam::Rotate(Point3D k, Point3D Ey)  
{
	Point3D e1;
	double cs = k*this->r;  
	SphCrd crd = GetSpherical(this->r);
	if(fabs(1.0-cs) <= DBL_EPSILON) e1 = -Ey; 
	else 
		if(fabs(1.0+cs) <= DBL_EPSILON) e1 = Ey; 
		else {
			e1 = Point3D(-sin(crd.fi),cos(crd.fi),0);    
		}
	Rot(this->e, e1, this->r, this->mt);

	this->e=e1;
	return *this;
}
//------------------------------------------------------------------------------

SphCrd   Beam::Spherical(void) const
{ 
	return GetSpherical(this->r); 
}
//------------------------------------------------------------------------------

Beam  Beam::RotatePlane(const Point3D& NewE)   {
	Beam res = *this;
	Rot(res.e, NewE, res.r, res.mt);
	res.e = NewE;
	return res;
}
//------------------------------------------------------------------------------




