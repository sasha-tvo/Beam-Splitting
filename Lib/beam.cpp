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
double  AreaOfBeam(const Beam &bm)
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

void Beam::SetCoefficients_abcd(Point3D& rv, Point3D& Tv, Point3D& Fv, Point3D& r0)
{
	const double E = 1e9*DBL_EPSILON;

	Point3D cnt = Proj(Tv,Fv,rv,r0);

	this->vpr.clear();

	for(std::list<Point3D>::const_iterator p = this->v.begin(); p!=this->v.end(); p++)
	{
		Point3D pt = Proj(Tv,Fv,rv,*p)-cnt;
		this->vpr.push_back(Point2D(pt.x,pt.y));
	}

	this->ab.clear();
	this->cd.clear();

	std::list<Point2D>::const_iterator pr = this->vpr.begin();
	Point2D p1 = *pr++, p2;

	for(unsigned int  i=1, imax=this->vpr.size(); i<=imax; i++)
	{
		p2 = (i!=imax?*pr++:this->vpr.front());
		if(fabs(p1.x-p2.x)>E)
		{
			double ai = (p1.y-p2.y)/(p1.x-p2.x), bi = p1.y - ai*p1.x;
			this->ab.push_back(Point2D(ai,bi));
		}
		else
			this->ab.push_back(Point2D(0,0));
		if(fabs(p1.y-p2.y)>E)
		{
			double ci = (p1.x-p2.x)/(p1.y-p2.y), di = p1.x - ci*p1.y;
			this->cd.push_back(Point2D(ci,di));
		}
		else
			this->cd.push_back(Point2D(0,0));
		p1 = p2;
	}

}
//------------------------------------------------------------------------------

complex   Beam::DiffractionInclinePr(const Point3D& pt, double lam) const
 {
  const double E = 1e9*DBL_EPSILON, E2 = 1e6*DBL_EPSILON; //eps for circle of exactly backscattering

  Point3D pt_proj = Proj(this->N,-pt+this->r);
  const double A = pt_proj.x, B = pt_proj.y;
  complex one(0,-1);
  const double k = m_2pi/lam;
  if(fabs(A)<E2 && fabs(B)<E2) return -one/lam*AreaOfBeam(*this);

  complex s(0,0);
  std::list<Point2D>::const_iterator p = this->vpr.begin();
  Point2D p1 = *p++, p2;

  if(fabs(B)>fabs(A)) {
	std::list<Point2D>::const_iterator ab = this->ab.begin();
	for(unsigned int  i=1, imax=this->vpr.size(); i<=imax; i++) {
	  p2 = (i!=imax?*p++:this->vpr.front());
	  if(fabs(p1.x-p2.x)<E) {ab++; p1 = p2; continue; }
	  const double ai = (*ab).x,
				   bi = (*ab).y,
				   Ci = A+ai*B;
	  ab++;
	  s += exp_im(k*B*bi)*
			  (fabs(Ci)<E ? complex(-k*k*Ci*(p2.x*p2.x-p1.x*p1.x)/2.0,k*(p2.x-p1.x)) :
							(exp_im(k*Ci*p2.x) - exp_im(k*Ci*p1.x))/Ci);
	  p1 = p2;
	 }
	s /= B;
   }
  else {
	std::list<Point2D>::const_iterator cd = this->cd.begin();
	for(unsigned int  i=1, imax=this->vpr.size(); i<=imax; i++) {
	  p2 = (i!=imax?*p++:this->vpr.front());
	  if(fabs(p1.y-p2.y)<E) {cd++; p1 = p2; continue; }
	  const double ci = (*cd).x,
				   di = (*cd).y,
				   Ei = A*ci+B;
	  cd++;
	  s += exp_im(k*A*di) *
		   (fabs(Ei)<E ? complex(-k*k*Ei*(p2.y*p2.y-p1.y*p1.y)/2.0,k*(p2.y-p1.y)) :
						 (exp_im(k*Ei*p2.y) - exp_im(k*Ei*p1.y))/Ei);
	  p1 = p2;
	 }
	s /= -A;
   }
  return one*lam*s/SQR(m_2pi);
 }
//------------------------------------------------------------------------------

complex   Beam::DiffractionShiftedPr(const Point3D& pt, double lam) const
 {
  const double E = 1e9*DBL_EPSILON, E2 = 1e6*DBL_EPSILON; //eps for circle of exactly backscattering

  Point3D pt_proj = Proj(this->T,this->F,this->r, pt);
  const double A = pt_proj.x, B = pt_proj.y;
  complex one(0,-1);
  const double k = m_2pi/lam;
  if(fabs(A)<E2 && fabs(B)<E2) return -one/lam*CrossSection(*this); // ????

  complex s(0,0);
  std::list<Point2D>::const_iterator p = this->vpr.begin();
  Point2D p1 = *p++, p2;

  if(fabs(B)>fabs(A)) {
	std::list<Point2D>::const_iterator ab = this->ab.begin();
	for(unsigned int  i=1, imax=this->vpr.size(); i<=imax; i++) {
	  p2 = (i!=imax?*p++:this->vpr.front());
	  if(fabs(p1.x-p2.x)<E) {ab++; p1 = p2; continue; }
	  const double ai = (*ab).x,
				   bi = (*ab).y,
				   Ci = A+ai*B;
	  ab++;
	  s += exp_im(k*B*bi)*
			  (fabs(Ci)<E ? complex(-k*k*Ci*(p2.x*p2.x-p1.x*p1.x)/2.0,k*(p2.x-p1.x)) :
							(exp_im(k*Ci*p2.x) - exp_im(k*Ci*p1.x))/Ci);
	  p1 = p2;
	 }
	s /= B;
   }
  else {
	std::list<Point2D>::const_iterator cd = this->cd.begin();
	for(unsigned int  i=1, imax=this->vpr.size(); i<=imax; i++) {
	  p2 = (i!=imax?*p++:this->vpr.front());
	  if(fabs(p1.y-p2.y)<E) {cd++; p1 = p2; continue; }
	  const double ci = (*cd).x,
				   di = (*cd).y,
				   Ei = A*ci+B;
	  cd++;
	  s += exp_im(k*A*di) *
		   (fabs(Ei)<E ? complex(-k*k*Ei*(p2.y*p2.y-p1.y*p1.y)/2.0,k*(p2.y-p1.y)) :
						 (exp_im(k*Ei*p2.y) - exp_im(k*Ei*p1.y))/Ei);
	  p1 = p2;
	 }
	s /= -A;
   }
  return one*lam*s/SQR(m_2pi);
 }
//------------------------------------------------------------------------------

complex   Beam::DiffractionShifted(const Point3D& pt, double lam) const
 {
  const double E = 1e9*DBL_EPSILON, E2 = 1e6*DBL_EPSILON; //eps for circle of exactly backscattering

  Point3D r0 = CenterOfBeam(*this),
		  cnt = Proj(this->T,this->F,this->r,r0),
		  pt_proj = Proj(this->T,this->F,this->r, pt);
  const double A = pt_proj.x, B = pt_proj.y;
  //cnt = Point3D(0,0,0);
  complex one(0,-1);
  const double k = m_2pi/lam;
  if(fabs(A)<E2 && fabs(B)<E2) return -one/lam*CrossSection(*this); // ????

  complex s(0,0);
  std::list<Point3D>::const_iterator p = this->v.begin();
  Point3D p1 = Proj(this->T,this->F,this->r, *p++)-cnt, p2;

  if(fabs(B)>fabs(A)) {
	for(unsigned int  i=1; i<=this->v.size(); i++) {
	  p2 = Proj(this->T,this->F,this->r, (i!=this->v.size()?*p++:this->v.front()))-cnt;
	  if(fabs(p1.x-p2.x)<E) { p1 = p2; continue; }
	  const double ai = (p1.y-p2.y)/(p1.x-p2.x),
				   bi = p1.y - ai*p1.x,
				   Ci = A+ai*B;
	  s += exp_im(k*B*bi)*
			  (fabs(Ci)<E ? complex(-k*k*Ci*(p2.x*p2.x-p1.x*p1.x)/2.0,k*(p2.x-p1.x)) :
							(exp_im(k*Ci*p2.x) - exp_im(k*Ci*p1.x))/Ci);
	  p1 = p2;
	 }
	s /= B;
   }
  else {
	for(unsigned int  i=1; i<=this->v.size(); i++) {
	  p2 = Proj(this->T,this->F,this->r, (i!=this->v.size()?*p++:this->v.front()))-cnt;
	  if(fabs(p1.y-p2.y)<E) { p1 = p2; continue; }
	  const double ci = (p1.x-p2.x)/(p1.y-p2.y),
				   di = p1.x - ci*p1.y,
				   Ei = A*ci+B;

		s += exp_im(k*A*di) *
			 (fabs(Ei)<E ? complex(-k*k*Ei*(p2.y*p2.y-p1.y*p1.y)/2.0,k*(p2.y-p1.y)) :
						   (exp_im(k*Ei*p2.y) - exp_im(k*Ei*p1.y))/Ei);
	  p1 = p2;
	 }
	s /= -A;
   }
  return one*lam*s/SQR(m_2pi);
 }
//==============================================================================

complex   Beam::DiffractionIncline(const Point3D& pt, double lam) const
 {
  const double E = 1e9*DBL_EPSILON, E2 = 1e6*DBL_EPSILON;
  Point3D k_k0=-pt+this->r;
  Point3D r0 = CenterOfBeam(*this),
		  pt_proj = Proj(this->N, k_k0), cnt = Proj(this->N, r0);

  const double A = pt_proj.x, B = pt_proj.y;
  complex one(0,-1);
  const double k = m_2pi/lam;
  if(fabs(A)<E2 && fabs(B)<E2) return  -one/lam*AreaOfBeam(*this);

  complex s(0,0);

  std::list<Point3D>::const_iterator p = this->v.begin();
  Point3D p1 = Proj(this->N, *p++)-cnt, p2; //переводим вершины в систему координат грани

  if(fabs(B)>fabs(A)) {
	for(unsigned int  i=1; i<=this->v.size(); i++) {
	  p2 = Proj(this->N, (i!=this->v.size()?*p++:this->v.front()))-cnt;
	  if(fabs(p1.x-p2.x)<E) { p1 = p2; continue; }
	  const double ai = (p1.y-p2.y)/(p1.x-p2.x),
				   bi = p1.y - ai*p1.x,
				   Ci = A+ai*B;
	  s += exp_im(k*B*bi)*
			  (fabs(Ci)<E ? complex(-k*k*Ci*(p2.x*p2.x-p1.x*p1.x)/2.0,k*(p2.x-p1.x)) :
							(exp_im(k*Ci*p2.x) - exp_im(k*Ci*p1.x))/Ci);
	  p1 = p2;
	 }
	s /= B;
   }
  else {
	for(unsigned int  i=1; i<=this->v.size(); i++) {
	  p2 = Proj(this->N, (i!=this->v.size()?*p++:this->v.front()))-cnt;
	  if(fabs(p1.y-p2.y)<E) { p1 = p2; continue; }
	  const double ci = (p1.x-p2.x)/(p1.y-p2.y),
				   di = p1.x - ci*p1.y,
				   Ei = A*ci+B;

		s += exp_im(k*A*di) *
			 (fabs(Ei)<E ? complex(-k*k*Ei*(p2.y*p2.y-p1.y*p1.y)/2.0,k*(p2.y-p1.y)) :
						   (exp_im(k*Ei*p2.y) - exp_im(k*Ei*p1.y))/Ei);
	  p1 = p2;
	 }
	s /= -A;
   }
  return one*lam*s/SQR(m_2pi);
 }
//------------------------------------------------------------------------------

Beam  Beam::RotatePlane(const Point3D& NewE)   {
	Beam res = *this;
	Rot(res.e, NewE, res.r, res.mt);
	res.e = NewE;
	return res;
}
//------------------------------------------------------------------------------




