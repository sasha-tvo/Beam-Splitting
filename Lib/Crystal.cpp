#include "Crystal.hpp"

//==============================================================================

void BeamSplitting::AllocMem(void) // allocation of memory for particle
{
	this->Facets = new unsigned int*[this->Itr+1];
	for(unsigned int i=0;i<this->Itr+1; i++)
		this->Facets[i] = new unsigned int[this->NoF];

	for(unsigned int i=0;i<this->Itr+1; i++)
		for(unsigned int j=0;j<this->NoF;j++)
			this->Facets[i][j] = 1;
}

void BeamSplitting::FreeMem(void) // allocation of memory for particle
{
	for(unsigned int i=0;i<this->Itr+1; i++)
		delete[] this->Facets[i];
	delete[]  this->Facets;

}
// members of "Frame"
void Frame::AllocMem(void) // allocation of memory for particle
{
	this->p = new Point3D[this->M];
	this->Gr = new int*[this->K];
	this->Gran = new double*[this->K];
	for(unsigned int i=0; i<this->K; i++) {
		this->Gr[i] = new int[this->Km+1];
		this->Gran[i] = new double[4];
    }
}

void Frame::FreeMem(void) // memory release
{
	delete[] this->p;
	for(unsigned int i=0; i<this->K; i++) {
		delete[] this->Gran[i];
		delete[] this->Gr[i];
	}
	delete[] this->Gran; delete[] this->Gr;
}

Frame Frame::operator=(const Frame &f)
{
	if(this == &f) return *this;
	this->FreeMem();
	this->M=f.M; this->K=f.K; this->Km=f.Km;
	this->AllocMem();
	for(unsigned int i=0; i<this->M; i++) this->p[i] = f.p[i];
	for(unsigned int i=0; i<this->K; i++) {
		for(unsigned int j=0; j<this->Km+1; j++) {
			this->Gr[i][j] = f.Gr[i][j];
			if(j<4) this->Gran[i][j] = f.Gran[i][j];
		}
	}
	return *this;
}
//------------------------------------------------------------------------------

// members of the class "Crystal"

// projection of the facet i to 2D plane k as a convex polygon
Polyg Crystal::Projection(unsigned int i, int k) const
{
	const double Eps = 1e3*DBL_EPSILON;
	Polyg g;
	#ifdef _DEBUG
	if(k > 2 || i > this->K) throw " Crystal::Projection : Error! ";
	#endif
	const Point2D 	p0 = Proj(this->p[this->Gr[i][0]],k),
					p1 = Proj(this->p[this->Gr[i][1]],k),
					p2 = Proj(this->p[this->Gr[i][2]],k),
					v1 = p1 - p0,
					v2 = p2 - p0;
	const double e = v1.x*v2.y - v1.y*v2.x; ;
	if(e < -Eps) {
		g.insert(p0);
		g.insert(p1);
		g.insert(p2);
		for(unsigned int j=3; this->Gr[i][j+1]!=-1 && j<this->Km; j++) // vertex loop
			g.insert(Proj(this->p[this->Gr[i][j]],k));
	}
	else
		for(unsigned int j=Km; j>=1; j--) // vertex loop
			if(this->Gr[i][j] == -1) continue;
			else g.insert(Proj(this->p[this->Gr[i][j]],k));
	return g;
}

// recovering the 3D result polygon from its 2D projection to the facet i
std::list<Point3D> 
Crystal::Retrieve(Polyg& p, unsigned int i, int k) const
{
	#ifdef _DEBUG
	if(k > 2 || i > this->K) throw " Crystal::Retrieve : Error! ";
	#endif
	std::list<Point3D> g;
	const double A = this->Gran[i][0], B = this->Gran[i][1],
							 C = this->Gran[i][2], D = this->Gran[i][3];
	for(unsigned int n=0; n<p.size(); n++, p.advance(CLOCKWISE)) {
		Point2D pt = p.point();
		Point3D res;
		switch(k) {
			case 0: res = Point3D((-D-B*pt.x-C*pt.y)/A, pt.x, pt.y); break; // yz
			case 1: res = Point3D(pt.x, (-D-A*pt.x-C*pt.y)/B, pt.y); break; // xz
			case 2: res = Point3D(pt.x, pt.y, (-D-A*pt.x-B*pt.y)/C); break; // xy
		}
		g.push_back(res);
	}
	return g;
}

// intersection of the incident beam and projection of facet
bool Crystal::Intersection
(const Beam& bm, unsigned int i, std::list<Point3D>& Result) const
{
	const Point3D n = this->NormToFacet(i);
	// kind of the projection
	const int fl = nMX(fabs(n.x*bm.r.x), fabs(n.y*bm.r.y), fabs(n.z*bm.r.z));
	Polyg 	facet = this->Projection(i, fl), // projection of the facet
			beam = bm.Projection(Facet(i), fl); // projection of the beam
	//----------------------------------------------------------------------------
	Point2D max_fct(facet.point()), min_fct(facet.point()),
			max_bm(beam.point()), min_bm(beam.point());
	for(unsigned int i=0; i<facet.size(); i++, facet.advance(CLOCKWISE)) {
		const Point2D& p = facet.point();
		if(p.x > max_fct.x) max_fct.x = p.x;
		if(p.y > max_fct.y) max_fct.y = p.y;
		if(p.x < min_fct.x) min_fct.x = p.x;
		if(p.y < min_fct.y) min_fct.y = p.y;
	}
	for(unsigned int i=0; i<beam.size(); i++, beam.advance(CLOCKWISE)) {
		const Point2D& p = beam.point();
		if(p.x > max_bm.x) max_bm.x = p.x;
		if(p.y > max_bm.y) max_bm.y = p.y;
		if(p.x < min_bm.x) min_bm.x = p.x;
		if(p.y < min_bm.y) min_bm.y = p.y;
	}
	if(max_fct.x<=min_bm.x || max_bm.x<=min_fct.x ||
		 max_fct.y<=min_bm.y || max_bm.y<=min_fct.y) return false;
	Polyg intersection;
	#ifdef _DEBUG
	const double e = AreaOfConvexPolygon(beam),
				 c = AreaOfConvexPolygon(facet);
	if(!CheckPolygon(facet)) {
		OutPolyg("Errors.log", "facet", facet);
		throw "Crystal::Intersection : Error!";
	}
	if(!CheckPolygon(beam)) {
		OutPolyg("Errors.log", "beam", beam);
		throw "Crystal::Intersection : Error!";
	} 
	#endif
    if(FastPolygonIntersect(beam, facet, intersection)) {
		#ifdef _DEBUG
		if(!CheckPolygon(intersection)) {
			OutPolyg("Errors.log", "beam", beam);
			OutPolyg("Errors.log", "facet", facet);
			OutPolyg("Errors.log", "intersection", intersection);
			throw "Crystal::Intersection : Error!";
		} 
		#endif
		const double e = AreaOfConvexPolygon(intersection);
        //if(e < S_eps) return false;
		Result = this->Retrieve(intersection, i, fl);
		//------------------------------------------------------------------------
		return true;
	}
	//----------------------------------------------------------------------------
	return false; // the intersection is empty
}

// center of gravity for facet ii
Point3D Crystal::CenterOfFacet(unsigned int ii) const
{
	Point3D pt(0,0,0);
	unsigned int j;
	for(j=0; this->Gr[ii][j+1]!=-1 && j<this->Km; j++)
		pt += this->p[this->Gr[ii][j]];
	return pt/j;
}

double Crystal::AreaOfFacet(unsigned int ii) const
{
	#ifdef _DEBUG
	if(ii > this->K || this->Km<3) throw " Crystal::Retrieve : Error! ";
	#endif
	double s = 0;
	const Point3D bgn = this->p[this->Gr[ii][0]];
	Point3D p1 = this->p[this->Gr[ii][1]]-bgn;
	unsigned int j=2;
	while(j<this->Km) {
		const Point3D p2 = this->p[this->Gr[ii][j]]- bgn;
		s += length(p2%p1);
		p1 = p2;
		j++;
	}
	return s/2;
}

// full recalculation of coordinates of all crystal vertices
// for new Euler angles tt, ps, f
void Crystal::ChangeVertices(double tt, double ps, double f)
{
	#ifdef __FASTMATH_H
	double c1, c2, c3, s1, s2, s3;
	sincos(tt, &s1, &c1); sincos(ps, &s2, &c2); sincos(f, &s3, &c3);
	#else
	const double 	c1 = cos(tt), c2 = cos(ps), c3 = cos(f),
					s1 = sin(tt), s2 = sin(ps), s3 = sin(f);
	#endif
	const double 	bf1 = c1*s2, bf2 = c1*c2,
					l1 = c2*c3-bf1*s3, l2 = -bf1*c3-c2*s3, l3 = s1*s2,
					m1 = s2*c3+bf2*s3, m2 = bf2*c3-s2*s3, m3 = -s1*c2,
					n1 = s1*s3, n2 = s1*c3, n3 = c1;
	Point3D buf; // recalculation of the vertices
	for(unsigned int i=0; i<this->M; i++) {
		buf = this->p[i];
		this->p[i].x = m2*buf.x + l2*buf.y + n2*buf.z;
		this->p[i].y = m1*buf.x + l1*buf.y + n1*buf.z;
		this->p[i].z = m3*buf.x + l3*buf.y + n3*buf.z;
	}
}

// reductive recalculation of coordinates of all crystal vertices
// for new Euler angles tt, ps and for fi = 0
void Crystal::ChangeVertices(double tt, double ps)
{
	#ifdef __FASTMATH_H
	double c1, c2, s1, s2;
	sincos(tt, &s1, &c1); sincos(ps, &s2, &c2);
	#else
	const double c1 = cos(tt), c2 = cos(ps),
				 s1 = sin(tt), s2 = sin(ps);
	#endif
	Point3D buf; 
	for(unsigned int i=0; i<this->M; i++) {
		buf = this->p[i];
		this->p[i].x = c1*(c2*buf.x-s2*buf.y)+s1*buf.z;
		this->p[i].y = s2*buf.x + c2*buf.y;
		this->p[i].z = s1*(s2*buf.y - c2*buf.x) + c1*buf.z; 
	}
}

// recalculation of the crystal facets
void Crystal::ChangeFacets(void)
{
	for(unsigned int i=0; i<this->K; i++) { // Fill out the array of coefficients of plane equation
		const int i0 = this->Gr[i][0], i1 = this->Gr[i][1], i2 = this->Gr[i][2];
		double* const q = this->Gran[i];
		q[0] = (this->p[i1].y-this->p[i0].y) * (this->p[i2].z-this->p[i0].z) -
				(this->p[i2].y-this->p[i0].y) * (this->p[i1].z-this->p[i0].z);
		q[1] = (this->p[i2].x-this->p[i0].x) * (this->p[i1].z-this->p[i0].z) -
				(this->p[i1].x-this->p[i0].x) * (this->p[i2].z-this->p[i0].z);
		q[2] = (this->p[i1].x-this->p[i0].x) * (this->p[i2].y-this->p[i0].y) -
				(this->p[i2].x-this->p[i0].x) * (this->p[i1].y-this->p[i0].y);
		q[3] = - q[0]*this->p[i0].x - q[1]*this->p[i0].y - q[2]*this->p[i0].z;
		double s = sqrt(SQR(q[0])+SQR(q[1])+SQR(q[2]));
		#ifdef _DEBUG
		if(fabs(s) <= DBL_EPSILON) throw " Error! ";
		#endif
		if(q[3]>0) s=-s;
		q[0] /= s; q[1] /= s; q[2] /= s; q[3] /= s;
	}
}

// returns true if the point lies on the facet number ii 
bool Crystal::PointInFacet(const Point3D& a, int ii) const
{
	double total = 0;
	const Point3D n = this->NormToFacet(ii);
	int fl = nMX(fabs(n.x),fabs(n.y),fabs(n.z));
	for(unsigned int j=0; this->Gr[ii][j+1]!=-1 && j<this->Km; j++) {
		Point2D pint(Proj(a,fl)),
				pd0(Proj(this->p[this->Gr[ii][j]],fl)),
				pd1(Proj(this->p[this->Gr[ii][j+1]],fl));
        double x = signedAngle(pint,pd1,pd0);
		if(x == M_PI) return false;
		total += x;
	}
	return (fabs(fabs(total)-m_2pi)<1e-5 ? true : false);
}
//------------------------------------------------------------------------------
