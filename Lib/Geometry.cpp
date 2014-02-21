//==============================================================================
#include "Geometry.hpp"


//==============================================================================
/// Position of the 2D point *this relatively the edge {org, dest}.
/// Originally it was written by M.Laszlo (p.76) and was refined a bit by us for our
/// beam-splitting code
position Point2D::classify(const Point2D& org, const Point2D& dst) const
{
	const Point2D a  = dst-org, b = *this-org;
	const double  la = length(a), lb = length(b),
				  sa = a.x*b.y - b.x*a.y; // = a*normalTo(b);
	// the order of the tests is highly IMPORTANT, please, don't change it!
	if(fabs(sa) > E_DBL)                                  return sa>0 ? LEFT : RIGHT;
	else
		if(lb < E_DBL*la)                                 return ORIGIN;
		else
			if(length(dst-*this) < E_DBL*la)              return DESTINATION;
			else
				if(a.x*b.x<-E_DBL || a.y*b.y<-E_DBL)      return BEHIND;
				else
					if(la < lb)                           return BEYOND;
					else                                  return BETWEEN;
}

//==============================================================================
/// distance between the edge e and the point with sign
double  Point2D::distance(const Edge &e) const
{
	Edge ab = e;
	ab.flip().rot90();
	Point2D n(ab.dest-ab.org);
	n /= length(n);
	Edge f(*this, *this+n);
	double t;
	f.intersect(e,t);
	return t;
}
//==============================================================================
position  Point2D::classify(const Edge &e) const 
{ 
	return this->classify(e.org, e.dest); 
}
//==============================================================================	
double    Point2D::dist(const Edge &e) const // unsigned distance
{ 
	const Point2D de = e.dest-e.org;
	return fabs(de.x*(this->y-e.dest.y)-de.y*(this->x-e.dest.x))/length(de); 
}

//==============================================================================
/// returns polar angle (radians) of the 2D vector a
double  polarAngle(const Point2D& a)
{
	if(!a.x)
	{
		if(!a.y) return -1.0;
		else     return a.y>0.0 ? M_PI_2 : 1.5*M_PI;
	}
	const double tt = atan(a.y/a.x);
	if(a.x>0) return (a.y>=0 ? tt : m_2pi + tt);
	else      return (tt+M_PI);
}
//==============================================================================
/// returns signed angle

double signedAngle(const Point2D& a, const Point2D &org, const Point2D &dst)
{
	const double va = polarAngle(org - a), wa = polarAngle(dst - a);
	if(va == -1.0 || wa == -1.0) return M_PI;
	const double x = wa - va;
	if(fabs(x) == M_PI)    return M_PI;
	else if(x < -M_PI)     return x+m_2pi;
		 else if(x > M_PI) return x-m_2pi;
			  else         return x;
}

double  Edge::signedAngle(const Point2D& a) const
   { return ::signedAngle(a, this->org, this->dest); }


//==============================================================================
/// Rotation of the point around the 3D axis \b a through the angle t < 17.10.2002 >.

/// see G.A.Korn, T.M.Korn, MATHEMATICAL HANDBOOK FOR SCIENTISTS AND ENGINEERS
/// - 5th edition, 1984, p. 448 (Russian)
Point3D  Point3D::rotate(double t, const Point3D& a) const
{
	t /= 2.0;
	const double sn = sin(t), ro = cos(t);
	// Euler parameters: lam, mu, nu
	const double lam = a.x*sn,        mu = a.y*sn,       nu = a.z*sn,
				 LamMu = lam*mu,      LamNu = lam*nu,    LamRo = lam*ro,
				 RoMu = ro*mu,        RoNu = ro*nu,      MuNu = mu*nu,
				 sqr_lam = SQR(lam),  sqr_mu = SQR(mu),
				 sqr_nu = SQR(nu),    sqr_ro = SQR(ro);
	// elements of the rotation matrix a[3X3]
	const double a11 = sqr_lam-sqr_mu-sqr_nu+sqr_ro,
				 a12 = LamMu-RoNu,
				 a13 = LamNu+RoMu,
				 a21 = LamMu+RoNu,
				 a22 = sqr_mu-sqr_nu-sqr_lam+sqr_ro,
				 a23 = MuNu-LamRo,
				 a31 = LamNu-RoMu,
				 a32 = MuNu+LamRo,
				 a33 = sqr_nu-sqr_lam-sqr_mu+sqr_ro;
	return Point3D(this->x*a11 + 2.0*(this->y*a12+this->z*a13),
				   this->y*a22 + 2.0*(this->x*a21+this->z*a23),
				   this->z*a33 + 2.0*(this->x*a31+this->y*a32));
}
//==============================================================================
Point3D  Point3D::rotateAxisY(double a) const
{
	Point3D b = *this;
	double sn = sin(a), cs = cos(a);
	return Point3D(b.x*cs-b.z*sn, b.y, b.x*sn+b.z*cs);
}
//==============================================================================
Point3D  Point3D::rotateAxisZ(double a) const
{
	Point3D b = *this;
	double sn = sin(a), cs = cos(a);
	return Point3D(b.x*cs-b.y*sn, b.x*sn+b.y*cs, b.z);
}
//==============================================================================
Point3D  Point3D::RotateEuler(double tt, double ps, double f)
{
	const double c1 = cos(tt), c2 = cos(ps), c3 = cos(f),
				 s1 = sin(tt), s2 = sin(ps), s3 = sin(f);
	const double  bf1 = c1*s2, bf2 = c1*c2,
				  l1 =  c2*c3-bf1*s3,  l2 = -bf1*c3-c2*s3,  l3 =  s1*s2,
				  m1 =  s2*c3+bf2*s3,  m2 =  bf2*c3-s2*s3,  m3 = -s1*c2,
				  n1 =  s1*s3,         n2 =  s1*c3,         n3 =  c1;
	Point3D buf = *this;
	return Point3D(m2*buf.x + l2*buf.y + n2*buf.z,
				   m1*buf.x + l1*buf.y + n1*buf.z,
				   m3*buf.x + l3*buf.y + n3*buf.z);
}


//==============================================================================
// members of the class "Edge"
/// intersection of 2D edges
state  Edge::intersect(const Edge &e, double &t) const
{
	const Point2D a = this->org, b = this->dest,
				  n = Point2D(e.dest.y-e.org.y, e.org.x-e.dest.x); // normal to e
	Point2D tmp = b-a;
	const double denom = n*tmp;
	if(fabs(denom) < E_DBL*length(n)*length(tmp)) {
		const position aclass = this->org.classify(e.org, e.dest);
		if(aclass==LEFT || aclass==RIGHT) return PARALLEL;
		else return COLLINEAR;
	}
	tmp = a-e.org;
	const double num = n*tmp;
	t = -num/denom;
	return SKEW;
}
// friend of the class "Edge"
direction  parallel(const Edge& e, const Edge& f)
{
	Point2D a = e.dest - e.org,
			b = f.dest - f.org;
	const double la = length(a), lb = length(b);
	#ifdef _DEBUG
	if(la < E_DBL || lb < E_DBL)
		throw "parallel : Incorrect arguments!";
	#endif
	a /= la; b /= lb;
	const double s = a*b;
	if(s < E_DBL-1.0)     return OPPOSITE;
	if(s > 1.0-E_DBL)     return ACCORDANT;
	if(fabs(s) < E_DBL) return PERPENDICULAR;
	return NON_PARALLEL;
}
//==============================================================================
// members of the class "Node"
Node*  Node::insert(Node *b)
{
	Node *c = this->_next;
	b->_next = c;
	b->_prev = this;
	this->_next = b;
	c->_prev = b;
	return b;
}
Node*  Node::remove(void)
{
	this->_prev->_next = this->_next;
	this->_next->_prev = this->_prev;
	this->_next = this->_prev = this;
	return this;
}
void  Node::splice(Node *b)
{
	Node *a = this,
		 *an = a->_next,
		 *bn = b->_next;
	a->_next = bn;
	b->_next = an;
	an->_prev = b;
	bn->_prev = a;
}
//==============================================================================
// members of the class "Vertex"
Vertex*  Vertex::split(Vertex *b)
{
	Vertex *bp = b->ccw()->insert(new Vertex(b->point())); // placing bp before b
	this->insert(new Vertex(this->point())); // placing ap after current vertex see Laszlo p. 86
	this->splice(bp);
	return bp;
}
//==============================================================================
// members of the class "Polyg"
void  Polyg::resize(void)
{
	if(!this->_v) this->_size=0;
	else {
		Vertex* v = this->_v->cw();
		for(this->_size=1; v!=this->_v; ++this->_size, v=v->cw());
	}
}
// copy-function for the class (copy constructor)
void  Polyg::Copy(const Polyg& p)
{
	this->_size = p._size;
	if(!this->_size) this->_v = NULL;
	else {
		this->_v = new Vertex(p.point());
		for(unsigned int i=1; i<this->_size; i++) {
			p.advance(CLOCKWISE);
			this->_v = this->_v->insert(new Vertex(p.point()));
		}
		p.advance(CLOCKWISE);
		this->_v = this->_v->cw();
	}
}
// destructor of the "Polyg"
 Polyg::~Polyg(void) 
{
	if(this->_v) {
		Vertex *w = this->_v->cw();
		while(this->_v != w) {
			delete w->remove();
			w = this->_v->cw();
		}
		delete this->_v;
	}
}
//----
bool  Polyg::Test(void)
{
	if(this->size()<3)  return true;
	const double Eps = 1e-5;
	for(unsigned int i=0; i<this->size(); i++, this->advance(CLOCKWISE)) {
		const Point2D& pt1 = this->point(),
					   pt2 = this->neighbor(CLOCKWISE)->point();
		if(length(pt1-pt2) < Eps) this->remove();
	}
	return this->size()<3 ? true : false; //
}
// friends of the "Polyg"
double  AreaOfConvexPolygon(const Polyg &p)
{
	double s = 0;
	const Point2D& bgn = p.point();
	p.advance(CLOCKWISE);
	unsigned int i = 1;
	Point2D p1 = p.point() - bgn;
	while(i < p.size()) {
		i++;
		p.advance(CLOCKWISE);
		const Point2D p2 = p.point() - bgn;
		s += fabs(p1.x*p2.y - p1.y*p2.x);
		p1 = p2;
	}
	return s/2.0;
}
//------
void  OutPolyg(char* name, char* comment, const Polyg &p)
{
	std::ofstream out(name, std::ios::app);
	out.precision(20);
	out << comment << ':' << std::endl;
	for(unsigned int i=0; i<p.size(); i++, p.advance(CLOCKWISE)) {
		const Point2D& pt = p.point();
		out << "a[" << i << "] = Point2D(" << pt.x << "," << pt.y << ");\n" << std::flush;
	}
}
//------
bool  CheckPolygon(const Polyg &p)
{
	for(unsigned int i=0; i<p.size(); i++, p.advance(CLOCKWISE)) {
		const Point2D& pt1 = p.point(),
					   pt2 = p.neighbor(CLOCKWISE)->point();
		const double ln = length(pt1-pt2);
		if(ln < E_DBL) return false;
	}
	return true;
}
//---------------------------------------------------------------------------
double  Polyg::PolygonIntersection (const Polyg& A, const Polyg& B)
{
	Point2D max_A(A.point()), min_A(A.point()),
			max_B(B.point()), min_B(B.point());
	for(unsigned int i=0; i<A.size(); i++, A.advance(CLOCKWISE)) {
		const Point2D& p = A.point();
		if(p.x > max_A.x) max_A.x = p.x;  if(p.y > max_A.y) max_A.y = p.y;
		if(p.x < min_A.x) min_A.x = p.x;  if(p.y < min_A.y) min_A.y = p.y;
	}
	for(unsigned int i=0; i<B.size(); i++, B.advance(CLOCKWISE)) {
		const Point2D& p = B.point();
		if(p.x > max_B.x) max_B.x = p.x;  if(p.y > max_B.y) max_B.y = p.y;
		if(p.x < min_B.x) min_B.x = p.x;  if(p.y < min_B.y) min_B.y = p.y;
	}
	if(max_A.x<=min_B.x || max_B.x<=min_A.x ||
		max_A.y<=min_B.y || max_B.y<=min_A.y) return 0;

	#ifdef _DEBUG
	if(!CheckPolygon(A)) {
		OutPolyg("Errors.log", "A", A);
		throw "Crystal::Intersection : Error!";
	}
	if(!CheckPolygon(B)) {
		OutPolyg("Errors.log", "B", B);
		throw "Crystal::Intersection : Error!";
	}
	#endif

	Polyg Result;
	if(FastPolygonIntersect(A, B, Result)) {
		#ifdef _DEBUG
		if(!CheckPolygon(Result)) {
			OutPolyg("Errors.log", "B", B);
			OutPolyg("Errors.log", "A", A);
			OutPolyg("Errors.log", "Result", Result);
			throw "Crystal::Intersection : Error!";
		}
		#endif
		const double e = AreaOfConvexPolygon(Result);
		return e;
	}
	return 0;
}

//==============================================================================
/// Checks if the 2D point \b s belongs to the 2D convex polygon \b p
bool  pointInConvexPolygon(const Point2D &s, const Polyg &p)
{
	if(p.size() == 1)
		return (s==p.point());
	if(p.size() == 2) {
		const position c = s.classify(p.edge().org, p.edge().dest);
		return c==BETWEEN || c==ORIGIN || c==DESTINATION;
	}
	Vertex* const org = p.v();
	for(unsigned int i=0; i<p.size(); i++, p.advance(CLOCKWISE))
		if(s.classify(p.edge().org, p.edge().dest) == LEFT) {
			const_cast<Polyg&> (p).setV(org);
			return false;
		}
	return true; //
}
//==============================================================================
// checks intersection of 2D edges \b e and \b f
state  crossingPoint(const Edge& e, const Edge& f, Point2D& p)
{
	// original version of the "crossingPoint" of M. Laszlo
/*  double s, t;
	const state classe = e.intersect(f,s);
	if(classe == COLLINEAR || classe == PARALLEL) return classe;
	const double lene = length(e);
	if(s < -E_DBL*lene || s > 1.0+E_DBL*lene) return SKEW_NO_CROSS;
	f.intersect(e,t);
	const double lenf = length(f);
	if((t >= -E_DBL*lenf) && (t <= 1.0+E_DBL*lenf)) {
		if(t <= E_DBL*lenf) p = f.org;
		else if(t >= 1.0-E_DBL*lenf) p = f.dest;
				else if(s <= E_DBL*lene) p = e.org;
							else if(s >= 1.0-E_DBL*lene) p = e.dest;
									else p = f.point(t);
		return SKEW_CROSS;
	}
	else return SKEW_NO_CROSS; */

	const Point2D de = e.dest-e.org,    df = f.dest-f.org;
	const double  le = length(de),      lf = length(df), lelf = le*lf;

	Point2D n = Point2D(df.y, -df.x); // normal to f
	const double denom = n*de;
	if(fabs(de*df) > lelf-E_DBL*lelf) {
		const position aclass = f.org.classify(e.org, e.dest);
		return aclass < 2 ? PARALLEL : COLLINEAR;
	} 
	Point2D tmp = e.org - f.org;
	if(length(tmp) < E_DBL) { p = f.org; return SKEW_CROSS;}

	const double s = -n*tmp/denom;
	if(s < -E_DBL || s > 1.0+E_DBL) return SKEW_NO_CROSS;

	n = Point2D(de.y, -de.x); // normal to e
	const double t = n*tmp/(n*df);

	if(t >= -E_DBL && t <= 1.0+E_DBL) {
		if(t <= E_DBL) p = f.org;
		else if(t >= 1.0-E_DBL) p = f.dest;
			 else if(s <= E_DBL) p = e.org;
				  else if(s >= 1.0-E_DBL) p = e.dest;
					   else p = f.point(t);
		return SKEW_CROSS;
	}
	else return SKEW_NO_CROSS; 
}
//==============================================================================
/// this function returns TRUE if and only if edge \b a aims at edge \b b. See p. 156
bool aimsAt(const Edge& a, const Edge& b, position aclass, state crossType)
{
	const Point2D va = a.dest - a.org,
				  vb = b.dest - b.org;
	return  (crossType != COLLINEAR) ?
			(va.x*vb.y>=vb.x*va.y ? aclass!=RIGHT : aclass!=LEFT) :
			(aclass != BEYOND); 
}
//==============================================================================
/// moving to the next edge (auxiliary routine for convexPolygonIntersect)
void  advance(const Polyg& A, Polyg& R, bool inside)
{
	A.advance(CLOCKWISE);
	if(inside)
		if(length(R.point()-A.point()) > E_DBL) R.insert(A.point());
} //
//==============================================================================
// "convexPolygonIntersect" < 19.10.2002 >
/// The routine of Laszlo was slightly refined for degenerated cases. It returns
/// true if the intersection R of the input polygons P and Q is not empty.
bool  convexPolygonIntersect(const Polyg &P, const Polyg& Q, Polyg& R)
{
	Point2D  iPnt, startPnt;
	crossing inflag = UNKNOWN;
	bool     phase = true,
			 flag = true;
	const    unsigned int maxItns = 2*(P.size()+Q.size());
	//----------------------------------------------------------------------------
	for(unsigned int i=1; i<=maxItns || !phase; i++) {
		#ifdef _DEBUG   // the safe return
		if(i>1000) {
			OutPolyg("P.dat", "beam", P);
			OutPolyg("Q.dat", "facet", Q);
			//throw "convexPolygonIntersect::Error!";
			return false;
		}
		#endif
		const Edge p = P.edge(), q = Q.edge();
		const position  pclass = p.dest.classify(q.org, q.dest),
						qclass = q.dest.classify(p.org, p.dest);
		const state crossType = crossingPoint(p, q, iPnt);
		switch(crossType) { // kind of the intersection
			case SKEW_CROSS: {
				if(phase) {
					phase = false;
					R.insert(startPnt=iPnt);
				}
				else {
					if(iPnt != R.point()) {
						if(iPnt != startPnt)  R.insert(iPnt);
						else return true;
						}
					else if(flag) flag = false; // for the case of intersection of the
						 else if(iPnt == startPnt) // polygons P and Q in one bound point
							  return true;
				}
				inflag = (pclass == RIGHT) ? P_IS_INSIDE :
						 (qclass == RIGHT  ? Q_IS_INSIDE : UNKNOWN);
			}
			break;
			// for the case of the intersection of the polygons by bound edge
			case COLLINEAR: {
				if(parallel(p,q) == OPPOSITE) {
					#ifdef _DEBUG
					  if(R.size() > 2)  throw "convexPolygonIntersect::Error!";
					#endif
					if(R.size() == 2) return true;
					if(pclass==BETWEEN && p.org.classify(q.org, q.dest)==BETWEEN) {
						if(!R.size()) R.insert(p.org);
						R.insert(p.dest);
						return true;
					}
					if(qclass==BETWEEN && q.org.classify(p.org, p.dest)==BETWEEN) {
						if(!R.size()) R.insert(q.org);
						R.insert(q.dest);
						return true;
					}
					if(pclass == qclass) {
						switch(pclass) {
							case BETWEEN:
								if(!R.size()) { R.insert(q.dest); R.insert(p.dest);}
								else if(R.point() != p.dest) R.insert(p.dest);
										else R.insert(q.dest);
								return true;
							case BEHIND:
								if(!R.size()) { R.insert(q.org); R.insert(p.org);}
								else if(R.point() != p.org) R.insert(p.org);
										else R.insert(q.org);
								return true;
							default : {
								#ifdef _DEBUG
								  throw "convexPolygonIntersect::Error!";
								#endif
							};
						}
					}
				}
				if(pclass!=BEHIND && qclass!=BEHIND)  inflag = UNKNOWN;
			}
			break;
			default : {
				#ifdef _DEBUG
				  throw "convexPolygonIntersect::Error!";
				#endif
			};
		}
		bool pAIMSq = aimsAt(p, q, pclass, crossType),
			 qAIMSp = aimsAt(q, p, qclass, crossType);
		if(pAIMSq && qAIMSp) {
			if((inflag == Q_IS_INSIDE) || (inflag==UNKNOWN && pclass==LEFT))
				advance(P, R, false);
			else advance(Q, R, false);
		}
		else
			if(pAIMSq)      { advance(P, R, inflag == P_IS_INSIDE);}
			else if(qAIMSp) { advance(Q, R, inflag == Q_IS_INSIDE);}
				 else {
					if((inflag == Q_IS_INSIDE) || (inflag==UNKNOWN && pclass==LEFT))
						advance(P, R, false);
					else advance(Q, R, false);
				 }
	} // end of loop
	//----------------------------------------------------------------------------
	if(pointInConvexPolygon(P.point(),Q))       { R = P; return true;}
	else if(pointInConvexPolygon(Q.point(),P))  { R = Q; return true;}
	return false;
}
inline position base_clsfy(const Point2D& d, const Point2D& n, double eps, bool cmp)
{
	return (d.x*n.x<-eps||d.y*n.y<-eps) ? BEHIND : (cmp?BEYOND:BETWEEN);
}
//==============================================================================
// eps > 0.0
position clsfy(const Point2D& d, const Point2D& n, double eps, double ld, double ln)
{
	const double e = d.x*n.y - d.y*n.x;
	//----------------------------------------------------------------------------
	if((ld*ln) < eps) { // relatively short edge
		const Point2D _d = d/ld, _n = n/ln;
		const double  e1 = length(_d-_n);
		if(e1>eps && e1<2.0-eps) return e>0 ? LEFT : RIGHT;
		else return base_clsfy(d,n,eps,ld<ln);
	}
	else // relatively long edge
		return fabs(e)>eps ? (e>0 ? LEFT : RIGHT) : base_clsfy(d,n,eps,ld<ln);
}

//==============================================================================
/// Checks if the 2D point \b s lies inside the 2D convex polygon \b p
/// The routine originally was written by M. Laszlo and slightly refined for our
/// beam-splitting code;
bool  InclusionInsideTest(const Point2D &s, const Polyg &p)
{
	#ifdef _DEBUG
	if(p.size()<3)  throw "InclusionInsideTest : Error!";
	#endif
	Vertex* const org = p.v();
	for(unsigned int i=0; i<p.size(); i++, p.advance(CLOCKWISE)) {
		const Edge& e = p.edge();
		const Point2D de = e.dest - e.org,
					  ds = s - e.org;
        if(de.x*ds.y > de.y*ds.x+E_DBL) {
			const_cast<Polyg&> (p).setV(org);
			return false;
		} 
	}
	return true;
}
//==============================================================================
/// It is just the same as previous "convexPolygonIntersect".
///
/// The routine was specially optimized for beam-splitting algorithm. It returns
/// true if and only if the intersection R of the input convex polygons P and Q
/// is not singular polyhedron.
bool  FastPolygonIntersect(const Polyg &P, const Polyg& Q, Polyg& R)
{
	const     unsigned int  maxItns = 2*(P.size()+Q.size());
    const     double   Eps = pow(0.1,10), E_Size = 0.005;
	Point2D   iPnt,    startPnt;
	crossing  inflag = UNKNOWN;
	bool      phase = true;
	//----------------------------------------------------------------------------
	for(unsigned int i=1; i<=maxItns || !phase; i++) {
		state crossType;
		position pclass, qclass;
		if(i>100) {
			#ifdef _DEBUG
				OutPolyg("P.dat", "beam", P);
				OutPolyg("Q.dat", "facet", Q);
				throw "FastPolygonIntersect : Error!";
			#endif
			return false;
		}

		const Edge&     p = P.edge(),         q = Q.edge();
		const Point2D   dp = p.dest-p.org,    dq = q.dest-q.org; // vectors p and q
		const double    lp = length(dp),      lq = length(dq);
		if(lp < E_Size) { P.advance(CLOCKWISE); continue;}
		if(lq < E_Size) { Q.advance(CLOCKWISE); continue;}
		const double    l = MAX(lp,lq),
						lplq = lp*lq,
						E = Eps*lplq, //MAX(Eps*lplq, Eps), //MAX(E_DBL*lplq, E_DBL),
						dpdq = dp*dq;
		const double    E2 = Eps*l; //MAX(Eps, Eps*l);
		const bool      opp = dpdq < E-lplq, // p and q are opposite
						acc = dpdq > lplq-E, // p and q are accordant
						par = opp || acc; // p and q are parallel
		//--------------------------------------------------------------------------
		if(length(q.dest-p.dest) < E2) { // q.dest == p.dest
			pclass = qclass = DESTINATION;
			crossType = par ? COLLINEAR : (iPnt = q.dest, SKEW_CROSS);
		}
		else {
			const Point2D q_p = q.dest-p.org,  p_q = p.dest-q.org;
			const double  lq_p = length(q_p),  lp_q = length(p_q);
			if(lq_p<E2 && lp_q<E2)
				return false;
			if(lq_p < E2) { // q.dest == p.org
				qclass = ORIGIN;
				pclass = clsfy(dq, p_q, E, lq, lp_q);
				crossType = par ? COLLINEAR : (iPnt = q.dest, SKEW_CROSS);
			}
			else {
				if(lp_q < E2) { // p.dest == q.org
						pclass = ORIGIN;
						qclass = clsfy(dp, q_p, E, lp, lq_p);
						crossType = par ? COLLINEAR : (iPnt = p.dest, SKEW_CROSS);
				}
				else {
					if(par) { // the case of the parallel edges
						const double pp = p.dest.dist(q),
									 qq = q.dest.dist(p);
						if(pp < E2 || qq < E2) { // p and q are collinear
							if(opp) return false;
							pclass = base_clsfy(dq, p_q, E, lq<lp_q);
							qclass = base_clsfy(dp, q_p, E, lp<lq_p);
							crossType = COLLINEAR;
						}
						else { // p and q are parallel
							pclass = clsfy(dq, p_q, E, lq, lp_q);
							qclass = clsfy(dp, q_p, E, lp, lq_p);
							if(pclass == LEFT && qclass == LEFT && opp)
								return false;
							crossType = PARALLEL;
						}
					}
					else { // the case of the non-parallel
						pclass = clsfy(dq, p_q, E, lq, lp_q);
						qclass = clsfy(dp, q_p, E, lp, lq_p);
						const Point2D tmp = p.org - q.org;  // p.org == q.org
						if(length(tmp) < E2) { crossType = SKEW_CROSS; iPnt = q.org;}
						else { // intersection point
							double E0 = Eps;
							Point2D n = Point2D(dq.y, -dq.x); // normal to q
							const double denom = n*dp,
										 s = -n*tmp/denom;
							if(s<-E0 || s>1.0+E0) crossType = SKEW_NO_CROSS;
							else {
								n = Point2D(dp.y, -dp.x); // normal to e
								const double t = n*tmp/(n*dq);
								if(t>=-E0 && t<=1.0+E0) {
									crossType = SKEW_CROSS;
									if(t <= E0) iPnt = q.org;
									else 
										if(t >= 1.0-E0) iPnt = q.dest;
										else 
											if(s <= E0) iPnt = p.org;
											else 
												if(s >= 1.0-E0) iPnt = p.dest;
												else iPnt = q.point(t); 
								}
								else crossType = SKEW_NO_CROSS; 
							}
						}
					}
				}
			}
		}
		bool pAIMSq, qAIMSp;
		switch(crossType) { // kind of the intersection
			case COLLINEAR: {
				if(opp)
					return false; // if p and q are opposite
				if(!phase) {
					if(pclass == BETWEEN) {
						inflag = P_IS_INSIDE; P.advance(CLOCKWISE); continue;
					}
					if(qclass == BETWEEN) {
						inflag = Q_IS_INSIDE; Q.advance(CLOCKWISE); continue;
					}
				}
				if(pclass!=BEHIND && qclass!=BEHIND)  inflag = UNKNOWN;
				pAIMSq = pclass!=BEYOND;
				qAIMSp = qclass!=BEYOND;
			} break;
			case PARALLEL: {
				if(!phase) {
					if(acc) // p and q are accordant
						inflag = (pclass == RIGHT) ? P_IS_INSIDE :
								 (qclass == RIGHT  ? Q_IS_INSIDE : UNKNOWN);
					else
						if(opp) inflag = UNKNOWN;
				}
				pAIMSq = pclass!=RIGHT;
				qAIMSp = qclass!=RIGHT;
			} break;
			case SKEW_CROSS: {
				if(phase) {
					phase = false;
					R.insert(startPnt=iPnt);
				}
				else {
					double s = length(iPnt - R.point());
					if(s > E2) {
						s = length(iPnt - startPnt);
						if(s > E2) R.insert(iPnt);
						else
							return R.size()<3 ? false : true;
					}
					else 
						if(i>maxItns)
							return R.size()<3 ? false : true;
				}
				if(pclass==ORIGIN)
					{ P.advance(CLOCKWISE); inflag=UNKNOWN; continue;}
				else 
					if(qclass==ORIGIN)
						{ Q.advance(CLOCKWISE); inflag=UNKNOWN; continue;}
					else 
						if(pclass==RIGHT)
							inflag = (qclass==BETWEEN?UNKNOWN:P_IS_INSIDE);
						else 
							if(qclass==RIGHT)
								inflag = (pclass==BETWEEN?UNKNOWN:Q_IS_INSIDE);
							else inflag = UNKNOWN; 
			} // There mustn't be "break"
			default: { // cases SKEW_NO_CROSS and SKEW_CROSS
				const bool  pNEQr = pclass!=RIGHT,
							qNEQr = qclass!=RIGHT;
				if(qclass == BEYOND) {
					switch(inflag) {
						case UNKNOWN:     (pNEQr ? P : Q).advance(CLOCKWISE);  break;
						case Q_IS_INSIDE: P.advance(CLOCKWISE);                break;
						case P_IS_INSIDE: {
							P.advance(CLOCKWISE);
							const Point2D& pp = P.point();
							if(length(pp-startPnt) < E2)
								return R.size()<3 ? false : true;
							if(length(pp-R.point()) > E2)  R.insert(pp);
							break;
						}
					}
					continue;
				}
				if(pclass == BEYOND) {
					switch(inflag) {
						case UNKNOWN:     (qNEQr ? Q : P).advance(CLOCKWISE);  break;
						case P_IS_INSIDE: Q.advance(CLOCKWISE);                break;
						case Q_IS_INSIDE: {
							Q.advance(CLOCKWISE);
							if(inflag == Q_IS_INSIDE) {
								const Point2D& qq = Q.point();
								if(length(qq-startPnt) < E2)
									return R.size()<3 ? false : true;
								if(length(qq-R.point()) > E2) R.insert(qq);
							}
							break;
						}
					}
					continue;
				} //*/
				const bool c = dp.x*dq.y-dq.x*dp.y >= 0; // dp*normTo(dq)
				pAIMSq =  c ? pNEQr        : pclass!=LEFT;
				qAIMSp =  c ? qclass!=LEFT : qNEQr;
			}
		}
		//--------------------------------------------------------------------------
		if(pAIMSq && qAIMSp) {
			if((inflag==Q_IS_INSIDE) || (inflag==UNKNOWN && pclass==LEFT))
				P.advance(CLOCKWISE); // advance(P, R, false, E2);
			else Q.advance(CLOCKWISE); // advance(Q, R, false, E2);
		}
		else
			if(pAIMSq) { //  advance(P, R, inflag == P_IS_INSIDE, E2);}
				P.advance(CLOCKWISE);
				if(inflag == P_IS_INSIDE) {
					const Point2D& pp = P.point();
					if(length(pp-startPnt) < E2)
						return R.size()<3 ? false : true;
					if(length(pp-R.point()) > E2)  R.insert(pp);
				}
			}
			else 
				if(qAIMSp) { // advance(Q, R, inflag == Q_IS_INSIDE, E2);}
					Q.advance(CLOCKWISE);
					if(inflag == Q_IS_INSIDE) {
						const Point2D& qq = Q.point();
						if(length(qq-startPnt) < E2)
							return R.size()<3 ? false : true;
						if(length(qq-R.point()) > E2) R.insert(qq);
					}
				}
				else {
					if((inflag == Q_IS_INSIDE) || (inflag==UNKNOWN && pclass==LEFT))
						P.advance(CLOCKWISE);     // advance(P, R, false, E2);
					else Q.advance(CLOCKWISE);  // advance(Q, R, false, E2);
				}
	} // end of loop
	//----------------------------------------------------------------------------
	if(InclusionInsideTest(P.point(),Q))
		{ R = P; return true;}
	else 
		if(InclusionInsideTest(Q.point(),P))
			{ R = Q; return true;}
	return false;
}
//==============================================================================
// removing from 2D polygon s all parts outside edge e
bool  clipPolygonToEdge(const Polyg &s, const Edge &e, Polyg& result)
{
	Polyg p;
	Point2D crossingPt;
	for(unsigned int i=0; i<s.size(); s.advance(CLOCKWISE), i++) {
		const Point2D org = s.point(),
					  dest = s.cw()->point();
		bool orgIsInside  = (org.classify(e.org, e.dest) != LEFT),
			 destIsInside = (dest.classify(e.org, e.dest) != LEFT);
		if(orgIsInside != destIsInside) {
			double t;
			e.intersect(s.edge(), t);
			crossingPt = e.point(t);
		}
		if(orgIsInside && destIsInside) p.insert(dest);
		else {
			if(orgIsInside && !destIsInside) {
				if(org != crossingPt) p.insert(crossingPt);
			}
			else {
				if(!orgIsInside && !destIsInside);
				else {
					p.insert(crossingPt);
					if(dest != crossingPt) p.insert(dest);
				}
			}
		}
	}
	result = p;
	return result.size()>0;
}
//==============================================================================
// the intersection of the polygons s and p; p must be convex
bool clipPolygon(const Polyg &s, const Polyg &p, Polyg& result)
{
	Polyg q(s), r;
	bool flag = true;
	for(unsigned int i=0; i<p.size(); i++, p.advance(CLOCKWISE)) {
		if(clipPolygonToEdge(q, p.edge(), r)) q = r;
		else { flag = false; break;}
	}
	return flag ? (result=q, true) : false;
}
//==============================================================================
Point3D  operator*(const Point3D& a, const matrix &mt)
{
	if(Str(mt)!=3 || Col(mt)!=3)  throw "operator*(Poin3d, matrix) : Error!";

	double  x = a.x*mt[0][0] + a.y*mt[1][0] + a.z*mt[2][0],
			y = a.x*mt[0][1] + a.y*mt[1][1] + a.z*mt[2][1],
			z = a.x*mt[0][2] + a.y*mt[1][2] + a.z*mt[2][2];
	return Point3D(x, y, z);
}

SphCrd   GetSpherical(const Point3D& r)
{
	if(fabs(r.z+1.0) < DBL_EPSILON) return SphCrd(0, M_PI); // forward
	if(fabs(r.z-1.0) < DBL_EPSILON) return SphCrd(0, 0);    // bacward
	double tmp = SQR(r.y);
	if(tmp < DBL_EPSILON) tmp = r.x>0?0:M_PI;
	else {
		tmp = acos(r.x/sqrt(SQR(r.x)+tmp));
		if(r.y<0) tmp = m_2pi-tmp;
	}
	return SphCrd(tmp<m_2pi?tmp:0, acos(r.z));
}
