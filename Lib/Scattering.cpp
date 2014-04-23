#ifndef CRYSTAL_H
#include "Crystal.hpp"
#endif


//const double Eps1 = 1.7453292431333680334067268304459e-4; //cos(89.99)
const double Eps1 = 1.7453292519943295769148298069306e-10; //cos(89.99999999)
const double Eps2 = 0.99999999998254670756866631966593; //1- cos(89.99999999)

//==============================================================================
// tracing of the separated beam bm from facet k INSIDE particle
int Crystal::TracingOfInternalBeam ( Beam& inc_beam, unsigned int i, Hand Handler) const
{
	static unsigned int cnt = 0;
	#ifdef _DEBUG
		static unsigned int mncnt = 0;
		mncnt++;
	#endif
	if(norm(inc_beam()) < Ep)  return 0;  
	if(cnt >= Itr)             return 1; // safe return
	cnt++;
	//----------------------------------------------------------------------------
	unsigned int *faces_cnt;
	faces_cnt = Facets[cnt];
	for(unsigned int l=0; l<this->K; l++) { // facet loop
		if(!faces_cnt[l]) continue; // if this facet is out of interest
		if(l == i) continue;
		const Point3D n = this->NormToFacet(l); 
		const double csa = n*inc_beam.r; 
		if(csa < Eps1) continue;
		//--------------------------------------------------------------------------
		std::list<Point3D> Result;
		if(!this->Intersection(inc_beam,l,Result)) continue;
		//--------------------------------------------------------------------------
		Beam Reflected(Result), // internal (reflected) beam
			 Refracted(Result); // outgoing beam
		std::list<unsigned int>::const_iterator fs = inc_beam.BeginP();
		for(;fs!=inc_beam.EndP(); fs++) {
			Reflected.PushBackP(*fs);
			Refracted.PushBackP(*fs);
		}
		Reflected.PushBackP(l);
		Refracted.PushBackP(l);

		Reflected.id = inc_beam.id;
		Refracted.id = inc_beam.id;
		Reflected.id *= (this->K+1); Reflected.id += (l+1);
		Refracted.id *= (this->K+1); Refracted.id += (l+1);
		double sqr_csa = SQR(csa),Nr = (this->mr+sqrt(SQR(this->mr)+this->mi/sqr_csa))/2.0;
		//--------------------------------------------------------------------------
		if(csa >= Eps2) { // case of the normal incidence
			const complex Tref = 2.0*this->RefI/(1.0+this->RefI);
			Refracted(0,0) = Tref*inc_beam(0,0);	Refracted(0,1) = Tref*inc_beam(0,1);
			Refracted(1,0) = Tref*inc_beam(1,0);	Refracted(1,1) = Tref*inc_beam(1,1);
			Refracted.e = inc_beam.e;
			Refracted.r = inc_beam.r;
			//--------------------------------------------------------------------------
			const complex Rref = (1.0-this->RefI)/(1.0+this->RefI); 
			Reflected(0,0) = Rref*inc_beam(0,0);	Reflected(0,1) = Rref*inc_beam(0,1);
			Reflected(1,0) = -Rref*inc_beam(1,0);	Reflected(1,1) = -Rref*inc_beam(1,1);
			Reflected.e =  inc_beam.e;//
			Reflected.r = -inc_beam.r;
			//------------------------------------------------------------------------
			if(this->Phase()) { // optical path
				const Point3D pnt = CenterOfBeam(Reflected);
				Refracted.D = -pnt*Refracted.r;
				Reflected.D = -pnt*Reflected.r;
				Reflected.lng =  inc_beam.lng + sqrt(Nr)*fabs(inc_beam.r*pnt+inc_beam.D);
				Refracted.lng = Reflected.lng + fabs(d+Refracted.D);
			}
			//--------------------------------------------------------------------
			std::list<Point3D>::const_iterator p = Refracted.Begin();
			Point3D p0 = *p, p1 = *++p, pcnt = CenterOfBeam(Refracted);
			Point3D pp0=p0-pcnt,pp1=p1-pcnt; 
			Point3D pv=pp1%pp0;
			if (pv*Refracted.r<0.0) {
				Refracted.Clear();
				p = Result.begin();
				for(unsigned int j=0; j<Result.size(); p++, j++)
					Refracted.PushFront(*p);
			} 
			//--------------------------------------------------------------------
			Refracted.N = n;
			Handler(Refracted.Rotate(k, Ey)); // treating of the outgoing beam
		}
		else { // case of the sloping incidence
			Point3D R = inc_beam.r/csa-n, T = R-n;
			T /= length(T);  
			double s = 1.0/(Nr*sqr_csa)-norm(R), csta;
			Point3D Eh = n%inc_beam.r; // perpendicular to the scattering plane
			Eh/=length(Eh);
			Beam bm = inc_beam.RotatePlane(Eh);
			Reflected.e = Eh;
			Reflected.r = T;
			if(s > DBL_EPSILON) { 
				R = R/sqrt(s)+n;
				R /= length(R); // direction of the escaping ray
				csta = n*R;   // refraction angle
				const complex 	bf = 2.0*this->RefI*csa,       
								bTv = this->RefI*csta + csa, 
								bTh = this->RefI*csa + csta; 
				
				complex Tv = bf/bTv, Th = bf/bTh;
				// outgoing beam
				Refracted(0,0) = Tv*bm(0,0);	Refracted(1,0) = Th*bm(1,0);
				Refracted(0,1) = Tv*bm(0,1);	Refracted(1,1) = Th*bm(1,1);
				Refracted.e = Eh;
				Refracted.r = R;
				//----------------------------------------------------------------------
				if(this->Phase()) { // optical path
					const Point3D pnt = CenterOfBeam(Reflected);
					Reflected.D = -pnt*T;
					Refracted.D = -pnt*R;
					Reflected.lng = inc_beam.lng + sqrt(Nr)*fabs(inc_beam.r*pnt+inc_beam.D);
					Refracted.lng = Reflected.lng + fabs(d+Refracted.D);
				}
				//----------------------------------------------------------------------
				std::list<Point3D>::const_iterator p = Refracted.Begin();
				Point3D p0 = *p, p1 = *++p, pcnt = CenterOfBeam(Refracted);
				Point3D pp0=p0-pcnt,pp1=p1-pcnt; 
				Point3D pv=pp1%pp0;
				if (pv*Refracted.r<0.0) {
					Refracted.Clear();
					p = Result.begin();
					for(unsigned int j=0; j<Result.size(); p++, j++)
						Refracted.PushFront(*p);
				} //
				//----------------------------------------------------------------------
				Refracted.N = n;
				Handler(Refracted.Rotate(k,Ey)); // treating of the escaping beam
				Tv = (csa-this->RefI*csta) / bTv;  
				Th = (this->RefI*csa-csta) / bTh;   // internal field
				Reflected(0,0) = Tv*bm(0,0);  Reflected(0,1) = Tv*bm(0,1);
				Reflected(1,0) = Th*bm(1,0);  Reflected(1,1) = Th*bm(1,1);
			}
			else { // case of the complete internal reflection
				const double bf = Nr*(1.0-sqr_csa)-1.0;
				#ifdef _DEBUG
					if(bf < 0.0) throw "TracingOfInternalBeam : Error!";
				#endif
				const complex sq(0, bf>0?sqrt(bf):0);
				complex Rv = (csa-this->RefI*sq)/(this->RefI*sq+csa), 
						Rh = (this->RefI*csa-sq)/(this->RefI*csa+sq); 
				Reflected(0,0) = Rv*bm(0,0);	Reflected(0,1) = Rv*bm(0,1);// reflected field
				Reflected(1,0) = Rh*bm(1,0);	Reflected(1,1) = Rh*bm(1,1);
				if(this->Phase()) { // optical path
					const Point3D pnt = CenterOfBeam(Reflected);
					Reflected.D = -pnt*Reflected.r;
					Reflected.lng = inc_beam.lng + sqrt(Nr)*fabs(inc_beam.r*pnt+inc_beam.D);
				}
			} // end of if... else...
		}
        //if(CrossSection(Reflected) < S_eps) continue;
		// treating of the internal field
		this->TracingOfInternalBeam(Reflected, l, Handler);
	} // end of loop
	cnt--;
	return 0;
}
//==============================================================================
// face-tracing code
// starting routine
double  Crystal::FTforConvexCrystal(Hand Handler) const
{
	double s = 0, csa, csta; // the area of the particle projection
	unsigned int *faces_0;
	faces_0 = Facets[0];
	for(unsigned int l=0; l<this->K; l++) { // facet loop
		if(!faces_0[l]) continue;
		const Point3D n = this->NormToFacet(l);
		csa = k*n;
		if(csa < Eps1) continue;
		Beam ext_beam, int_beam; //
		//--------------------------------------------------------------------------
		if(csa >= Eps2) { // case of the normal incidence
			int_beam(0,0) = int_beam(1,1) = 2.0/(this->RefI+1.0); // internal field
			int_beam.e = Ey;
			int_beam.r = -k;
			ext_beam(0,0)= (this->RefI-1.0)/(this->RefI+1.0);
			ext_beam(1,1)= -(this->RefI-1.0)/(this->RefI+1.0); 
			ext_beam.e = Ey;
			ext_beam.r = k;	//
		}
		//--------------------------------------------------------------------------
		else {
			double  sqr_csa = SQR(csa), 
					Nr = (this->mr+1.0-sqr_csa+sqrt(SQR(this->mr+sqr_csa-1.0)+this->mi))/2.0;
			Point3D T = n-k/csa, 
					R = n+T;
			R /= length(R); // direction of the external beam
			T = (T/sqrt((Nr/sqr_csa)-norm(T)))-n;
			T /= length(T);  // internal ray
			csta = -n*T;
			Point3D Eh = -n%k; 
			const double NormEh = norm(Eh);
			int_beam.e = Ey;
			int_beam.r = -k; 
			Eh /= sqrt(NormEh);
			int_beam = int_beam.RotatePlane(Eh); 
			complex bfTv = this->RefI*csa + csta,  
					bfTh = this->RefI*csta + csa,
					Tv = (this->RefI*csa-csta) / bfTv, 
					Th = (csa-this->RefI*csta) / bfTh; 

			ext_beam(0,0) = Tv*int_beam(0,0);	ext_beam(0,1) = Tv*int_beam(0,1);
			ext_beam(1,0) = Th*int_beam(1,0);	ext_beam(1,1) = Th*int_beam(1,1);
			ext_beam.r = R; // direction of the external beam
			ext_beam.e = int_beam.e; // direction of the external beam
			Tv = 2.0*csa/bfTv;		Th = 2.0*csa/bfTh;
			int_beam(0,0) *= Tv;	int_beam(0,1) *= Tv; // internal field
			int_beam(1,0) *= Th;	int_beam(1,1) *= Th; 
			int_beam.r = T; // direction of the internal beam
		}
		//--------------------------------------------------------------------------
		// vertexes of the external and internal beams
		for(unsigned int j=0; this->Gr[l][j+1]!=-1 && j<this->Km; j++) { // edge loop
			ext_beam.PushFront(this->p[this->Gr[l][j]]);
			int_beam.PushFront(this->p[this->Gr[l][j]]);
		}
		if(this->Phase()) {  // taking into account the optical path of the beams
			Point3D pnt = CenterOfBeam(int_beam);
			int_beam.D = -int_beam.r*pnt;
			ext_beam.D = -ext_beam.r*pnt;
			int_beam.lng = d-k*pnt;
			ext_beam.lng = int_beam.lng + fabs(d+ext_beam.D);
		}
		ext_beam.PushFrontP(l);
		int_beam.PushFrontP(l);
		ext_beam.id = l+1;
		int_beam.id = l+1;
		s += AreaOfBeam(ext_beam)*csa; // equal to call of CrossSection(ext_beam)
		ext_beam.N = n;
        Handler(ext_beam.Rotate(k,Ey)); // handling of the outgoing beam
		// handling of the internal beam
		this->TracingOfInternalBeam(int_beam, l, Handler);
	} // end of facet loop
	//----------------------------------------------------------------------------
	return s;
}
//==============================================================================
