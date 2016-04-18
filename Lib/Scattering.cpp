#ifndef CRYSTAL_H
#include "Crystal.hpp"
#endif

#include "Intersection.h"
#include <QDebug>

struct WorkBeam {
	unsigned int facet_index;
	unsigned int level;
	Beam beam;
};

//const double Eps1 = 1.7453292431333680334067268304459e-4; //cos(89.99)
const double Eps1 = 1.7453292519943295769148298069306e-10; //cos(89.99999999)
const double Eps2 = 0.99999999998254670756866631966593; //1- cos(89.99999999)

const unsigned MAX_BEAM_VERTEX_NUM = 64;

using namespace std;

inline bool Crystal::projection(const Beam& inputBeam, __m128 *_output_points, __m128 _normal, unsigned int &facetIndex) const
{
	__m128 _direction = _mm_setr_ps(inputBeam.r.x, inputBeam.r.y, inputBeam.r.z, 0.0);

	__m128 _d_param = _mm_set_ps1((float)(Gran[facetIndex][3]));
	__m128 _dp0 = _mm_dp_ps(_direction, _normal, 0x7F);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_PROJECTION) {
		return false; /// beam is parallel to facet
	}

	int i = 0;

	for (auto it = inputBeam.v.begin(); it != inputBeam.v.end(); ++it)
	{
		__m128 point = _mm_setr_ps((*it).x, (*it).y, (*it).z, 0.0);
		__m128 dp1 = _mm_dp_ps(point, _normal, 0x7F);
		__m128 add = _mm_add_ps(dp1, _d_param);
		__m128 t = _mm_div_ps(add, _dp0);
		__m128 mul = _mm_mul_ps(t, _direction);

		_output_points[i++] = _mm_sub_ps(point, mul);
	}

	return true;
}

inline int Crystal::calcFacetSize(unsigned int facetIndex) const
{
	int facetSize = 0;

	while (Gr[facetIndex][facetSize] != -1) {
		++facetSize;
	}

	return --facetSize;
}

inline void Crystal::setOutputBeam(__m128 *_output_points, int outputSize, std::list<Point3D>& outputBeam) const
{
	Point3D p;

	__m128 eps = _mm_load_ps1(&EPS_INTERSECTION);
	__m128 sign_mask = _mm_set1_ps(-0.f);

	__m128 p0 = _output_points[outputSize-1];

	for (int i = 0; i < outputSize; ++i)
	{
		__m128 difference = _mm_sub_ps(_output_points[i], p0);
		__m128 abs = _mm_andnot_ps(sign_mask, difference);
		__m128 cmp = _mm_cmplt_ps(eps, abs);

		int res = _mm_movemask_ps(cmp) & 0b111;

		if (res != 0) {
			p.x = _output_points[i][0];
			p.y = _output_points[i][1];
			p.z = _output_points[i][2];
			outputBeam.push_back(p);
		}

		p0 = _output_points[i];
	}
}

bool Crystal::Intersection2(const Beam& inputBeam, unsigned int facetIndex, std::list<Point3D>& outputBeam) const
{
	__m128 _output_points[MAX_BEAM_VERTEX_NUM];
	int outputSize = inputBeam.v.size();

	__m128 _normal_to_facet = _mm_setr_ps(Gran[facetIndex][0], Gran[facetIndex][1], Gran[facetIndex][2], 0.0);

	bool isProjected = projection(inputBeam, _output_points, _normal_to_facet, facetIndex);

	if (!isProjected) {
		return false;
	}

	__m128 _buffer[MAX_BEAM_VERTEX_NUM];
	int bufferSize;

	int facetSize = calcFacetSize(facetIndex);

	__m128 _p1, _p2; /// points of facet
	__m128 _s_point, _e_point; /// points of projection
	bool isInsideE, isInsideS;

	Point3D lastPoint = p[Gr[facetIndex][facetSize-1]];
	_p2 = _mm_setr_ps(lastPoint.x, lastPoint.y, lastPoint.z, 0.0);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;

		Point3D nextPoint = p[Gr[facetIndex][i]];
		_p2 = _mm_setr_ps(nextPoint.x, nextPoint.y, nextPoint.z, 0.0);

		bufferSize = outputSize;
		outputSize = 0;

		for (int j = 0; j < bufferSize; ++j)
		{
			_buffer[j] = _output_points[j];
		}

		_s_point = _buffer[bufferSize-1];
		isInsideS = is_inside_i(_s_point, _p1, _p2, _normal_to_facet);

		bool isIntersected;

		for (int j = 0; j < bufferSize; ++j)
		{
			_e_point = _buffer[j];
			isInsideE = is_inside_i(_e_point, _p1, _p2, _normal_to_facet);

			if (isInsideE)
			{
				if (!isInsideS)
				{
					__m128 x = computeIntersection_i(_s_point, _e_point, _p1, _p2, _normal_to_facet, isIntersected);

					if (isIntersected)
					{
						_output_points[outputSize++] = x;
					}
				}

				_output_points[outputSize++] = _e_point;
			}
			else if (isInsideS)
			{
				__m128 x = computeIntersection_i(_s_point, _e_point, _p1, _p2, _normal_to_facet, isIntersected);

				if (isIntersected)
				{
					_output_points[outputSize++] = x;
				}
			}

			_s_point = _e_point;
			isInsideS = isInsideE;
		}
	}

	setOutputBeam(_output_points, outputSize, outputBeam);
	return outputBeam.size() > 2;
}

//==============================================================================
// face-tracing code
double Crystal::FTforConvexCrystal(Hand Handler) const
{
	unsigned int cnt = 0;
	#ifdef _DEBUG
		static unsigned int mncnt = 0;
		mncnt++;
	#endif

	// start beam
	Beam i_beam;
	i_beam.isInternal = false;

	double square = 0, s = 0, csa, csta; // the area of the particle projection
	unsigned int *faces_0;
	faces_0 = Facets[0];

	WorkBeam beams[128];
	unsigned int dept = 0;

	dept ++;
	beams[dept-1].facet_index = -1;
	beams[dept-1].level = cnt;
	beams[dept-1].beam = /*std::move*/(i_beam);

	Beam incedence_beam;
	unsigned int facet;
	unsigned int *faces_cnt;
	std::list<Point3D> Result;
	std::list<unsigned int>::const_iterator fs;
	std::list<Point3D>::const_iterator p;
	Point3D Eh, n, pnt, pcnt, pv, R, T;

	const complex Tref = 2.0*RefI/(1.0+RefI);
	const complex Rref = (1.0-RefI)/(1.0+RefI);

	do { // tree
		facet = beams[dept-1].facet_index;
		cnt = beams[dept-1].level;
		incedence_beam = beams[dept-1].beam;
		dept --;

		if (!incedence_beam.isInternal)
		{
			for (unsigned int l = 0; l< K; l++)
			{ // facet loop
				if(!faces_0[l]) continue;
				const Point3D n = this->NormToFacet(l);
				csa = k*n;

				if(csa < Eps1)
					continue;

				Beam ext_beam, int_beam;  //
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
				square += AreaOfBeam(ext_beam)*csa; // equal to call of CrossSection(ext_beam)
				ext_beam.N = n;
				Handler(ext_beam.Rotate(k,Ey)); // handling of the outgoing beam

				beams[dept].facet_index = l;
				beams[dept].level = cnt;
				beams[dept].beam = int_beam;
				dept ++;
			}
		}
		else { // internal beam
			if (norm(incedence_beam()) < Ep) {
				continue;
			}
			if (cnt >= Itr) {
				continue;
			}

			cnt++;
			faces_cnt = Facets[cnt];

			for (unsigned int l = 0; l < K; l ++)
			{ // facet loop
				if(!faces_cnt[l]) {
					continue; // if this facet is out of interest
				}
				if(l == facet) {
					continue;
				}

				n = std::move(NormToFacet(l));
				csa = n*incedence_beam.r;

				if(csa < Eps1) {
					continue;
				}

				//--------------------------------------------------------------------------
				Result.clear();

//				if(!Intersection(incedence_beam, l, Result)) {
//					continue;
//				}
				if(!Intersection2(incedence_beam, l, Result)) {
					continue;
				}

				Beam Reflected(Result), // internal (reflected) beam
						Refracted(Result); // outgoing beam

				for (fs = incedence_beam.BeginP(); fs != incedence_beam.EndP(); fs++) {
					Reflected.PushBackP(*fs);
					Refracted.PushBackP(*fs);
				}

				Reflected.PushBackP(l);
				Refracted.PushBackP(l);

				double sqr_csa = SQR(csa),
						Nr = (mr+sqrt(SQR(mr)+mi/sqr_csa))/2.0;
				//--------------------------------------------------------------------------
				if(csa >= Eps2) { // case of the normal incidence
					Refracted(0,0) = Tref*incedence_beam(0,0);
					Refracted(0,1) = Tref*incedence_beam(0,1);
					Refracted(1,0) = Tref*incedence_beam(1,0);
					Refracted(1,1) = Tref*incedence_beam(1,1);
					Refracted.e = incedence_beam.e;
					Refracted.r = incedence_beam.r;
					//--------------------------------------------------------------------------
					Reflected(0,0) = Rref*incedence_beam(0,0);
					Reflected(0,1) = Rref*incedence_beam(0,1);
					Reflected(1,0) = -Rref*incedence_beam(1,0);
					Reflected(1,1) = -Rref*incedence_beam(1,1);
					Reflected.e =  incedence_beam.e;//
					Reflected.r = -incedence_beam.r;
					//------------------------------------------------------------------------
					if (Phase()) { // optical path
						pnt = std::move(CenterOfBeam(Reflected));
						Refracted.D = -pnt*Refracted.r;
						Reflected.D = -pnt*Reflected.r;
						Reflected.lng = incedence_beam.lng + sqrt(Nr)*fabs(incedence_beam.r*pnt+incedence_beam.D);
						Refracted.lng = Reflected.lng + fabs(d+Refracted.D);
					}
					//--------------------------------------------------------------------
					p = Refracted.Begin();
					pcnt = std::move(CenterOfBeam(Refracted));
					pv = std::move(((*++p)-pcnt)%((*p)-pcnt));

					if (pv*Refracted.r < 0.0)
					{
						Refracted.Clear();
						p = Result.begin();

						for(unsigned int j=0; j<Result.size(); p++, j++) {
							Refracted.PushFront(*p);
						}
					}

					//--------------------------------------------------------------------
					Refracted.N = n;
					Handler(Refracted.Rotate(k, Ey)); // treating of the outgoing beam
				}
				else { // case of the sloping incidence
					R = std::move(incedence_beam.r/csa-n);
					T = std::move(R-n);
					T /= length(T);
					s = 1.0/(Nr*sqr_csa)-norm(R);
					Eh = std::move(n%incedence_beam.r); // perpendicular to the scattering plane
					Eh /= length(Eh);
					Beam bm = incedence_beam.RotatePlane(Eh);
					Reflected.e = Eh;
					Reflected.r = T;

					const complex RefI_csa = RefI*csa;

					if(s > DBL_EPSILON)
					{
						R = std::move(R/sqrt(s)+n);
						R /= length(R); // direction of the escaping ray
						csta = n*R;   // refraction angle

						const complex bf = 2.0*RefI_csa,
								bTv = RefI*csta + csa,
								bTh = RefI_csa + csta;

						complex Tv = bf/bTv,
								Th = bf/bTh;
						// outgoing beam
						Refracted(0,0) = Tv*bm(0,0);
						Refracted(0,1) = Tv*bm(0,1);
						Refracted(1,0) = Th*bm(1,0);
						Refracted(1,1) = Th*bm(1,1);
						Refracted.e = Eh;
						Refracted.r = R;
						//----------------------------------------------------------------------
						if(Phase()) { // optical path
							pnt = std::move(CenterOfBeam(Reflected));
							Reflected.D = -pnt*T;
							Refracted.D = -pnt*R;
							Reflected.lng = incedence_beam.lng + sqrt(Nr)*fabs(incedence_beam.r*pnt+incedence_beam.D);
							Refracted.lng = Reflected.lng + fabs(d+Refracted.D);
						}
						//----------------------------------------------------------------------
						p = Refracted.Begin();
						pcnt = std::move(CenterOfBeam(Refracted));
						pv = std::move(((*++p)-pcnt)%((*p)-pcnt));

						if (pv*Refracted.r < 0.0) {
							Refracted.Clear();
							p = Result.begin();
							for(unsigned int j=0; j < Result.size(); p++, j++)
								Refracted.PushFront(*p);
						} //

						Refracted.N = n;
						Handler(Refracted.Rotate(k,Ey)); // treating of the escaping beam

						if (false /** Вышедший пучок посылаем в стек, а не в Handler */) {
							beams[dept].facet_index = l;
							beams[dept].level = cnt;
							beams[dept].beam = Reflected;
							dept ++;
						}

						Tv = (csa-RefI*csta) / bTv;
						Th = (RefI_csa-csta) / bTh;   // internal field
						Reflected(0,0) = Tv*bm(0,0);
						Reflected(0,1) = Tv*bm(0,1);
						Reflected(1,0) = Th*bm(1,0);
						Reflected(1,1) = Th*bm(1,1);
					}
					else { // case of the complete internal reflection
						const double bf = Nr*(1.0-sqr_csa)-1.0;
	#ifdef _DEBUG
						if(bf < 0.0) throw "TracingOfInternalBeam : Error!";
	#endif
						const complex sq(0, bf>0?sqrt(bf):0);
						complex Rv = (csa-RefI*sq)/(RefI*sq+csa),
								Rh = (RefI_csa-sq)/(RefI_csa+sq);
						Reflected(0,0) = Rv*bm(0,0);
						Reflected(0,1) = Rv*bm(0,1);// reflected field
						Reflected(1,0) = Rh*bm(1,0);
						Reflected(1,1) = Rh*bm(1,1);

						if(Phase()) { // optical path
							pnt = std::move(CenterOfBeam(Reflected));
							Reflected.D = -pnt*Reflected.r;
							Reflected.lng = incedence_beam.lng + sqrt(Nr)*fabs(incedence_beam.r*pnt+incedence_beam.D);
						}
					} // end of if... else...
				}
				beams[dept].facet_index = l;
				beams[dept].level = cnt;
				beams[dept].beam = Reflected;
				dept ++;
			} // end of loop
		}
	}
	while (dept != 0);

	return square;
}
//==============================================================================
