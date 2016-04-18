#pragma once

#include "intrinsics.h"
#include <QVector>

const double EPS_PROJECTION = 0.00174532836589830883577820272085;
const float EPS_INTERSECTION = 0.01;

bool inside(const Point3 &x, const Point3 &p1, const Point3 &p2, const Point3 &normal)
{
	__m128 _x = _mm_load_ps(x.point);
	__m128 _p1 = _mm_load_ps(p1.point);
	__m128 _p2 = _mm_load_ps(p2.point);
	__m128 _n = _mm_load_ps(normal.point);

	__m128 eps = _mm_set_ss(EPS_INTERSECTION);
	__m128 m_eps = _mm_set_ss(-EPS_INTERSECTION);

	__m128 a = _mm_sub_ps(_p2, _p1);
	__m128 b = _mm_sub_ps(_x, _p1);

	__m128 dir = _mm_dp_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1)),
												 _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2))),
									  _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2)),
												 _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1)))),
						   _n, 0x71);

	return _mm_ucomigt_ss(dir, m_eps) && _mm_ucomilt_ss(dir, eps);
}

inline bool is_inside_i(__m128 x, __m128 p1, __m128 p2, __m128 normal)
{
	__m128 m_eps = _mm_set_ss(-EPS_INTERSECTION);

	__m128 p1_p2 = _mm_sub_ps(p2, p1);
	__m128 p1_x = _mm_sub_ps(x, p1);

	__m128 dir = _mm_dp_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(p1_p2, p1_p2, _MM_SHUFFLE(3,0,2,1)),
												 _mm_shuffle_ps(p1_x, p1_x, _MM_SHUFFLE(3,1,0,2))),
									  _mm_mul_ps(_mm_shuffle_ps(p1_p2, p1_p2, _MM_SHUFFLE(3,1,0,2)),
												 _mm_shuffle_ps(p1_x, p1_x, _MM_SHUFFLE(3,0,2,1)))),
						   normal, 0x71);

	return _mm_ucomigt_ss(dir, m_eps);
}

inline __m128 computeIntersection_i(__m128 _a1, __m128 _a2, __m128 _b1, __m128 _b2, __m128 _normal_to_facet, bool &ok)
{
	__m128 _length_a = _mm_sub_ps(_a2, _a1);
	__m128 _length_b = _mm_sub_ps(_b2, _b1);

	/// normal of new plane (cross product)
	__m128 _normal_to_line = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(_length_b, _length_b, _MM_SHUFFLE(3,0,2,1)),
												   _mm_shuffle_ps(_normal_to_facet, _normal_to_facet, _MM_SHUFFLE(3,1,0,2))),
										_mm_mul_ps(_mm_shuffle_ps(_length_b, _length_b, _MM_SHUFFLE(3,1,0,2)),
												   _mm_shuffle_ps(_normal_to_facet, _normal_to_facet, _MM_SHUFFLE(3,0,2,1))));
	/// normalize normal
	__m128 _normal_n = _mm_mul_ps(_normal_to_line, _mm_rsqrt_ps(_mm_dp_ps(_normal_to_line, _normal_to_line, 0x7F)));

	/// intersection vector and new plane
	__m128 _dp0 = _mm_dp_ps(_length_a, _normal_n, 0x7F);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_INTERSECTION)
	{
		ok = false;
		return _dp0;
	}

	__m128 _dp1 = _mm_dp_ps(_a1, _normal_n, 0x7F);
	__m128 m_d_param = _mm_dp_ps(_b1, _normal_n, 0x7F);

	__m128 _add = _mm_sub_ps(_dp1, m_d_param);
	__m128 _t = _mm_div_ps(_add, _dp0);

	__m128 _mul = _mm_mul_ps(_t, _length_a);

	ok = true;
	return _mm_sub_ps(_a1, _mul);
}

void computeIntersection(const Point3 &s, const Point3 &e,
						 const Point3 &p1, const Point3 &p2, const Point3 &normal,
						 Point3 &x)
{
	__m128 _a = _mm_load_ps(s.point);
	__m128 _b = _mm_load_ps(e.point);

	__m128 _p1 = _mm_load_ps(p1.point);
	__m128 _p2 = _mm_load_ps(p2.point);

	__m128 n = _mm_load_ps(normal.point);

	__m128 ab = _mm_sub_ps(_b, _a);
	__m128 p12 = _mm_sub_ps(_p2, _p1);

	/// normal of new plane (cross product)
	__m128 n2 = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(p12, p12, _MM_SHUFFLE(3,0,2,1)),
									  _mm_shuffle_ps(n, n, _MM_SHUFFLE(3,1,0,2))),
						   _mm_mul_ps(_mm_shuffle_ps(p12, p12, _MM_SHUFFLE(3,1,0,2)),
									  _mm_shuffle_ps(n, n, _MM_SHUFFLE(3,0,2,1))));

	/// normalize normal
	__m128 nn2 = _mm_mul_ps(n2, _mm_rsqrt_ps(_mm_dp_ps(n2, n2, 0x7F)));

	/// D
	__m128 m_d_param = _mm_dp_ps(_p1, nn2, 0x7F);

	/// intersection vector and new plane
	__m128 dp0 = _mm_dp_ps(ab, nn2, 0x7F);
	__m128 dp1 = _mm_dp_ps(_a, nn2, 0x7F);

	__m128 add = _mm_sub_ps(dp1, m_d_param);
	__m128 t = _mm_div_ps(add, dp0);

	__m128 mul = _mm_mul_ps(t, ab);
	__m128 p = _mm_sub_ps(_a, mul);

	x.point[0] = p[0];
	x.point[1] = p[1];
	x.point[2] = p[3];
}

void intersection_concave(const Facet &face, const Polygon &beam, Polygon &out_beam)
{
	QVector<Point3> inputList, outputList;

	for (int i = 0; i < beam.size; ++i)
	{
		outputList.append(Point3());
		outputList[i].point[0] = beam.points[i].point[0];
		outputList[i].point[1] = beam.points[i].point[1];
		outputList[i].point[2] = beam.points[i].point[2];
	}

	Point3 p1, p2;
	p2 = face.vertices[face.size-1];

	bool e_inside, s_inside;

	for (int i = 0; i < face.size; ++i)
	{
		p1 = p2;
		p2 = face.vertices[i];

		inputList = outputList;
		outputList.clear();

		Point3 s = inputList.last();
		s_inside = inside(s, p1, p2, face.normal);

		for (int j = 0; j < inputList.count(); ++j)
		{
			Point3 e = inputList[j];

			e_inside = inside(e, p1, p2, face.normal);

			if (e_inside)
			{
				if (!s_inside)
				{
					Point3 x;
					computeIntersection(s, e, p1, p2, face.normal, x);
					outputList.append(x);
				}

				outputList.append(e);
			}
			else if (s_inside)
			{
				Point3 x;
				computeIntersection(s, e, p1, p2, face.normal, x);
				outputList.append(x);
			}

			s.point[0] = e.point[0];
			s.point[1] = e.point[1];
			s.point[2] = e.point[2];
		}
	}

	out_beam.size = outputList.size();

	for (int i = 0; i < outputList.size(); ++i)
	{
		out_beam.points[i].point[0] = outputList[i].point[0];
		out_beam.points[i].point[1] = outputList[i].point[1];
		out_beam.points[i].point[2] = outputList[i].point[2];
	}
}

