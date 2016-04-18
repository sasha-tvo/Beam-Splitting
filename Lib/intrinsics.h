#pragma once

#include <xmmintrin.h>
#include <smmintrin.h>
#include <vector>

#include "types.h"

#define EPS_IN_LINE 0.001
#define EPS_IN_POLYGON 0.001

/**
 * @brief The PointPosition struct
 * Position of point on the facet (uses in 'inFacet' function)
 */
struct PointPosition {
	int position; /// 1 - inside, 0 - on the side or vertex, -1 - outside
	int facet_side_index_1; /// index of side that point belongs
	int facet_side_index_2; /// index of second side (point belongs to vertex between the sides)
};

/**
 * @brief inLine_i Is point belongs to line
 * @param x Point
 * @param a First point of line
 * @param b Second point of line
 * @return True if belongs
 */
bool inLine_i(const Point3 &x, const Point3 &a, const Point3 &b)
{
	__m128 _x = _mm_load_ps(x.point);
	__m128 _a = _mm_load_ps(a.point);
	__m128 _b = _mm_load_ps(b.point);

	__m128 ab = _mm_sub_ps(_a, _b);
	__m128 ax = _mm_sub_ps(_a, _x);
	__m128 bx = _mm_sub_ps(_b, _x);

	__m128 sqr_len_ab = _mm_dp_ps(ab, ab, 0x71);
	__m128 sqr_len_ax = _mm_dp_ps(ax, ax, 0x71);
	__m128 sqr_len_bx = _mm_dp_ps(bx, bx, 0x71);

	return (sqr_len_ax[0] + sqr_len_bx[0] < sqr_len_ab[0] + EPS_IN_LINE);
}

/** Dot product */
//__m128 res = _mm_dp_ps(_mm_load_ps(p1.point), _mm_load_ps(p2.point), 0x71);

/** Cross product */
inline __m128 crossProduct_i(__m128 a, __m128 b)
{
	return _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1)),
								 _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2))),
					  _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2)),
								 _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1))));
}

/**
 * @brief intersection_i Intersection of vector and plane
 * @param source_point Start point of source vector
 * @param source_vector Source vector
 * @param d_param D-param of plane
 * @param normal_vector Normal of plane
 * @return Intersection point
 */
__m128 intersection_i(const Point3 &source_point, const Point3 &source_vector,
					  float d_param, const Point3 &normal_vector)
{
	__m128 sp = _mm_load_ps(source_point.point);
	__m128 sv = _mm_load_ps(source_vector.point);
	__m128 nv = _mm_load_ps(normal_vector.point);
	__m128 d = _mm_load_ps1(&d_param);

	__m128 dp0 = _mm_dp_ps(sv, nv, 0x7F);
	/// TODO: проверить случай, когда нормаль и вектор перпендикулярны (dp0 == 0)

	__m128 dp1 = _mm_dp_ps(sp, nv, 0x7F);

	__m128 add = _mm_add_ps(dp1, d);
	__m128 t = _mm_div_ps(add, dp0);

	__m128 mul = _mm_mul_ps(t, sv);
	return _mm_sub_ps(sp, mul);
}

/**
 * @brief inFacet_i Position of point on the facet
 * @param x Point
 * @param facet Facet
 * @param pos Position of point (result)
 */
void inFacet_i(const Point3 &x, const Facet &facet, PointPosition &pos)
{
	pos.facet_side_index_1 = -1;
	pos.facet_side_index_2 = -1;
	pos.position = 1;

	int size = facet.size;

	__m128 _x = _mm_load_ps(x.point);
	__m128 _n = _mm_load_ps(facet.normal.point);

	__m128 dir;
	__m128 p1, p2;

	__m128 eps = _mm_set_ss(EPS_IN_POLYGON);
	__m128 m_eps = _mm_set_ss(-EPS_IN_POLYGON);

	p2 = _mm_load_ps(facet.vertices[size-1].point);

	for (int i = 0; i < size; ++i)
	{
		p1 = p2;
		p2 = _mm_load_ps(facet.vertices[i].point);

		__m128 a = _mm_sub_ps(p2, p1);
		__m128 b = _mm_sub_ps(_x, p1);

		dir = _mm_dp_ps(_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1)),
											  _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2))),
								   _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2)),
											  _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1)))), _n, 0x71);

		if (_mm_ucomilt_ss(dir, m_eps))
		{
			pos.position = -1;
			return;
		}
		else if (_mm_ucomilt_ss(dir, eps))
		{
			if (pos.position != 0)
			{
				pos.position = 0;
				pos.facet_side_index_1 = i-1;
			}
			else {
				pos.facet_side_index_2 = i-1;
			}
		}
	}
}

__m128 intersection_line_facet_i(const Point3 &a, const Point3 &b, const Facet &facet)
{
	int size = facet.size;

	__m128 p1, p2;
	__m128 n = _mm_load_ps(facet.normal.point);

	__m128 _a = _mm_load_ps(a.point);
	__m128 _b = _mm_load_ps(b.point);
	__m128 source_vector = _mm_sub_ps(_b, _a);

	__m128 eps = _mm_set_ss(EPS_IN_LINE);

	p2 = _mm_load_ps(facet.vertices[size-1].point);

	for (int i = 0; i < size; ++i)
	{
		p1 = p2;
		p2 = _mm_load_ps(facet.vertices[i].point);

		__m128 p12 = _mm_sub_ps(p2, p1);

		/// normal of new plane (cross product)
		__m128 n2 = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(p12, p12, _MM_SHUFFLE(3,0,2,1)),
										  _mm_shuffle_ps(n, n, _MM_SHUFFLE(3,1,0,2))),
							   _mm_mul_ps(_mm_shuffle_ps(p12, p12, _MM_SHUFFLE(3,1,0,2)),
										  _mm_shuffle_ps(n, n, _MM_SHUFFLE(3,0,2,1))));

		/// normalize normal
		__m128 nn2 = _mm_mul_ps(n2, _mm_rsqrt_ps(_mm_dp_ps(n2, n2, 0x7F)));

		/// D
		__m128 m_d_param = _mm_dp_ps(p1, nn2, 0x7F);

		/// intersection vector and new plane
		__m128 dp0 = _mm_dp_ps(source_vector, nn2, 0x7F);
		/// TODO: проверить случай, когда нормаль и вектор перпендикулярны (dp0 == 0)

		__m128 dp1 = _mm_dp_ps(_a, nn2, 0x7F);

		__m128 add = _mm_sub_ps(dp1, m_d_param);
		__m128 t = _mm_div_ps(add, dp0);

		__m128 mul = _mm_mul_ps(t, source_vector);
		__m128 p = _mm_sub_ps(_a, mul);

		/// is point in line (p1, p2)
		__m128 ab = _mm_sub_ps(p1, p2);
		__m128 ax = _mm_sub_ps(p1, p);
		__m128 bx = _mm_sub_ps(p2, p);

		__m128 sqr_len_ab = _mm_dp_ps(ab, ab, 0x71);
		__m128 sqr_len_ax = _mm_dp_ps(ax, ax, 0x71);
		__m128 sqr_len_bx = _mm_dp_ps(bx, bx, 0x71);

		if (_mm_add_ss(sqr_len_ax, sqr_len_bx)[0] < _mm_add_ss(sqr_len_ab, eps)[0])
		{
			/// is point in line (a, b)
			ab = _mm_sub_ps(_a, _b);
			ax = _mm_sub_ps(_a, p);
			bx = _mm_sub_ps(_b, p);

			sqr_len_ab = _mm_dp_ps(ab, ab, 0x71);
			sqr_len_ax = _mm_dp_ps(ax, ax, 0x71);
			sqr_len_bx = _mm_dp_ps(bx, bx, 0x71);

			if (_mm_add_ss(sqr_len_ax, sqr_len_bx)[0] < _mm_add_ss(sqr_len_ab, eps)[0])
			{
				return p;
			}
		}
	}
}
