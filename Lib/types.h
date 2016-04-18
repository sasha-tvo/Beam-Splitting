#pragma once

/**
 * @brief The Point3 struct
 * 3D coordinate point
 */
struct Point3
{
	float point[4]; /// coordinates

	Point3() {}
	Point3(float x, float y, float z) {
		point[0] = x;
		point[1] = y;
		point[2] = z;
	}
} __attribute__ ((aligned (16)));

/**
 * @brief The Facet struct
 * Facet of polyhedron
 */
struct Facet
{
	Point3 vertices[32];
	Point3 normal;
	int size;

	Facet() {}
	Facet(int p_size) {
		size = p_size;
	}

} __attribute__ ((aligned (16)));

struct Polygon
{
	Point3 points[32];
	int size;

	Polygon() {}
	Polygon(int p_size) {
		size = p_size;
	}

} __attribute__ ((aligned (16)));
