
#ifndef CRYSTAL_H
#define CRYSTAL_H

#ifndef GEOMETRY_HPP
#include "Geometry.hpp"
#endif

#ifndef BEAM_HPP
#include "Beam.hpp"
#endif

#ifndef ScatteringH
//#include "Scattering.hpp"
#endif

/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */

/** @addtogroup particles Crystal particles
 * @{
 */

/** 
@brief Function prototype witch should be used for handling the outgoing beam in main program.
@return none
*/
typedef void (*Hand) (Beam&); 
//==============================================================================

	/**
	@brief This class is a prototype for Crystal class. It contains parameters for BeamSplitting algorithm.
	*/
class BeamSplitting{
	void AllocMem(void);
	void FreeMem(void);
protected:
	unsigned int 	Itr, ///<Number of internal reflection
					NoF; ///<Number of facets

	Point3D k,  ///< Direction to the incident wave
			Ey; ///< Normal to the incident plane
    double Ep, ///< the level of negligibility energy
            d, ///< distance to far zone
            S_eps;	///<Minimal area of a beam
public:
	unsigned int **Facets; ///<Mask for calculation only specified trajectories
	BeamSplitting(unsigned int _Itr, unsigned int _NoF, Point3D _k, Point3D _Ey) :
        Itr(_Itr), NoF(_NoF), k(_k), Ey(_Ey), Ep(2e-12), d(10000.0), S_eps(0.01)
		{ this->AllocMem();};
	~BeamSplitting (void) 
		{ this->FreeMem();};
};

	/**
	@brief This class is a prototype for Crystal class.
	*/
class Frame {
	void AllocMem(void); ///< allocation of memory
	void FreeMem(void); ///< deallocation of the memory
protected:
	/**
	@brief Numbers of the vertexes from which the facet consist of.

	For example: GR[0][1]=5, means that the vertex number 1 of the facet number 0 is defined in p[5]
	*/
	int** Gr; 
    unsigned int	M, ///< Total count of the crystal vertexes
					K, ///< Total count of the facets
					Km; ///< Maximal count of vertexes at any facet
	/**
	@brief The coefficients of an equation of plane which contains the facet of the crystal

	For example: Gran[0]={A,B,C,D}, means that the facet 0 lies in a plane with equation Ax+By+Cz+D=0
	*/
	double** Gran;
    Point3D* p; ///< Coordinates of vertexes of the crystal
	Frame(int m, int k, int km) : M(m), K(k), Km(km)
		{ this->AllocMem(); }
	Frame(const Frame& f)
		{ this->M=f.M; this->K=f.K; this->Km=f.Km; this->AllocMem(); }
public:
	// members
	Frame operator=(const Frame&);
    ~Frame(void) ///< destructor
		{ this->FreeMem(); } 
	/// Returns a parameter of the size of the crystal, depends on the realization
	virtual double First(void) { throw " Frame::First(): Error! "; } 
	/// Returns a parameter of the size of the crystal, depends on the realization
	virtual double Second(void) { throw " Frame::Second(): Error! "; }
	/// Returns a parameter of the size of the crystal, depends on the realization
	virtual double Third(void) { throw " Frame::Third(): Error! "; }
	/// Returns a parameter of the size of the crystal, depends on the realization
	virtual double Forth(void) { throw " Frame::Forth(): Error! "; }
};
//------------------------------------------------------------------------------



/** 
@brief This class contains all information about particle's geometry and the parameters of tracing.
*/
class Crystal : public Frame, public BeamSplitting {
	/// This function checks if the point \b a lies on the facet number \b ii.
    bool PointInFacet(const Point3D& a, int ii) const;
	/// Beam splitting algorithm for internal beam. The realization is in scattering.cpp
    int TracingOfInternalBeam(Beam &, unsigned int, Hand) const;
protected:
	/// This flag specifies whether the optical path is taking into account or not
	bool PhaseDelay; // taking into account the optical path
	double 	mi,		/// Imaginary part of refraction index
			mr;		///Real part of refraction index
	complex RefI;	/// The complex value of refraction index of the particle
	Crystal(const complex& r, unsigned int m, unsigned int k, unsigned int km, unsigned int Itr, Point3D _k, Point3D _Ey) :
		Frame(m,k,km),
        BeamSplitting(Itr,k, _k, _Ey),
		PhaseDelay(false),
		mi(4.0*SQR(real(r)*imag(r))),
		mr(SQR(real(r))-SQR(imag(r))), RefI(r) {}
	/// This function defines the coordinates of the vertexes.
	virtual void SetVertices(void)
	{ throw " Crystal::SetVertices(): Error! "; }
	// reductive recalculation of coordinates of all crystal vertices
	// for new Euler angles tt, ps and for fi = 0
	/// The function rotates the vertex of particle by Euler's angles: Betta, Gamma. (Alpha = 0)
	void ChangeVertices(double bt, double gm);
	/// The function rotates the vertex of particle by Euler's angles: Betta, Gamma, Alpha.
	void ChangeVertices(double bt, double gm , double al);
	/// Sets up facet's normals
	void ChangeFacets(void);    
public:
    virtual ~Crystal(void) {}
	// members
	/// The function sets whether the optical path is taking into account or not
	bool& Phase(void) 
		{ return this->PhaseDelay; }
	/// The function returns if the optical path is taking into account or not
	bool Phase(void) const 
		{ return this->PhaseDelay; }
	/// The function returns complex value of the refraction index
	complex RefIn(void) const 
		{ return this->RefI; }
	/// The function returns a normal to the face number i
    Point3D NormToFacet(unsigned int i) const // i - facet number
        { return Point3D(this->Gran[i][0], this->Gran[i][1], this->Gran[i][2]); }
	/// The function returns coordinates of vertex number j on the face number i 
    Point3D VertexOfFacet(unsigned int i, unsigned int j) const
		{ return this->p[this->Gr[i][j]]; } // i - facet number, j - vertex number
	/// The function returns a pointer on an array with coefficients of an equation of plane which contains the facet number \b i of the crystal 
	const double* Facet(unsigned int i) const 
		{ return this->Gran[i]; }
	/// The function returns a projection of a facet number \b i to the one of 2D planes: YZ (k=0), XZ (k=1), XY (k=2).
	Polyg Projection(unsigned int i, int k) const;
	// The function returns a list of vertex of Polygon in 3D by it's \b ProjectionPolygon, which has been projected from facet number i to on of 2D planes: YZ (k=0), XZ (k=1), XY (k=2)
    std::list<Point3D> Retrieve(Polyg& ProjectionPolygon, unsigned int i, int k) const;
	/// The function rotates the particle by Euler's angles: Betta, Gamma, Alpha.
	virtual void ChangePosition(double t, double p, double f) //bt gm alph
		{ this->SetVertices(); this->ChangeVertices(t,p,f); this->ChangeFacets(); }
	/// The function returns the center of gravity of the facet number \b i
    Point3D CenterOfFacet(unsigned int i) const; // center of facet
	/// The function returns the area of the facet number \b i
	double AreaOfFacet(unsigned int i) const;
	/**
	@brief This function runs the beam splitting algorithm. The realization is in scattering.cpp

	@param Handler: Function to handle outgoing beam
	@return Scattering Cross-Section
	*/
    double FTforConvexCrystal(Hand Handler) const;
	/// The function calculates an intersection between a beam \b Bm and the facet number \b i. The result of intersection is in \b pl. 
	///@return \b true, if there is an intersection
    bool Intersection(const Beam& Bm, unsigned int i, std::list<Point3D> & pl) const;
};
//------------------------------------------------------------------------------

/** @}
*/

/** @}
*/

#endif
