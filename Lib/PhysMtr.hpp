#ifndef PRHYSMTR_HPP
#define PRHYSMTR_HPP

#ifndef MATRIX_HPP
 #include "matrix.hpp"
#endif
//==============================================================================

/** @defgroup Auxiliary_classes Auxiliary classes
  * @{
  */

/**
 @brief The array with (N-rows x M-columns) dimensions of
 small \b real-value  matrixes with (n x m) dimensions.
 
 Example of using:
 @code
 ...
 Arr2D Arr(0,0,0,0); //creating an array
 main() 
 {
   ...
   Arr = Arr2D(1, 1, 2, 2); //changing the size
   Arr.ClearArr();   //clearing the array
   matrix mt(2,2);  
   Arr.insert(0,0, mt); //adding matrix mt to the array
   matrix M=Arr(0,0); //taking matrix M from the array
   ... 
 }
 @endcode
 
*/
class Arr2D {
private:
	unsigned int
		N, M,	///< N,M - Dimensions of array
		n, m;	///< n,m - Dimensions of matrixes
	double*** ptr;	///< Pointer to the array
	void AllocMem(double*** = NULL);
	void FreeMem(void);
public:
	// constructors
	///Creates the array with (_N-rows x _M-columns) dimensions of small real-value matrixes with (_n x _m) dimensions
	 Arr2D(unsigned int _N, unsigned int _M, unsigned int _n, unsigned int _m)
		{ this->N=_N; this->M=_M; this->n=_n; this->m=_m; this->AllocMem(); }
	 Arr2D(void)
		{ this->N=0; this->M=0; this->n=0; this->m=0; this->AllocMem(); }
	 Arr2D(const Arr2D &);
	// destructor
	 ~Arr2D(void) { this->FreeMem(); }
	// members
	matrix operator() (unsigned int, unsigned int) const; ///< Returns the matrix stored in the array by address (_N,_M)

	///Returns a reference to the matrix element (_n,_m) stored in the array by address (_N,_M)
	double& operator() (unsigned int _N, unsigned int _M, unsigned int _n, unsigned int _m)
		{ return this->ptr[_N][_M][_n*this->n+_m]; }

	/// Returns the matrix element (_n,_m) stored in the array by address (_N,_M)
	double operator() (unsigned int _N, unsigned int _M, unsigned int _n, unsigned int _m) const
		{ return this->ptr[_N][_M][_n*this->n+_m]; }
	Arr2D operator=(const Arr2D &);
	Arr2D operator+=(const Arr2D &);
	Arr2D operator*(double);

	/**
	@brief The function clears all elements of all matrixes in the array.		
	*/
	void ClearArr(void) const;

	/**
	@brief The function adds the matrix mt to existing matrix, located in the array by address (_N,_M)
	@param _N,_M: number of a cell in the array
	@param mt: adding matrix
	@return none
	*/
	void insert(unsigned int _N, unsigned int _M, const matrix & mt);

	/**
	@brief The function replaces existing matrix, located in the array by address (_N,_M), by the matrix mt
	@param _N,_M: number of a cell in the array
	@param mt: new matrix
	@return none
	*/
	void replace(unsigned int _N, unsigned int _M, const matrix & mt);

	/**
	@brief The function returns the maximal value of all elements contained in the array
	@param Arr: pointer to the array
	@return maximal element of the array
	*/
	friend double Max(const Arr2D & Arr);
	
	/**
	@brief The function sums up all matrixes contained in the array
	@param Arr: pointer to the array
	@return summed matrix
	*/
	friend matrix SumArr(const Arr2D &);
	
	/**
	@brief The function returns row counts
	@param Arr: pointer to the array
	@return row counts, (Arr.N)
	*/
	friend unsigned int StrArr(const Arr2D& Arr) 
		{ return Arr.N; }
	
	/**
	@brief The function returns column counts
	@param Arr: pointer to the array
	@return column counts, (Arr.M)
	*/
	friend unsigned int ColArr(const Arr2D& Arr) 
		{ return Arr.M; }
};
//------------------------------------------------------------------------------
/**
 @brief The array with (N-rows x M-columns) dimensions of
 small \b complex-value  matrixes with (n x m) dimensions.
 
 Example of using:
 @code
 ...
 Arr2DC Arr(0,0,0,0); //creating an array
 main() 
 {
   ...
   Arr = Arr2DC(1, 1, 2, 2); //changing the size
   Arr.ClearArr();   //clearing the array
   matrixC mt(2,2);  
   Arr.insert(0,0, mt); //adding matrix mt to the array
   matrixC M=Arr(0,0); //taking matrix M from the array
   ... 
 }
 @endcode
 */
class Arr2DC {
	unsigned int 
		N, M,	///< N,M - Dimensions of array
		n, m;	///< n,m - Dimensions of matrixes
	complex*** ptr;	///< Pointer to the array
	void AllocMem(complex*** = NULL);
	void FreeMem(void);
public:
	// constructors
	///Creates the array with (_N-rows x _M-columns) dimensions of small real-value matrixes with (_n x _m) dimensions
	Arr2DC(unsigned int _N, unsigned int _M, unsigned int _n, unsigned int _m)
		{ this->N=_N; this->M=_M; this->n=_n; this->m=_m; this->AllocMem(); }
	Arr2DC(void)
		{ this->N=0; this->M=0; this->n=0; this->m=0; this->AllocMem(); }
	Arr2DC(const Arr2DC &);
	// destructor
	~Arr2DC(void) { this->FreeMem(); }
	// members
	/// Returns the matrix stored in the array by address (_N,_M)
	matrixC operator() (unsigned int, unsigned int) const;
	
	///Returns a reference to the matrix element (_n,_m) stored in the array by address (_N,_M)
	complex& operator() (unsigned int _N, unsigned int _M, unsigned int _n, unsigned int _m) 
		{ return this->ptr[_N][_M][_n*this->n+_m]; }
	
	/// Returns the matrix element (_n,_m) stored in the array by address (_N,_M)
	complex operator() (unsigned int _N, unsigned int _M, unsigned int _n, unsigned int _m) const
		{ return this->ptr[_N][_M][_n*this->n+_m]; }
	Arr2DC operator=(const Arr2DC &);
	Arr2DC operator+=(const Arr2DC &);
	Arr2DC operator/=(double);
	
	/**
	@brief The function clears all elements of all matrixes in the array.
	*/
	void ClearArr(void) const;
	
	/**
	@brief The function adds the matrix mt to existing matrix, located in the array by address (_N,_M)
	@param _N,_M: number of a cell in the array
	@param mt: adding matrix
	@return none
	*/
	void insert(unsigned int, unsigned int, const matrixC &);
	
	/**
	@brief The function replaces existing matrix, located in the array by address (_N,_M), by the matrix mt
	@param _N,_M: number of a cell in the array
	@param mt: new matrix
	@return none
	*/
	void replace(unsigned int, unsigned int, const matrixC &);
	// friends
	
	/**
	@brief The function sums up all matrixes contained in the array
	@param Arr: pointer to the array
	@return summed matrix
	*/
	friend matrixC SumArr(const Arr2DC &);
	
	/**
	@brief The function returns row counts
	@param Arr: pointer to the array
	@return row counts, (Arr.N)
	*/
	friend unsigned int StrArr(const Arr2DC& Arr) { return Arr.N; }
	
	/**
	@brief The function returns column counts
	@param Arr: pointer to the array
	@return column counts, (Arr.M)
	*/
	friend unsigned int ColArr(const Arr2DC& Arr) { return Arr.M; }
};
//------------------------------------------------------------------------------
///@}
#endif
