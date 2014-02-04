#ifndef MATRIX_HPP
#define MATRIX_HPP

#ifndef COMPL_HPP_
 #include "compl.hpp"
#endif

#ifndef _GLIBCXX_FSTREAM
 #include <fstream>
#endif

class matrix;
class matrixC;

/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */
/**
 @brief The array with (n-rows x m-columns) dimensions of \b real values.
 Size of the array can't be changed.

*/

class matrix {
	unsigned int n, m; ///< n, m - Dimensions of array
	double** ptr; ///< Pointer to the array
	void AllocMem(void); // memory operation
	void FreeMem(void);
	matrix(void) {} ///< This constructor MUST BE private!
public:
	/// Creates and initializes new matrix
	matrix(unsigned int _n, unsigned int _m)
		{ this->n=_n; this->m=_m; this->AllocMem(); }
	matrix(const matrix &); // copy constructor
	// destroys matrix
	virtual	~matrix() { this->FreeMem(); }
	// methods
	matrix	operator=(const matrix &);
	matrix	operator+(const matrix &) const;
	matrix	operator+=(const matrix &m) { return *this = m+*this; }
	matrix	operator-(const matrix &) const;
	matrix	operator-=(const matrix &m) { return *this = m-*this; }

	/// The operator returns a pointer to the i-th row of the array.
	double*	operator[] (unsigned int i) const { return this->ptr[i]; }
	matrix	operator*(double) const;
	matrix	operator*(const matrix &) const;
	matrix	operator*=(double);
	matrix	operator/(double) const;
	matrix	operator/=(double x) { return *this = *this/x; }
	bool	operator==(const matrix &) const;
	bool	operator!=(const matrix &m) const { return !(*this == m); }
	/// The function checks if the dimensions of the array are equal (n==m).
	bool	isSquare(void) const { return this->n == this->m; }
	/// The function sets all values of the array to x.
	void	Fill(double x);

	/// The function sets all the elements of the principal diagonal to one and all other elements to zero. The function do not check if the matrix is square.
	void	Identity(void);
	/// The function exchanges two rows
	void	Exchange(unsigned int i, unsigned int j) // exchanging of rows
		{ double* const p = this->ptr[i];
		this->ptr[i] = this->ptr[j]; this->ptr[j] = p; }

	// friends
	friend matrix operator*(double x, const matrix& m)
		{ return m*x; }
	/// The function returns a sum of squares of all elements
	friend double norm(const matrix &);
	/// The function returns the square of maximal value of the array
	friend double Max(const matrix &);
	/**
	@brief The function returns row counts
	@param m: the array
	@return row counts, (m.n)
	*/
	friend unsigned int Str(const matrix& m) { return m.n; }
	/**
	@brief The function returns column counts
	@param m: the array
	@return row counts, (m.m)
	*/
	friend unsigned int Col(const matrix& m) { return m.m; }

	/// The operator outputs all elements of the array to the file separated by space
	friend std::ofstream& operator<<(std::ofstream&, const matrix &);

};
//------------------------------------------------------------------------------

/**
@brief The array with (n-rows x m-columns) dimensions of \b complex values.
Size of the array can't be changed.

*/
class matrixC {
	unsigned int n, m; ///< n, m - Dimensions of array
	complex** ptr; ///< Pointer to the array
	void	AllocMem(void); // memory operation
	void	FreeMem(void);
	matrixC(void) {} ///< This constructor MUST BE private!
public:
	/// Creates and initializes new matrix
	matrixC(unsigned int _n, unsigned int _m)
		{ this->n=_n; this->m=_m; this->AllocMem(); }
	matrixC(const matrixC &); // copy constructor
	// destroys matixC
	virtual		~matrixC() { this->FreeMem(); }
	// methods
	matrixC		operator=(const matrixC &);
	matrixC		operator+(const matrixC &) const;
	matrixC		operator+=(const matrixC &m) { return *this = m+*this; }
	matrixC		operator-(const matrixC &) const;
	matrixC		operator-=(const matrixC &m) { return *this = m-*this; }

	/// The operator returns a pointer to the i-th row of the array.
	complex*	operator[](unsigned int N) const { return this->ptr[N]; }
	matrixC		operator*(const complex &) const;
	matrixC		operator*(const matrix &) const;
	matrixC		operator*(const matrixC &) const;
	matrixC		operator*=(const complex &z) { return *this = *this*z; }
	matrixC		operator/(const complex &) const;
	matrixC		operator/=(const complex &);
	bool		operator==(const matrixC &) const;
	bool		operator!=(const matrixC &m) const
	{ return !(*this == m); }
	/// The function sets all values of the array to z.
	void		Fill(const complex & z);
	/// The function sets all the elements of the principal diagonal to one and all other elements to zero. The function do not check if the matrix is square.
	void		Identity(void);
	/// The function exchanges two rows
	void		Exchange(unsigned int i, unsigned int j) // exchanging of rows
	{ complex* const p = this->ptr[i];
	this->ptr[i] = this->ptr[j]; this->ptr[j] = p; }
	// friends
	/// The function returns a sum of squares of all elements
	friend double norm(const matrixC &);
	/// The function returns the square of maximal value of the array
	friend double Max(const matrixC &);
	/**
	@brief The function returns row counts
	@param m: the array
	@return row counts, (m.n)
	*/
	friend unsigned int Str(const matrixC& m) { return m.n; }
	/**
	@brief The function returns column counts
	@param m: the array
	@return row counts, (m.m)
	*/
	friend unsigned int Col(const matrixC& m) { return m.m; }
	friend matrixC operator*(const complex& z, const matrixC &m)
	{ return m*z; }
	/// The operator outputs all elements of the array to the file separated by space. Real and Imagen parts are separated by comma.
	friend std::ofstream& operator<<(std::ofstream&, const matrixC &);
};
//------------------------------------------------------------------------------

#endif
/** @}
*/
