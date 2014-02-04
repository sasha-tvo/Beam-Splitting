//==============================================================================
// base routines for treating Mueller and Johns matrices < 11.10.2002 >
//==============================================================================
#ifndef __Mueller_hpp_
#define __Mueller_hpp_

#ifndef MATRIX_HPP
#include "matrix.hpp"
#endif

/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */

/** @addtogroup AxFunc Axillary functions
 * @{
 */

	/**
	@brief The function returns Mueller matrix calculated from Jones matrix.

	The dimension of Mueller matrix is 4x4.
	The dimension of Jones matrix must be 2x2.
	@param in: complex-value matrix (matrixC)
	@return real-value matrix
	*/
matrix Mueller(const matrixC & in);

/**
 @brief Right multiplication of matrix \b m by <i> rotation matrix </i> with cos(f)=cs, sin(f)=sn.

 Rotation matrix
  @code
	| 1	0	0	0 |
	| 0	cs	sn	0 |
	| 0	-sn	cs	0 |
	| 0	0	0	1 |
 @endcode
*/
void RightRotateMueller(matrix&, double, double);

/**
 @brief Left multiplication of matrix \b m by <i> rotation matrix </i> with cos(f)=cs, sin(f)=sn.

 Rotation matrix
  @code
	| 1	0	0	0 |
	| 0	cs	sn	0 |
	| 0	-sn	cs	0 |
	| 0	0	0	1 |
 @endcode
*/
void LeftRotateMueller(matrix&, double, double);

/**
 @brief Multiplication of matrix \b m by left and right <i> rotation matrixes </i>

 Right rotation matrix
  @code
	| 1	0	0	0 |
	| 0	cs	sn	0 |
	| 0	-sn	cs	0 |
	| 0	0	0	1 |
 @endcode

  Left rotation matrix
  @code
	| 1	 0		0	0 |
	| 0	 _cs	_sn	0 |
	| 0	-(_sn)	_cs	0 |
	| 0	 0		0	1 |
 @endcode
*/
void RotateMueller(matrix& m, double _cs, double _sn, double cs, double sn);

/**
 @brief Calculate Mueller matrix in forward direction for randomly oriented particle. <b>See equation  below.</b>

 Mueller matrix in forward direction for randomly oriented particle
  @code
	| M11	0			0			0 	|
	| 0		(M22+M33)/2	0			0 	|
	| 0		0			(M22+M33)/2	0 	|
	| 0		0			0			M44 |
 @endcode

*/
void ForwardScattering(matrix &);

/**
 @brief Calculate Mueller matrix in backward direction for randomly oriented particle. <b>See equation  below. </b>

 Mueller matrix in backward direction for randomly oriented particle
  @code
	| M11	0			0			0 	|
	| 0		(M22-M33)/2	0			0 	|
	| 0		0			(M33-M22)/2	0 	|
	| 0		0			0			M44 |
 @endcode

*/
void BackwardScattering(matrix &);

#endif

/** @}
*/

/** @}
*/

