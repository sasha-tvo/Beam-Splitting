#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#ifndef _GLIBCXX_LIST
 #include<list>
#endif
/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */

/** @addtogroup AxFunc Axillary functions
 * @{
 */

//==============================================================================
	/**
	@brief The class represents chain of number of facets.

	This class allows to contain the trajectory that had been defined in input data file.
	This allows us to calculate only the trajectories we are interested in, thereby highly decrease calculation time.
	*/
class Chain {
public:
	std::list<unsigned int> Ch;
	unsigned int sz;
	Chain (const std::list<unsigned int>& _Ch) :
		   Ch(_Ch) { sz = this->Ch.size(); };
	virtual	~Chain() {};	
	unsigned int  Size(void) const { return this->sz; }
	std::list<unsigned int>::const_iterator  Begin(void) const
		{ return this->Ch.begin(); }
	std::list<unsigned int>::const_iterator  End(void) const
		{ return this->Ch.end(); }
};

/** @}
*/

/** @}
*/

#endif // TRAJECTORY_H
