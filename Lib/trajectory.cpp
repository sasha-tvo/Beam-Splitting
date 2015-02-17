#include "trajectory.hpp"

std::string Chain::PathChain(void)
{
	std::string tr = "";
	for(std::list<unsigned int>::const_iterator it = this->Ch.begin(); it!=this->Ch.end(); it++)
	  tr += "_"+std::to_string(*it);
	return tr;
}

