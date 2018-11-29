/*
 *  Random number generator adapted from:
 *
 *    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/MTARCOK/mt19937ar-cok.c
 *
 *  See http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 */
#ifndef MT19937_HH
#define MT19937_HH
#include "random.hh"
#include <hdf5.h>
#include <string>
class MT19937 :
	virtual public Random
{
	unsigned long * state;
	int left;
	bool initf;
	unsigned long * next;
	void next_state();
public:
	MT19937();
	MT19937(const MT19937 & rng);
	~MT19937();
	void init(unsigned long seed = 0UL);
	double uniform();
#ifdef HAVE_HDF5
	void h5_save(hid_t group, const std::string & name = "mt19937") const;
	void h5_load(hid_t group, const std::string & name = "mt19937");
#endif
};
#endif
