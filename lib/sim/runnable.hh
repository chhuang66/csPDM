/*
 * A runnable system
 */
#ifndef RUNNABLE_HH
#define RUNNABLE_HH
#include "with_param.hh"
#include <prm/name_set.hh>
#include <hdf5.h>
class Runnable :
	virtual public WithParam // running parameter manager
{
public:
	virtual ~Runnable() {}
	// essentials
	virtual void init() = 0; // initialize system
	virtual void dash() = 0; // short run
	// HDF5
	virtual void h5_save(hid_t group) const = 0;
	virtual void h5_load(hid_t group) = 0;
	virtual void h5_fr_buffer_load(hid_t group) = 0;
	// run meter
	virtual double count() const = 0;
	// optional abilities
	std::string description;
	virtual void initInputs() = 0; // initialize external inputs
};

/*
class Runnable_Original :
	virtual public WithParam // running parameter manager
{
public:
	virtual ~Runnable() {}
	// essentials
	virtual void init() = 0; // initialize system
	virtual void dash() = 0; // short run
	virtual void initInputs() = 0; // initialize external inputs 
	// HDF5
	virtual void h5_save(hid_t group) const = 0;
	virtual void h5_load(hid_t group) = 0;
	// run meter
	virtual double count() const = 0;
	// optional abilities
	NameSet dump_set;*/
	//virtual void dump(int /*what*/, std::ostream & /*output*/) const {} // dump states to a stream
	//std::string description;
//};

#endif // RUNNABLE_HH
