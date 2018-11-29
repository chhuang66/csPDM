#ifndef HDF5_HH
#define HDF5_HH 1
#include "prm.hh"
extern "C" {
#include <hdf5.h>
}
namespace prm {
	class H5Error
	{
	};

	class TypeUnknownError :
		public H5Error
	{
	protected:
		std::string info;
	public:
		TypeUnknownError(const std::string & info);
		std::string get_info();
	};

	class DataMissingError :
		public H5Error
	{
	protected:
		std::string name;
	public:
		DataMissingError(const std::string & name);
		std::string get_name();
	};

	// for Param
	void h5_save(hid_t g, const Param & p);
	size_t h5_load(hid_t g, Param & p);
	// for Array
	void h5_save(hid_t g, const Array & a, const std::string & n);
	bool h5_load(hid_t g, Array & a, const std::string & n);
}
#endif // HDF5_HH
