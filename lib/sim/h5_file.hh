/* Take care of openning and backing up HDF5 files */
#include <hdf5.h>
#include <string>
class H5File
{
	std::string name;
	bool is_open;
	hid_t file;
	hid_t group;
public:
	struct Error
	{
		std::string info;
	};
	H5File(std::string name = "") : name(name), is_open(false) {}
	~H5File();
	void set_name(std::string n) {name = n;}
	hid_t open_save(std::string path = "/");
	void done_save();
	hid_t open_load(std::string path = "/");
	void close();
};
