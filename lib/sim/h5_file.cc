#include "h5_file.hh"
#include <iostream>
#include <cerrno>
#include <cstring>
extern "C" {
#include <fcntl.h>
#if defined (_WIN32)
#include <io.h>
#define F_OK (0)
#define R_OK (2)
#else
#include <unistd.h>
#endif
}
H5File::~H5File()
{
	close();
}

hid_t H5File::open_save(std::string path)
{
	close();
	std::string tmpfile = name + ".new~";
	file = H5Fcreate(tmpfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file < 0) {
		std::cerr << "fail to open hdf5 file: " << tmpfile << '\n';
		throw Error();
	}
	group = H5Gopen(file, path.c_str(), H5P_DEFAULT);
	is_open = true;
	return group;
}

void H5File::done_save()
{
	close();
	if (access(name.c_str(), F_OK) == 0) { // old data file exists?
		std::string bak = name + "~";
		int ret = rename(name.c_str(), bak.c_str()); // backup old file
		if (ret && access(bak.c_str(), F_OK) == 0) { // error renaming file, and old backup exists
			unlink(bak.c_str()); // try removing old backup first...
			rename(name.c_str(), bak.c_str()); // do again...
		}
	}
	rename((name + ".new~").c_str(), name.c_str());
}

hid_t H5File::open_load(std::string path)
{
	close();
	if (access(name.c_str(), R_OK) != 0) {
		int en = errno;
		std::cerr << "Can not read file `" << name
			  << "' due to error: \"" << strerror(en) << "\".\n";
		throw Error();
	}
	file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file >= 0) group = H5Gopen(file, path.c_str(), H5P_DEFAULT);
	else throw;
	is_open = true;
	return group;
}

void H5File::close()
{
	if (is_open) {
		H5Gclose(group);
		H5Fclose(file);
		is_open = false;
	}
}
