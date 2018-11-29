/*
 * Executable program
 */
#ifndef EXECUTABLE_HH
#define EXECUTABLE_HH
#include <arg.hh>
class Executable
{
protected:
	int argc;
	char ** argv;
	char ** envp;

	arg::Parser parser; // the common parser
	std::string exec_name;

	//virtual void on_construct() {}
	virtual void set_parser() {}
	virtual void on_initialize() {}
	virtual int on_execute() {return 0;}
public:
	Executable(int argc, char ** argv, char ** envp) : argc(argc), argv(argv), envp(envp)
	{
		exec_name = argc ? argv[0] : "exec";
	}

	virtual ~Executable() {}

	int execute()
	{
		//on_construct();
		set_parser();
		parser.parse(argc, argv);
		on_initialize();
		return on_execute();
	}
};
#endif // EXECUTABLE_HH
