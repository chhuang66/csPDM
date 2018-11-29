#include "EEG_run.h"
#include <sim/runner.hh>
class RunnerEEGRun :
	private EEGrun,
	virtual public Runner
{
public:
	RunnerEEGRun(int argc, char ** argv, char ** envp) :
		Runner(argc, argv, envp)
	{}
};
int main(int argc, char * argv[], char * envp[]) {
	RunnerEEGRun obj(argc, argv, envp);
	return obj.execute();
}
