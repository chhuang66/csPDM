#include "lif.h"

LIF::LIF(){
    	using namespace prm::tag;
	param.add_var("C", C = 2)
		<< Name("C")
		<< Desc("capacitance (uF)")
		<< Range(0, 10)
		<< Save()
		<< CmdLine();
	param.add_var("Gl", Gl = 0.1) 
		<< Name("Gl")
		<< Desc("conductance of leaky channel (mS)")
		<< Range(0, 100)
		<< Save()
		<< CmdLine();
	param.add_var("Vr", Vr = -65)
		<< Name("Vr")
		<< Desc("resetting voltage (mV)")
		<< Range(-90, 0)
		<< Save()
		<< CmdLine();
	param.add_var("Vl", Vl = -65)
		<< Name("Vl")
		<< Desc("reversal voltage for leaky input (mV)")
		<< Range(-90, 0)
		<< Save()
		<< CmdLine();
	param.add_var("Ve", Ve = 0)
		<< Name("Ve")
		<< Desc("reversal voltage for excitatory input (mV)")
		<< Range(-30, 30)
		<< Save()
		<< CmdLine();
	param.add_var("Vi", Vi = -70)
		<< Name("Vi")
		<< Desc("reversal voltage for inhibitory input (mV)")
		<< Range(-90, -30)
		<< Save()
		<< CmdLine();
	param.add_var("Vth", Vth = -55)
		<< Name("Vth")
		<< Desc("threshold voltage (mV)")
		<< Range(-90, 0)
		<< Save()
		<< CmdLine();
	param.add_var("t_ref", t_ref = 3)
		<< Name("t_ref")
		<< Desc("absolutely refractory period (msec)")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("mExci", mExci = 0.008)
		<< Name("mExci")
		<< Desc("average voltage jump for excitatory input")
		<< Range(0, 10)
		<< Save()
		<< CmdLine();
	param.add_var("mInhi", mInhi = 0.027)
		<< Name("mInhi")
		<< Desc("average voltage jump for inhibitrory input")
		<< Range(0, 10)
		<< Save()
		<< CmdLine();

	tau = C/Gl;
		
}