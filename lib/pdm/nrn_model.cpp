#include "nrn_models.h"


void aEIF::set_C(double c){
	C = c;
	tcl = C/Gl;
}

void aEIF::set_Gl(double g){
	Gl = g;
	tcl = C/Gl;
}
aEIF::aEIF(){
    	using namespace prm::tag;
	param.add_var("C", C = 1)
		<< Name("C")
		<< Desc("capacitance (uF)")
		<< Range(0, 100)
		<< Save()
		<< CmdLine();
	param.add_var("Gl", Gl = 0.05) 
		<< Name("Gl")
		<< Desc("conductance of leaky channel (mS)")
		<< Range(0, 100)
		<< Save()
		<< CmdLine();
	param.add_var("Vr", Vr = -60)
		<< Name("Vr")
		<< Desc("resetting voltage (mV)")
		<< Range(-90, 0)
		<< Save()
		<< CmdLine();
	param.add_var("Vl", Vl = -60)
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
	param.add_var("Vi", Vi = -80)
		<< Name("Vi")
		<< Desc("reversal voltage for inhibitory input (mV)")
		<< Range(-100, -30)
		<< Save()
		<< CmdLine();
	param.add_var("Vsi", Vsi = -100)
		<< Name("Vsi")
		<< Desc("reversal voltage for slow inhibitory input (mV)")
		<< Range(-140, -30)
		<< Save()
		<< CmdLine();
	param.add_var("Vc", Vc = -40)
		<< Name("Vc")
		<< Desc("cutting voltage (mV)")
		<< Range(-90, 0)
		<< Save()
		<< CmdLine();
	param.add_var("Vt", Vt = -50)
		<< Name("Vt")
		<< Desc("threshold voltage (mV)")
		<< Range(-90, 0)
		<< Save()
		<< CmdLine();
	param.add_var("Dt", Dt = 2.5)
		<< Name("Dt")
		<< Desc("slope factor in exponential current (mV)")
		<< Range(0, 50)
		<< Save()
		<< CmdLine();
	param.add_var("Vlb", Vlb = -100)
		<< Name("Vlb")
		<< Desc("lower bound in V (mV)")
		<< Range(-200, -50)
		<< Save()
		<< CmdLine();
	param.add_var("t_ref", t_ref = 0)
		<< Name("t_ref")
		<< Desc("absolutely refractory period (msec)")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("jGe", jGe = 0.005)
		<< Name("jGe")
		<< Desc("average conductance jump for excitatory input")
		<< Range(0, 10)
		<< Save()
		<< CmdLine();
	param.add_var("jGi", jGi = 0.01)
		<< Name("jGi")
		<< Desc("average conductance jump for inhibitrory input")
		<< Range(0, 10)
		<< Save()
		<< CmdLine();
	param.add_var("jGsi", jGsi = 0.001)
		<< Name("jGsi")
		<< Desc("average conductance jump for slow inhibitrory input")
		<< Range(0, 10)
		<< Save()
		<< CmdLine();
	param.add_var("a", a = 0.005)
		<< Name("a")
		<< Desc("parameter a (subthreshold adaptation strength (mS))")
		<< Range(0, 50)
		<< Save()
		<< CmdLine();
	param.add_var("b", b = 0.05)
		<< Name("b")
		<< Desc("parameter b (spike adaptation strength (uA))")
		<< Range(0, 50)
		<< Save()
		<< CmdLine();
	param.add_var("taue", taue = 5)
		<< Name("taue")
		<< Desc("time constant for excitatory conductance (msec)")
		<< Range(0, 500)
		<< Save()
		<< CmdLine();
	param.add_var("taui", taui = 10)
		<< Name("taui")
		<< Desc("time constant for inhibitory conductance (msec)")
		<< Range(0, 500)
		<< Save()
		<< CmdLine();
	param.add_var("tausi", tausi = 100)
		<< Name("tausi")
		<< Desc("time constant for slow inhibitory conductance (msec)")
		<< Range(0, 500)
		<< Save()
		<< CmdLine();
	param.add_var("tauw", tauw = 600)
		<< Name("tauw")
		<< Desc("time constant for subthreshold adaptation (msec)")
		<< Range(0, 5000)
		<< Save()
		<< CmdLine();
	param.add_var("mN", mN = 0)
		<< Name("mN")
		<< Desc("mean value of current noise (uA)")
		<< Range(-100, 100)
		<< Save()
		<< CmdLine();
	param.add_var("stdN", stdN = 0)
		<< Name("stdN")
		<< Desc("std value of current noise (uA)")
		<< Range(0, 100)
		<< Save()
		<< CmdLine();
	tcl = C/Gl;    
}
