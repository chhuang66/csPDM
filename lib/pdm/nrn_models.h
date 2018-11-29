#ifndef MODEL_H
#define MODEL_H
#include <sim/with_param.hh>
    
class aEIF:
    virtual protected WithParam {
protected:
    double C;		// capacitance
    double Gl;		// contantance of leaky channel 
    double Vr;		// resetting voltage (mV)
    double Vl;		// reversal voltage for leaky input (mV)
    double Ve;		// reversal voltage for excitatory input (mV)
    double Vi;		// reversal voltage for inhibitory input (mV)
    double Vsi;		// reversal voltage for slow inhibitory input (mV)
    double Vc;		// cutting voltage (mV)
    double Vt;		// threshold voltage (mV)
    double Dt;		// slope factor for exponential term
    double Vlb;		// lower bound in V
    double tcl;		// time constant (msec)
    double t_ref;	// absolutely refractory period (msec)
    double jGe;		// conductance jump for excitatory input
    double jGi;		// conductance jump for inhibitrory input
    double jGsi;		// conductance jump for slow inhibitrory input
    double a;		// parameter a for subthreshold adaptation
    double b;		// parameter b for spike adaptation
    double taue;	// time constant for excitatory
    double taui;		// time constant for inhibitory
    double tausi;	// time constant for slow inhibitory
    double tauw;	// time constant for subthreshold adaptation
    double mN;		// mean value of current noise
    double stdN;	// std value of current noise
    
    size_t Vr_node;
    size_t Vl_node;
public:
    aEIF();
    virtual ~aEIF(){}
    void set_C(double c);
    void set_Vr(double v){Vr=v;}
    void set_Ve(double v){Ve=v;}
    void set_Vl(double v){Vl=v;}
    void set_Vi(double v){Vi=v;}
    void set_Gl(double g);
    void set_Vc(double v){Vc=v;}
    void set_Vt(double v){Vt=v;}
    void set_Dt(double v){Dt=v;}
    void set_t_ref(double t){t_ref=t;}
    void set_jGe(double c){jGe=c;}
    void set_jGi(double c){jGi=c;}
    void set_jGsi(double c){jGsi=c;}
    void set_a(double c){a=c;}
    void set_b(double c){b=c;}
    void set_taue(double c){taue=c;}
    void set_taui(double c){taui=c;}
    void set_tausi(double c){tausi=c;}
    void set_tauw(double c){tauw=c;}
    
    double get_C() const{return C;}
    double get_Vr() const{return Vr;}
    double get_Vi() const{return Vi;}
    double get_Ve() const{return Ve;}
    double get_Vl() const{return Vl;}
    double get_Vc() const{return Vc;}
    double get_Vt() const{return Vt;}
    double get_Dt() const{return Dt;}
    double get_Gl() const{return Gl;}
    double get_tcl() const{return tcl;}
    double get_t_ref() const{return t_ref;}
    double get_jGe() const{return jGe;}
    double get_jGi() const{return jGi;}	
    double get_jGsi() const{return jGsi;}	
    double get_a() const{return a;}
    double get_b() const{return b;}	
    double get_taue() const{return taue;}
    double get_taui() const{return taui;}
    double get_tausi() const{return tausi;}	
    double get_tauw() const{return tauw;}	
};
#endif