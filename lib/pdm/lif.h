#ifndef LIF_H
#define LIF_H
#include <sim/with_param.hh>
class LIF:
	virtual protected WithParam {
protected:
    double C;		// capacitance
    double Gl;		// contantance of leaky channel 
    double Vr;		// resetting voltage (mV)
    double Vl;		// reversal voltage for leaky input (mV)
    double Ve;		// reversal voltage for excitatory input (mV)
    double Vi;		// reversal voltage for inhibitory input (mV)
    double Vth;		// threshold voltage (mV)
    double tau;		// time constant (msec)
    double t_ref;	// absolutely refractory period (msec)
    double mExci;	// average voltage jump for excitatory input
    double mInhi;	// average voltage jump for inhibitrory input
    size_t Vr_node;
    size_t Vl_node;
public:
    LIF();
    virtual ~LIF(){}
    void set_Vr(const double v){Vr=v;}
    void set_Ve(const double v){Ve=v;}
    void set_Vl(const double v){Vl=v;}
    void set_Vi(const double v){Vi=v;}
    void setVth(const double v){Vth=v;}
    void set_tau(const double t){tau=t;}
    void set_t_ref(const double t){t_ref=t;}
    void set_mExci(const double a){mExci=a;}
    void set_mInhi(const double a){mInhi=a;}
    double get_Vr() const{return Vr;}
    double get_Vi() const{return Vi;}
    double get_Ve() const{return Ve;}
    double get_Vl() const{return Vl;}
    double get_Vth() const{return Vth;}
    double get_tau() const{return tau;}
    double get_t_ref() const{return t_ref;}
    double get_mExci() const{return mExci;}
    double get_mInhi() const{return mInhi;}
};
#endif
