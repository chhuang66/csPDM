#ifndef INPUT_MODULE_H
#define INPUT_MODULE_H
#include <vector>
#include <cmath>
#include <hdf5.h>
#include <sigc.hh>
#define pi M_PI
#include <prm/argu.hh>
#include <arg.hh>
#include <sim/mt19937.hh>
#include <sim/with_param.hh>
#include <cassert>
class baseModule{
    friend class InputModule;
protected:
	std::string id;
	MT19937 rng;
	size_t seed;
	prm::Param param;
	//prm::Argu ag;
public:
	baseModule(std::string const & id);
	std::string get_id();
	virtual ~baseModule() {}
	void initializeRng() {rng.init(seed);}
	virtual double get_value(double ct, double ts=0) = 0;
	//bool alter(){ag.alter();}
	//virtual void h5_save(hid_t group) = 0;
	void addto_parser(arg::Parser& parser);
	//virtual void h5_load(hid_t group) = 0;
};

class FixedConst:
    virtual public baseModule{
    protected:
	double amp;
    public:
	FixedConst();
	double get_value(double ct, double ts=0) override;
 };
 
class Sinusoidal:
    virtual public baseModule{
    protected:
	double amp;
	double freq;
	double phase;
	double mbg;
    public:
	Sinusoidal();
	double get_value(double ct, double ts=0) override; 
 };

class GaussianWhite:
    virtual public baseModule{
    protected:
	double m;
	double scale;
    public:
	GaussianWhite();
	double get_value(double ct, double ts=0) override; 
 };

class OUProcess:
    virtual public baseModule{
    protected:
        double m;
        double timeconst;
        double diffu;
        double cval;
    public:
	OUProcess();
	double get_value(double ct, double ts=0) override; 
	
	
 };
 
 class InputModule{
     friend class extInputs;
 protected:
     std::vector<baseModule *> mlist;
     size_t active;
     prm::Argu ag;
     std::string id;
 public:
     InputModule(std::string const & id);
     ~InputModule();
    //bool alter();
    //void initialize();
   bool set_module(const std::string & mod_line);
   void h5_save(hid_t group) const;
   void h5_load(hid_t group);
 };
 
 class extInputs{
 protected:
     std::vector<InputModule*> inputlist;
 public:
     extInputs();
     ~extInputs();
     double get_input(std::string name, double ct, double ts);
     bool alter();
     void initialize();
     static bool call_input_line(int, const std::string & mod_line, void * data);
     void h5_save(hid_t group) const;
     void h5_load(hid_t group);
};
#endif