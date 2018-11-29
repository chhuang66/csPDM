#ifndef EEG_RUN_H
#define EEG_RUN_H
#include <pdm/pdm.h>
#include <sim/runnable.hh>
#include <sim/parser_add.hh>
#include <prm/argu.hh>
#include <sigc.hh>
#include "input_mod.h"
#include "csvparser.hh"

// define delivery delay of action potentials between different cortical areas
class TransDelay{
    friend class EEGrun;
protected:
    double dThaCor;
    double dThaTha;
    double dCorCor;
    prm::Argu ag;
    prm::Param param;
    
public:
    TransDelay();
    ~TransDelay(){}
    void h5_save(hid_t group) const;
    void h5_load(hid_t group);
    prm::Param & get_param() {return param;}
    prm::Argu & get_argu() {return ag;}
    
    
};
// define connection strength of all possible connections between different populations
class ConnConst{
    friend class EEGrun;
protected:
    double extAMPA;
    double extGABAa;
    double extGABAb;

    
    prm::Argu ag;
    prm::Param param;
public:
    ConnConst();
    ~ConnConst(){}
    void h5_save(hid_t group) const;
    void h5_load(hid_t group);
    prm::Param & get_param() {return param;}
    prm::Argu & get_argu() {return ag;}
};

class EEGrun:
	virtual public Runnable,
	virtual protected ParserAdd{
protected:
    double time_step;
    double recT;
    double recIntv;
    double age;
    std::vector<InitModule *> nrnlist;
    ConnConst connC;
    extInputs inputs;
    TransDelay td;
    bool instant_jump_run;
    bool white_approx_run;
    bool mf_approx_run;
    bool var_mat_run;
    bool rec_every_step;

    // external input
    std::string ext_in_file; 	// file storing external inputs
    int ext_in;			// indicate having external inputs or not
    CsvParser * csvparser;
    CsvRow * row;
    double * extIns;
    bool isEOF;
    size_t nins;
public:
	// overrides
	void init() override;
	void dash() override;
	
	void h5_save(hid_t group) const override;
	void h5_load(hid_t group) override;
	void h5_fr_buffer_load(hid_t group) override;
	double count() const override;
	void addto_parser(arg::Parser & parser) override;
	bool alter() override;
	void initInputs() override;
	EEGrun();
	~EEGrun();
	
    void go_next(){age += time_step;}
    double get_age() const {return age;}
    double get_ts() const {return time_step;}
    double get_recT() const {return recT;}
    double get_recIntv() const {return recIntv;}
    void set_ts(double t) {time_step = t;}
    void set_recT(double t) {recT = t;}
    void set_recIntv(double t) {recIntv = t;}   
};



#endif
