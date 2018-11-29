#ifndef PDM_HH
#include <vector>
#include <cmath>
#include <pdm/dgm.h>
#include <pdm/nrn_models.h>
#include <hdf5.h>
#include <sigc.hh>
#define pi M_PI
#include <prm/argu.hh>
#include <arg.hh>
#include <sim/mt19937.hh>
#include <sim/with_param.hh>
#include "Vectors.h"
#include "Matrices.h"
#include "filt.h"
#include <lapacke.h>


/* DGESV prototype */
extern void dgbsv_( int *n, int *kl, int *ku, int *
	nrhs, double *ab, int *ldab, int *ipiv, double *b, 
	int *ldb, int *info );
class InitModule:
    virtual protected WithParam{
	friend class EEGrun;
protected:
	std::string id;
	prm::Argu ag;
	double * frBuffer_intra;	// firing rate buffer for transport delay in the intra-areas
	size_t n_frBuffer_intra, i_frBuffer_intra;
	double * frBuffer_inter;	// firing rate buffer for transport delay in the inter-areas
	size_t n_frBuffer_inter, i_frBuffer_inter;
	
public:
	InitModule(std::string const & id);
	std::string get_id();
	virtual ~InitModule() {}
	virtual void initialize(double ts, double recIntv) = 0;
	bool alter(){return ag.alter();}
	virtual void h5_save(hid_t group) = 0;
	void addto_parser(arg::Parser& parser);
	virtual void h5_load(hid_t group) = 0;
	virtual void h5_fr_den_load(hid_t group) = 0;
	virtual void make_buffer(double ts, double inter_delay, double intra_delay) = 0;
	virtual void init_states() = 0;
	double get_intraFr();
	double get_interFr();
	virtual void dash(double dt, double exS, double inS, double sinS, bool run_instant_jump, bool run_white_approx, bool run_mf_approx, bool run_var_mat, double extExciIn,double extInhiIn, double extSlowInhiIn) = 0;
	virtual void store_state() = 0;
	void h5_save_doubleArr(hid_t group, const double * arr, const int nx, char const * str);
	void h5_save_doubleMat(hid_t group, double ** mat, const int nx, const int ny, char const * str);
	void h5_save_intScalar(hid_t group, const size_t * arr, char const * str);
	void h5_load_doubleArr(hid_t group, double * arr, const int nx, char const * str);
	void h5_load_doubleMat(hid_t group, double ** mat, const int nx, const int ny, char const * str);
	void h5_load_intScalar(hid_t group, size_t * arr, char const * str);
};


class pdmAEIF:
    virtual public aEIF,
    virtual public LGM,
    virtual public InitModule
{
protected:
	double fr;		// firing rate
	double init_V_loc;		// setting Gaussian pdf for initial values, loc: mean; scale: std
	double init_V_scale;
	double mean_w;		// mean adaptation current
	double mean_ge;		// mean excitatory synaptic conductance
	double mean_gi;		// mean inhibitory synaptic conductance
	double mean_gsi;		// mean inhibitory synaptic conductance
	double std_ge, std_gi, std_gsi;    // standard deviations of all synaptic conductance variables
	double * frFluxBuffer_ref;	// firing rate buffer for refractoriness
	double * frVarBuffer_ref;	// firing rate buffer for refractoriness
	size_t n_frBuffer_ref, i_frBuffer_ref;

	
	double Ae, Ai, Asi;		// Total change in conductance 
	typedef struct {
		double firingrate;
		double meanV;
        double stdV;
		double meanW;
		double meanGe;
		double meanGi;
		double meanGsi;
		double stdGe;
		double stdGi;
		double stdGsi;
		double extExciInput;
		double extInhiInput;
		double extSlowInhiInput;
	}state;
	
	state filtered_state;
	std::vector<state> hist_states;
	std::vector<Filter *> all_filters;
public:
    pdmAEIF(std::string const & id);
    ~pdmAEIF();
    void initialize(double ts, double recIntv) override;
    void set_initDen() override;
    void set_loc() override;
    void cal_MassStiffMatrix() override;
    void h5_save(hid_t group) override;
    void h5_load(hid_t group) override;
    void h5_fr_den_load(hid_t group) override;
    double step(double dt, double exS, double inS, double sinS, double refFlux, double refVar, double & frFlux_ref, double & frVar_ref, bool run_instant_jump, bool run_white_approx,bool run_mf_approx, bool run_var_mat);
    void slope_limiter(Vector2 * density);
    size_t get_lnum() const {return N;}
    void set_lnum(size_t n);
    size_t get_dnum() const {return N;}
    void set_dnum(size_t n);
    void set_init_V_loc(double v) {init_V_loc = v;}
    void set_init_V_scale(double v) {init_V_scale = v;}
    double get_init_V_loc() const {return init_V_loc;}
    double get_init_V_scale() const {return init_V_scale;}
    void make_buffer(double ts, double inter_delay, double intra_delay) override;
    void init_states() override;
    void dash(double dt, double exS, double inS, double sinS, bool run_instant_jump, bool run_white_approx, bool run_mf_approx, bool run_var_mat, double extExciIn,double extInhiIn, double extSlowInhiIn) override;
    void store_state() override;
    double get_meanV();
    double get_stdV();
    void w_dash(double dt);
    void ge_dash(double dt, double exS,bool run_white_approx);
    void gi_dash(double dt, double  inS,bool run_white_approx);
    void gsi_dash(double dt, double sinS,bool run_white_approx);
protected:
	prm::Compound cloc;
	prm::Compound cden;
	prm::Strt * res_cloc(size_t n) {return new prm::Stt<Vector2>(loc[n]);}
	prm::Strt * res_cden(size_t n) {return new prm::Stt<Vector2>(den[n]);}
};


#endif
