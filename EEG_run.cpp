#include "EEG_run.h"
#include <prm/hdf5.hh>
#include <vector>
#include <cassert>
#include <string.h>
#include <prm/ctrl_tag.hh>

void TransDelay::h5_save(hid_t group) const
{
    hid_t g = H5Gcreate(group,"TransDelay", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    prm::h5_save(g, param);
    H5Gclose(g);
}
void TransDelay::h5_load(hid_t group)
{
    hid_t g = H5Gopen(group,"TransDelay", H5P_DEFAULT);
    prm::h5_load(g, param);
    H5Gclose(g);

}

TransDelay::TransDelay()
{
	using namespace prm::tag;
	param.add_var("dThaCor", dThaCor = 8)
		<< Name("dThaCor")
		<< Desc("Transport delay between thalamus and cortex (msec)")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("dThaTha", dThaTha = 0)
		<< Name("dThaTha")
		<< Desc("Transport delay within thalamus (msec)")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("dCorCor", dCorCor = 0)
		<< Name("dCorCor")
		<< Desc("Transport delay within cortex (msec)")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
}

void ConnConst::h5_save(hid_t group) const
{
    hid_t g = H5Gcreate(group,"ConnConstants", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    prm::h5_save(g, param);
    H5Gclose(g);
}
void ConnConst::h5_load(hid_t group)
{
    hid_t g = H5Gopen(group, "ConnConstants", H5P_DEFAULT);
    prm::h5_load(g, param);
    H5Gclose(g);

}

ConnConst::ConnConst()
{
	using namespace prm::tag;
	param.add_var("extAMPA", extAMPA = 10)
		<< Name("extAMPA")
		<< Desc("Connection constant from the external via AMPA")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("extGABAa", extGABAa = 10)
		<< Name("extGABAa")
		<< Desc("Connection constant from the external via GABAa")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("extGABAb", extGABAb = 10)
		<< Name("extGABAb")
		<< Desc("Connection constant from the external via GABAb")
		<< Range(0, 1000)
		<< Save()
		<< CmdLine();

}

void EEGrun::dash()
{
    if (instant_jump_run+white_approx_run+mf_approx_run > 1){
	std::cerr << "More than one run type are set! Only one is required!\n";
	abort();
    }
    double exciIn, inhiIn, sinhiIn;
    double PY_toInter, PY_toIntra;
    
    for (size_t i = 0; i < nrnlist.size(); i++){
	 if (nrnlist[i] -> id == "PY"){
	     PY_toInter = nrnlist[i] -> frBuffer_inter[nrnlist[i] ->i_frBuffer_inter];
	     PY_toIntra = nrnlist[i] -> frBuffer_intra[nrnlist[i] ->i_frBuffer_intra];
	 }
     }
     

	if (ext_in){
	    //std::cerr << extIns[1] << std::endl;
	    exciIn = extIns[0];
	    inhiIn = extIns[1];
	    sinhiIn = extIns[2];
	    if (!(isEOF)){
		CsvParser_destroy_row(row);
		row = CsvParser_getRow(csvparser);
		if (!(row)){
		    if (!(strcmp(CsvParser_getErrorMessage(csvparser),"Reached EOF"))){
			isEOF = true;
		    }
		}
		else{
		    const char **rowFields = CsvParser_getFields(row);
		    for (size_t i = 0 ; i < nins ; i++) {
			extIns[i] = atof(rowFields[i]);
		    }
		}
	    }
	}
	else{
	    exciIn = inputs.get_input("PYexci", age, time_step);
	    inhiIn = inputs.get_input("PYinhi", age, time_step);
	    sinhiIn = inputs.get_input("PYslowinhi", age, time_step);
	}
     

    for (size_t i = 0; i < nrnlist.size(); i++){
	 if (nrnlist[i] -> get_id() == "PY"){
	     nrnlist[i]  -> dash(time_step, exciIn*connC.extAMPA, inhiIn*connC.extGABAa, sinhiIn*connC.extGABAb, instant_jump_run,white_approx_run, mf_approx_run, var_mat_run, exciIn, inhiIn, sinhiIn);
	 }
     }
     age += time_step;
     if (rec_every_step){
	  for (size_t i = 0; i < nrnlist.size(); i++)
	      nrnlist[i] -> store_state();
	  recT += recIntv;
        
    }
     else{
     if (age >= recT){
	  for (size_t i = 0; i < nrnlist.size(); i++)
	      nrnlist[i] -> store_state();
	  recT += recIntv;
     }
     }
}

void EEGrun::h5_save(hid_t group) const
{
    connC.h5_save(group);
    td.h5_save(group);
    inputs.h5_save(group);
    for (size_t i =0; i< nrnlist.size(); i++)
	nrnlist[i]->h5_save(group);
    prm::h5_save(group,param);
}

void EEGrun::h5_load(hid_t group)
{
    connC.h5_load(group);
    td.h5_load(group);
    inputs.h5_load(group);
    for (size_t i =0; i< nrnlist.size(); i++)
	nrnlist[i]->h5_load(group);
    prm::h5_load(group,param);
}

void EEGrun::h5_fr_buffer_load(hid_t group)
{
    for (size_t i =0; i< nrnlist.size(); i++)
	nrnlist[i]->h5_fr_den_load(group);

}

double EEGrun::count() const
{
    return get_age();
}

void EEGrun::init()
{
    for (size_t i = 0; i < nrnlist.size(); i++){
	nrnlist[i] -> initialize(time_step,recIntv);
	if (nrnlist[i] -> get_id() == "TC" || nrnlist[i] -> get_id() == "RE"){
	    nrnlist[i] -> make_buffer(time_step, td.dThaCor, td.dThaTha);
	    nrnlist[i] -> init_states();}
	else{
	    nrnlist[i] -> make_buffer(time_step, td.dThaCor, td.dCorCor);
	    nrnlist[i] -> init_states();}
    }
    inputs.initialize();
}

bool EEGrun::alter()
{
    bool al = false;
    for (size_t i = 0; i < nrnlist.size(); i++){
	al =al |  nrnlist[i] -> alter();
    }
    al = al | connC.ag.alter();
    al = al | td.ag.alter();
    al = al | inputs.alter();
    
    return al;
}


void EEGrun::addto_parser(arg::Parser& parser)
{
    
	for (size_t i = 0; i < nrnlist.size(); i++)
	    nrnlist[i] -> addto_parser(parser);
	connC.ag.clear();
	connC.ag.add(connC.param);
	parser.add_opt('c', "connection")
		.store(connC.ag.make_parser())
		.help("connection constants, `help' for a list", "PARAM=XX[,PARAM=XX...]");
	td.ag.clear();
	td.ag.add(td.param);
	parser.add_opt('d', "delay")
		.store(td.ag.make_parser())
		.help("transport delay, `help' for a list", "PARAM=XX[,PARAM=XX...]");
	
	parser.add_opt('i',"input")
		.store()
		.call(& extInputs::call_input_line, &inputs)
		.help("set external input, 'help' for usage", "Nrn=Module:PARAM=XX[,PARAM=XX...]");

	parser.add_opt("ext_inFile")
		.stow(ext_in_file)
		.set(& ext_in, true)
		.help("specify the file name for storing external inputs in csv format; delimiter=\\t ", "filename")
		.show_default();

}

void EEGrun::initInputs()
{
    if (ext_in){
	csvparser = CsvParser_new(ext_in_file.c_str(),  "\t", 0);
	if (!(row = CsvParser_getRow(csvparser))){
	    std::cerr << CsvParser_getErrorMessage(csvparser) << std::endl;
	    CsvParser_destroy(csvparser);
	    abort();
	}
	nins = CsvParser_getNumFields(row);
	extIns = new double[nins];
        const char **rowFields = CsvParser_getFields(row);
        for (size_t i = 0 ; i < nins ; i++) {
            extIns[i] = atof(rowFields[i]);
        }
    }
}

EEGrun::EEGrun()
{
	using namespace prm::tag;
	param.add_var("time_step", time_step = 0.1)
		<< Name("Time Step")
		<< Desc("Time step size for numerical integration (ms)")
		<< Range(0.0001, 10)
		<< Save()
		<< CmdLine();
 	param.add_var("recIntv", recIntv = 1.0)
		<< Name("Recording Iinterval")
		<< Desc("Time interval for storing states (ms)")
		<< Range(0.0001, 10)
		<< Save()
		<< CmdLine();
	param.add_var("age", age = 0.0)
		<< Name("age")
		<< Desc("starting time (ms)")
		<< Save()
		<< CmdLine();
	param.add_var("recT", recT = age+recIntv)
		<< Name("recT")
		<< Desc("next time instant for recording data")
		<< Save();
	param.add_var("instant_jump_run", instant_jump_run = false)
		<< Name("instant_jump_run")
		<< Desc("run instant conductance jump simulation (1st piority)")
		<< Control()
		<< Save()
		<< CmdLine();
	param.add_var("white_approx_run", white_approx_run = false)
		<< Name("white_approx_run")
		<< Desc("run white-approximation simulation")
		<< Control()
		<< Save()
		<< CmdLine();
	param.add_var("mf_approx_run", mf_approx_run = false)
		<< Name("mf_approx_run")
		<< Desc("run mean-field-approximation simulation")
		<< Control()
		<< Save()
		<< CmdLine();
	param.add_var("var_mat_run", var_mat_run = false)
		<< Name("var_mat_run")
		<< Desc("run variance matched")
		<< Control()
		<< Save()
		<< CmdLine();
	param.add_var("rec_every_step", rec_every_step = false)
		<< Name("rec_every_step")
		<< Desc("record states at every single time step")
		<< Control()
		<< Save()
		<< CmdLine();
	
	nrnlist.push_back(new pdmAEIF("PY"));

    if (rec_every_step){
        recT = time_step;
        recIntv = time_step;
    }
	ext_in_file = "";
	ext_in = false;
	isEOF = false;
	
	//recT = 0.0;
}

EEGrun::~EEGrun()
{
	for (std::vector<InitModule *>::iterator i = nrnlist.begin(); i != nrnlist.end(); i ++) {
		delete * i;
	}

    if (ext_in){
	CsvParser_destroy(csvparser);
	delete [] extIns;
    }

}

