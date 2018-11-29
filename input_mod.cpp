    #include "input_mod.h"
#include <prm/hdf5.hh>
baseModule::baseModule(const std::string& id): id(id)
{
	using namespace prm::tag;
	param.add_var("seed", seed = 0)
		<< Name("seed")
		<< Desc("random seed")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();

}

double FixedConst::get_value(double ct, double ts)
{
    return amp;
}

FixedConst::FixedConst(): baseModule("Constant")
{
	using namespace prm::tag;
	param.add_var("amp", amp = 0)
		<< Name("amp")
		<< Desc("amplitude")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();

}

double Sinusoidal::get_value(double ct, double ts)
{
    return amp*(sin(2*pi*freq*ct/1000.+phase)+1.0)+mbg;
}


Sinusoidal::Sinusoidal(): baseModule("Sinusoidal")
{
	using namespace prm::tag;
	param.add_var("amp", amp = 0)
		<< Name("amp")
		<< Desc("amplitude")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	param.add_var("mbg", mbg = 0)
		<< Name("mbg")
		<< Desc("background input")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	param.add_var("freq", freq = 0)
		<< Name("freq")
		<< Desc("frequency")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	param.add_var("phase", phase = 0)
		<< Name("phase")
		<< Desc("phase")
		<< Range(0, 2*pi)
		<< CmdLine()
		<< Save();

}

double GaussianWhite::get_value(double ct, double ts)
{
    double x = m+rng.gaussian()*scale;
    if (x>=0.)
	return x;
    else
	return 0.0;
}

GaussianWhite::GaussianWhite(): baseModule("GaussianWhite")
{
	using namespace prm::tag;
	param.add_var("m", m = 0)
		<< Name("m")
		<< Desc("mean value")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	param.add_var("scale", scale = 0)
		<< Name("scale")
		<< Desc("standard deviation")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
}

double OUProcess::get_value(double ct, double ts)
{
    double x = cval * (1.0 - ts / timeconst) +  sqrt(diffu * ts) * rng.gaussian();
    cval = x;
    if (x+m>=0.)
	return x+m;
    else
	return 0.0;    
}

OUProcess::OUProcess(): baseModule("OUProcess")
{
	using namespace prm::tag;
	param.add_var("m", m = 0)
		<< Name("m")
		<< Desc("mean value")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	param.add_var("timeconst", timeconst = 1)
		<< Name("timeconst")
		<< Desc("time constant (msec)")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	param.add_var("diffu", diffu = 1)
		<< Name("diffu")
		<< Desc("diffusion constant")
		<< Range(0, 100000)
		<< CmdLine()
		<< Save();
	cval = 0.0;
}

bool InputModule::set_module(const std::string& mod_line)
{
    std::string::size_type k = mod_line.find(':');
    std::string mod_name = mod_line.substr(0, k);
    for (size_t i = 0; i < mlist.size(); i ++) if (mlist[i]->id == mod_name) {
	active = i;
	ag.clear();
	ag.add(mlist[i]->param);
	if (k != std::string::npos) {
		arg::SubParser * sp =ag.make_parser();
		sp->set(mod_line.substr(k+1));
		delete sp;
	}
	ag.alter();
	return true;
    }
    // show help
    bool req_help = mod_name == "help";
    if (! req_help) std::cout << "Unknown module: " << mod_name << '\n';
    std::cout << "\nList of available modules\n\n";
    for (size_t i = 0; i < mlist.size(); i ++) {
	std::cout << '\t' << mlist[i]->id << '\n';
    }
    std::cout << '\n';
    exit(req_help ? 0 : 1);
    
}

InputModule::InputModule(std::string const & id) :
	id(id), active(0)
{
	mlist.push_back(new FixedConst);
	mlist.push_back(new Sinusoidal);
	mlist.push_back(new GaussianWhite);
	mlist.push_back(new OUProcess);
}

InputModule::~InputModule()
{
    	for (std::vector<baseModule *>::iterator i = mlist.begin(); i != mlist.end(); i ++) {
		delete * i;
	}

}

void InputModule::h5_save(hid_t group) const
{
	char const * s = id.c_str();
	hid_t g = H5Gcreate(group, s, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t stype = H5Tcopy(H5T_C_S1);
	H5Tset_size(stype, H5T_VARIABLE);
	hid_t space = H5Screate(H5S_SCALAR);
	hid_t dataset = H5Dcreate(g, "mod_name", stype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	std::string mn = mlist[active]->id;
	char const * s1 = mn.c_str();
	H5Dwrite(dataset, stype,  H5S_ALL, H5S_ALL, H5P_DEFAULT, & s1);
	
	prm::h5_save(g, mlist[active]->param);
	H5Dclose(dataset);
	H5Sclose(space);
	H5Tclose(stype);
	H5Gclose(g);

}

void InputModule::h5_load(hid_t group)
{
	char const * s = id.c_str();
	hid_t g = H5Gopen(group, s, H5P_DEFAULT);
	hid_t dataset = H5Dopen(g, "mod_name", H5P_DEFAULT);
	hid_t filetype = H5Dget_type (dataset);
	size_t sdim = H5Tget_size (filetype);
	hid_t space = H5Dget_space(dataset);
	hsize_t dims[1] ={sdim};
	
	char * mn= new char[sdim];
	hid_t stype = H5Tcopy(H5T_C_S1);
	herr_t status = H5Tset_size (stype, H5T_VARIABLE);
	status = H5Dread (dataset, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, & mn);
	
	for (size_t i = 0; i < mlist.size(); i ++) if (mlist[i]->id == std::string(mn)) {
		active = i;
		ag.clear();
		ag.add(mlist[i]->param);
	}
	//H5Dvlen_reclaim(stype, space, H5P_DEFAULT, & mn);
	prm::h5_load(g, mlist[active]->param);

	H5Sclose(space);
	H5Tclose(stype);	
	H5Gclose(g);
	delete [] mn;

}


double extInputs::get_input(std::string name, double ct, double ts)
{
    for (size_t i = 0; i < inputlist.size(); i ++){
	if (inputlist[i] -> id == name){
	    double value = inputlist[i] -> mlist[inputlist[i]->active] -> get_value(ct, ts);
	    return value;
	}
    }
}

bool extInputs::alter()
{
    bool al = false;
    for (size_t i = 0; i < inputlist.size(); i ++){
	al = al || inputlist[i] -> ag.alter();
    }
    return al;
}

void extInputs::initialize()
{
    for (size_t i = 0; i < inputlist.size(); i ++){
	inputlist[i] -> mlist[inputlist[i]->active] -> initializeRng();
    }
    
}

extInputs::extInputs()
{
    inputlist.push_back(new InputModule("TCexci"));
    inputlist.push_back(new InputModule("TCinhi"));
    inputlist.push_back(new InputModule("REexci"));
    inputlist.push_back(new InputModule("REinhi"));
    inputlist.push_back(new InputModule("PYexci"));
    inputlist.push_back(new InputModule("PYinhi"));
    inputlist.push_back(new InputModule("PYslowinhi"));
}

extInputs::~extInputs()
{
    	for (std::vector<InputModule *>::iterator i = inputlist.begin(); i != inputlist.end(); i ++) {
		delete * i;
	}

}

bool extInputs::call_input_line(int, const std::string& mod_line, void * data)
{
    extInputs * in = static_cast<extInputs *>(data);
    std::string::size_type k = mod_line.find('=');
    std::string nrn_name = mod_line.substr(0, k);
    
	for (size_t i = 0; i < in -> inputlist.size(); i ++)
	    if (in ->inputlist[i] -> id == nrn_name)
		return in ->inputlist[i] -> set_module(mod_line.substr(k+1));
	// show help
	bool req_help = nrn_name == "help";
	if (! req_help) std::cout << "Unknown neural population: " << nrn_name << '\n';
	std::cout << "\nList of available populations\n\n";
	for (size_t i = 0; i < in ->inputlist.size(); i ++) {
		std::cout << '\t' << in ->inputlist[i]->id << '\n';
	}
	std::cout << '\n';
	exit(req_help ? 0 : 1);
    
}

void extInputs::h5_save(hid_t group) const
{
    hid_t g = H5Gcreate(group,"extInputs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (size_t i = 0; i < inputlist.size(); i++)
	inputlist[i] ->h5_save(g);
    H5Gclose(g);

}

void extInputs::h5_load(hid_t group)
{
    hid_t g = H5Gopen(group,"extInputs", H5P_DEFAULT);
    for (size_t i = 0; i < inputlist.size(); i++)
	inputlist[i] ->h5_load(g);
    H5Gclose(g);

}

