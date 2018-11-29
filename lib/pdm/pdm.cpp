#include <iostream>
#include<iomanip>
#include <cmath>
#include <vector>
#include "pdm.h"
#include <cstdlib>

#include <prm/hdf5.hh>

double InitModule::get_interFr()
{
    return frBuffer_inter[i_frBuffer_inter];
}

double InitModule::get_intraFr()
{
    return frBuffer_intra[i_frBuffer_intra];
}

void InitModule::addto_parser(arg::Parser& parser)
{
    ag.clear();
    ag.add(param);
    parser.add_opt(id+"param")
	    .store(ag.make_parser())
	    .help("set parameters of "+id+" population, `help' for usage", "PARAM=XX[,PARAM=XX,...]");  
}

std::string InitModule::get_id()
{
    return id;
}

void InitModule::h5_save_doubleArr(hid_t group, const double * arr, const int nx, char const * str)
{
	hsize_t dimf[] = {nx};
	hid_t stateS = H5Screate_simple(1, dimf, NULL);
	hid_t stateT = H5Tcopy(H5T_NATIVE_DOUBLE);
	hid_t stateD = H5Dcreate(group, str, stateT, stateS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status = H5Dwrite(stateD, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, arr);
	H5Tclose(stateT);	 
	H5Sclose(stateS);
	H5Dclose(stateD);
	return;
}

void InitModule::h5_load_doubleArr(hid_t group, double* arr, const int nx, const char* str)
{
        hid_t stateD = H5Dopen (group, str, H5P_DEFAULT);
	hid_t stateS = H5Dget_space (stateD);
	hsize_t dims[1];
 	int ndims = H5Sget_simple_extent_dims (stateS, dims, NULL);
	if (nx != dims[0]){
	    std::cerr << "size of buffer  has been changed. Abort()!!";
	    abort();
	}
	hid_t memtype = H5Tcopy(H5T_NATIVE_DOUBLE);
	double * data = new double [dims[0]];
	herr_t status = H5Dread (stateD, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,data);
	
	for (size_t i = 0 ; i < dims[0]; i++)
	    arr[i] = data[i];

	H5Tclose(memtype);	 
	H5Sclose(stateS);
	H5Dclose(stateD);
	
	delete [] data;
	
}

void InitModule::h5_save_intScalar(hid_t group, const size_t * arr, char const * str)
{
	//char const * s = str.c_str();
    	hsize_t dimf[] = {1};
	hid_t stateS = H5Screate_simple(1, dimf, NULL);
	hid_t stateT = H5Tcopy(H5T_NATIVE_INT);
	hid_t stateD = H5Dcreate2(group, str, stateT, stateS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status = H5Dwrite(stateD, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, arr);
	H5Tclose(stateT);	 
	H5Sclose(stateS);
	H5Dclose(stateD);
	return;

}

void InitModule::h5_load_intScalar(hid_t group, size_t* arr, const char* str)
{
        hid_t stateD = H5Dopen (group, str, H5P_DEFAULT);
	hid_t stateS = H5Dget_space (stateD);
	hsize_t dims[] ={1};
	hid_t memtype = H5Tcopy(H5T_NATIVE_INT);
	int * data = new int [1];
	herr_t status = H5Dread (stateD, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,data);
	
	arr[0] = data[0];

	H5Tclose(memtype);	 
	H5Sclose(stateS);
	H5Dclose(stateD);
	
	delete [] data;
	
}

void InitModule::h5_save_doubleMat(hid_t group, double** mat, const int nx, const int ny, const char* str)
{
	double data[nx][ny];
	for (size_t i = 0; i < nx; i++)
	    for (size_t j = 0; j < ny ; j++)
		data[i][j] = mat[i][j];
	hsize_t dimf[2] = {nx,ny};
	hid_t stateS = H5Screate_simple(2, dimf, NULL);
	hid_t stateT = H5Tcopy(H5T_NATIVE_DOUBLE);
	hid_t stateD = H5Dcreate2(group, str, stateT, stateS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status = H5Dwrite(stateD, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, data);
	H5Tclose(stateT);	 
	H5Sclose(stateS);
	H5Dclose(stateD);
	
	return;
}

void InitModule::h5_load_doubleMat(hid_t group, double** mat, const int nx, const int ny, const char* str)
{
        hid_t stateD = H5Dopen (group, str, H5P_DEFAULT);
	hid_t stateS = H5Dget_space (stateD);
	hsize_t dims[2];
 	int ndims = H5Sget_simple_extent_dims (stateS, dims, NULL);
	if (nx != dims[0]){
	    std::cerr << "size of buffer  has been changed. Abort()!!";
	    abort();
	}
	if (ny  != dims[1]){
	    std::cerr << "size of buffer  has been changed. Abort()!!";
	    abort();
	}
	hid_t memtype = H5Tcopy(H5T_NATIVE_DOUBLE);
	double data[nx][ny];
	herr_t status = H5Dread (stateD, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,data);
	
	for (size_t i = 0; i < dims[0]; i++)
	    for (size_t j = 0; j < dims[1] ; j++)
		mat[i][j] = data[i][j];

	H5Tclose(memtype);	 
	H5Sclose(stateS);
	H5Dclose(stateD);
	
	
}

InitModule::InitModule(std::string const & id) :
	id(id)
{
}

void pdmAEIF::initialize(double ts, double recIntv){
    if (allocMesh) {
	delete [] loc;
	delete [] den;
	delete []  kMat_l;
	delete []  kMat_exp;
	delete []  kMat_w;
	delete []  kMat_n;
	delete []  kMat_ad_ge;
	delete []  kMat_ad_gi;
	delete []  kMat_ad_gsi;
	delete []  kMat_diff_ge;
	delete []  kMat_diff_gi;
	delete []  kMat_diff_gsi;
	delete []  kMat_diff_n;

	delete []  jVec_l;
	delete []  jVec_exp;
	delete []  jVec_w;
	delete []  jVec_n;
	delete []  jVec_ad_ge;
	delete []  jVec_ad_gi;
	delete []  jVec_ad_gsi;
	delete []  jVec_diff_ge;
	delete []  jVec_diff_gi;
	delete []  jVec_diff_gsi;
	delete []  jVec_diff_n;
	allocMesh = false;
    }
    tcl = C/Gl;
    make_mesh();
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
    all_filters.push_back(new Filter(LPF,51,1./ts, 1./recIntv/5.));
}

void pdmAEIF::set_initDen()
{
    // used to setup initial values of density functions
    double dv = (Vc-Vlb)/N;
    den = new Vector2[N];
    double sumD = 0.0;
    for (size_t i = 0; i < N; i ++){
	//den[i][0] = 1.0 / sqrt(2*init_V_scale*init_V_scale*pi) * exp(-pow((loc[i][0]-init_V_loc),2)/2.0/init_V_scale/init_V_scale);
	//den[i][1] = 1.0 / sqrt(2*init_V_scale*init_V_scale*pi) * exp(-pow((loc[i][1]-init_V_loc),2)/2.0/init_V_scale/init_V_scale);
	//sumD += (den[i][0] + den[i][1])/2.0;
	den[i][0] = 0.0;
	den[i][1] = 0.0; 
    }
    den[iVl][0] = 1./dv;
    den[iVl][1] = 1./dv;
    /*
    for (size_t i = 0; i < N; i ++){
	den[i][0] = den[i][0]/sumD/dv;
	den[i][1] = den[i][1]/sumD/dv;
    }
    */
    
}

void pdmAEIF::set_loc()
{
    double dv = (Vc-Vlb)/N;
    loc = new Vector2[N];
    for (size_t i = 0; i < N; i ++){
	double v1 = i * dv + Vlb;
	double v2 = (i + 1) * dv + Vlb;
	if (Vr >= v1 && Vr < v2){
	    set_Vr(v1);
	    iVr = i;
	}
	if (Vl >= v1 && Vl < v2){
	    set_Vl(v1);
	    iVl = i;
	}
	loc[i][0] = v1;
	loc[i][1] = v2;
    }
}

void pdmAEIF::cal_MassStiffMatrix()
{
    double dv = (Vc-Vlb)/N;
    // mass matrix
    mMat = Matrix2(dv/3.0, dv/6.0, dv/6.0, dv/3.0);
    
    // Stiffness matrix
    kMat_l = new Matrix2[N];
    kMat_exp = new Matrix2[N];
    kMat_w = new Matrix2[N];
    kMat_n = new Matrix2[N];
    kMat_ad_ge = new Matrix2[N];
    kMat_ad_gi = new Matrix2[N];
    kMat_ad_gsi = new Matrix2[N];
    kMat_diff_ge = new Matrix2[N];
    kMat_diff_gi = new Matrix2[N];
    kMat_diff_gsi = new Matrix2[N];
    kMat_diff_n = new Matrix2[N];
    
    jVec_l = new Vector2[N];
    jVec_exp = new Vector2[N];
    jVec_w = new Vector2[N];
    jVec_n = new Vector2[N];
    jVec_ad_ge = new Vector2[N];
    jVec_ad_gi = new Vector2[N];
    jVec_ad_gsi = new Vector2[N];
    jVec_diff_ge = new Vector2[N];
    jVec_diff_gi = new Vector2[N];
    jVec_diff_gsi = new Vector2[N];
    jVec_diff_n = new Vector2[N];
    
    double e1,e2;
    double v1,v2;
    for (size_t i = 0; i <N; i++){
	v1 = loc[i][0];
	v2 = loc[i][1];
	// kMat_l
	e1 = -(Vl-v1)*dv/2.0+dv*dv/6.0;
	e2 = -(Vl-v2)*dv/2.0-dv*dv/6.0;
	kMat_l[i] = Matrix2(Gl/C/dv*e1, -Gl/C/dv*e1, Gl/C/dv*e2, -Gl/C/dv*e2);
	//std::cout << " kMat_l [" << i <<"]\n" <<  kMat_l[i] << std::endl;
	//kMat_exp
	e1 = -Dt*Dt*exp((v2-Vt)/Dt)+Dt*(dv+Dt)*exp((v1-Vt)/Dt);
	e2 = -Dt*Dt*exp((v1-Vt)/Dt)+Dt*(Dt-dv)*exp((v2-Vt)/Dt);
	kMat_exp[i] = Matrix2(Gl/C*Dt/dv/dv*e1, -Gl/C*Dt/dv/dv*e1,
			      Gl/C*Dt/dv/dv*e2, -Gl/C*Dt/dv/dv*e2);
	//std::cout << " kMat_exp [" << i <<"]\n" <<  kMat_exp[i] << std::endl;
	
	//kMat_w
	kMat_w[i] = Matrix2(1./C/2.0,-1./C/2.0,1./C/2.0,-1./C/2.0);
	//std::cout << " kMat_w [" << i <<"]\n" <<  kMat_w[i] << std::endl;

	//kMat_n
	kMat_n[i] = Matrix2(1./C/2.0,-1./C/2.0,1./C/2.0,-1./C/2.0);
	//std::cout << " kMat_n [" << i <<"]\n" <<  kMat_n[i] << std::endl;
	
	//kMat_ad_ge
	e1 = -(Ve-v1)*dv/2.0+dv*dv/6.0;
	e2 = -(Ve-v2)*dv/2.0-dv*dv/6.0;
	kMat_ad_ge[i] = Matrix2(1./C/dv*e1, -1./C/dv*e1, 1./C/dv*e2, -1./C/dv*e2);
	//std::cout << " kMat_ad_ge [" << i <<"]\n" <<  kMat_ad_ge[i] << std::endl;

	//kMat_ad_gi
	e1 = -(Vi-v1)*dv/2.0+dv*dv/6.0;
	e2 = -(Vi-v2)*dv/2.0-dv*dv/6.0;
	kMat_ad_gi[i] = Matrix2(1./C/dv*e1, -1./C/dv*e1, 1./C/dv*e2, -1./C/dv*e2);
	//std::cout << " kMat_ad_gi [" << i <<"]\n" <<  kMat_ad_gi[i] << std::endl;

	//kMat_ad_gsi
	e1 = -(Vsi-v1)*dv/2.0+dv*dv/6.0;
	e2 = -(Vsi-v2)*dv/2.0-dv*dv/6.0;
	kMat_ad_gsi[i] = Matrix2(1./C/dv*e1, -1./C/dv*e1, 1./C/dv*e2, -1./C/dv*e2);
	//std::cout << " kMat_ad_gsi [" << i <<"]\n" <<  kMat_ad_gsi[i] << std::endl;
	
	//kMat_diff_ge
	e1 = dv/3.*pow(Ve-v1,3)+1./12.*(pow(Ve-v2,4)-pow(Ve-v1,4));
	e2 = -dv/3.*pow(Ve-v2,3)-1./12.*(pow(Ve-v2,4)-pow(Ve-v1,4));
	kMat_diff_ge[i] = Matrix2(1./2./C/C/dv/dv*e1,-1./2./C/C/dv/dv*e1,
						1./2./C/C/dv/dv*e2,-1./2./C/C/dv/dv*e2);
	//std::cout << " kMat_diff_ge [" << i <<"]\n" <<  kMat_diff_ge[i] << std::endl;

	//kMat_diff_gi
	e1 = dv/3.*pow(Vi-v1,3)+1./12.*(pow(Vi-v2,4)-pow(Vi-v1,4));
	e2 = -dv/3.*pow(Vi-v2,3)-1./12.*(pow(Vi-v2,4)-pow(Vi-v1,4));
	kMat_diff_gi[i] = Matrix2(1./2./C/C/dv/dv*e1,-1./2./C/C/dv/dv*e1,
						1./2./C/C/dv/dv*e2,-1./2./C/C/dv/dv*e2);
	//std::cout << " kMat_diff_gi [" << i <<"]\n" <<  kMat_diff_gi[i] << std::endl;

	//kMat_diff_gsi
	e1 = dv/3.*pow(Vsi-v1,3)+1./12.*(pow(Vsi-v2,4)-pow(Vsi-v1,4));
	e2 = -dv/3.*pow(Vsi-v2,3)-1./12.*(pow(Vsi-v2,4)-pow(Vsi-v1,4));
	kMat_diff_gsi[i] = Matrix2(1./2./C/C/dv/dv*e1,-1./2./C/C/dv/dv*e1,
						1./2./C/C/dv/dv*e2,-1./2./C/C/dv/dv*e2);
	//std::cout << " kMat_diff_gsi [" << i <<"]\n" <<  kMat_diff_gsi[i] << std::endl;

	//kMat_diff_n
	kMat_diff_n[i] = Matrix2(1./2./C/C/2.,-1./2./C/C/2.,
						1./2./C/C/2.,-1./2./C/C/2.);
	//std::cout << " kMat_diff_n [" << i <<"]\n" <<  kMat_diff_n[i] << std::endl;
	
	//jVec_l
	e1 = Vl-v1;
	e2 = Vl-v2;
	jVec_l[i] = Vector2(Gl/C*e1,Gl/C*e2);
	//std::cout << " jVec_l [" << i <<"]\n" <<  jVec_l[i] << "\n" << std::endl;
	
	//jVec_exp
	e1 = exp((v1-Vt)/Dt);
	e2 = exp((v2-Vt)/Dt);
	jVec_exp[i] = Vector2(Gl/C*Dt*e1,Gl/C*Dt*e2);
	//std::cout << " jVec_exp [" << i <<"]\n" <<  jVec_exp[i] << "\n"<< std::endl;
	
	//jVec_ad_ge
	e1 = Ve-v1;
	e2 = Ve-v2;
	jVec_ad_ge[i] = Vector2(1./C*e1,1./C*e2);
	//std::cout << " jVec_ad_ge [" << i <<"]\n" <<  jVec_ad_ge[i] << "\n" << std::endl;
	
	//jVec_ad_gi
	e1 = Vi-v1;
	e2 = Vi-v2;
	jVec_ad_gi[i] = Vector2(1./C*e1,1./C*e2);
	//std::cout << " jVec_ad_gi [" << i <<"]\n" <<  jVec_ad_gi[i] << "\n" << std::endl;

	//jVec_ad_gsi
	e1 = Vsi-v1;
	e2 = Vsi-v2;
	jVec_ad_gsi[i] = Vector2(1./C*e1,1./C*e2);
	//std::cout << " jVec_ad_gsi [" << i <<"]\n" <<  jVec_ad_gsi[i] << "\n" << std::endl;
	
	//jVec_diff_ge
	e1 = -pow(Ve-v1,2);
	e2 = -pow(Ve-v2,2);
	jVec_diff_ge[i] = Vector2(1./2./C/C*e1,1./2./C/C*e2);
	//std::cout << " jVec_diff_ge [" << i <<"]\n" <<  jVec_diff_ge[i] << "\n" << std::endl;

	//jVec_diff_gi
	e1 = -pow(Vi-v1,2);
	e2 = -pow(Vi-v2,2);
	jVec_diff_gi[i] = Vector2(1./2./C/C*e1,1./2./C/C*e2);
	//std::cout << " jVec_diff_gi [" << i <<"]\n" <<  jVec_diff_gi[i] << "\n" << std::endl;
	
	//jVec_diff_gsi
	e1 = -pow(Vsi-v1,2);
	e2 = -pow(Vsi-v2,2);
	jVec_diff_gsi[i] = Vector2(1./2./C/C*e1,1./2./C/C*e2);
	//std::cout << " jVec_diff_gsi [" << i <<"]\n" <<  jVec_diff_gsi[i] << "\n" << std::endl;
	
	//jVec_diff_n
	e1 = -1.;
	e2 = -1.;
	jVec_diff_n[i] = Vector2(1./2./C/C*e1,1./2./C/C*e2);
	//std::cout << " jVec_diff_gsi [" << i <<"]\n" <<  jVec_diff_gsi[i] << "\n" << std::endl;
	
    }
    Ae = 1.-exp(-jGe/C);
    Ai = 1.-exp(-jGi/C);
    Asi = 1.-exp(-jGsi/C);
}


double pdmAEIF::step(double dt, double exS, double inS, double sinS, double refFlux, double refVar, double & frFlux_ref, double & frVar_ref, bool run_instant_jump, bool run_white_approx, bool run_mf_approx, bool run_var_mat)
{

    int temp_ku=0,temp_kl=0;
    double ** K = new double*[2*N];
    for(int i = 0; i < 2*N; i++){
        K[i] = new double[2*N];
        for(int j = 0; j < 2*N; j++){
            K[i][j] = 0.0;
        }
    }
    
    double h = (Vc-Vlb)/N;
    Matrix2 temp_K = Matrix2(1./h,-1./h,1./h,-1/h);
   
    double var_ge, var_gi, var_gsi;
    if (run_white_approx){
	var_ge = std_ge*std_ge;
	var_gi = std_gi*std_gi;
	var_gsi = std_gsi*std_gsi;
   }
    else{
	double g = Gl+mean_ge+mean_gi+mean_gsi;
	double tauv = C/Gl;
	double taum = C/g;
    double mv = get_meanV();
    double Dt_factor_ge = 1.-exp((mv-Vt)/Dt)+(mean_ge+mean_gi+mean_gsi)/Gl;
    double Dt_factor_gi = 1.-exp((mv-Vt)/Dt)+(mean_ge+mean_gi+mean_gsi)/Gl;
    double Dt_factor_gsi = 1.-exp((mv-Vt)/Dt)+(mean_ge+mean_gi+mean_gsi)/Gl;
	if (run_var_mat){
	    double scale_ge = tauv*taue/(tauv+taue*Dt_factor_ge)*(1.-a/(a+Gl)*taum/(taum+tauw));
	    double scale_gi = tauv*taui/(tauv+taui*Dt_factor_gi)*(1.-a/(a+Gl)*taum/(taum+tauw));
	    double scale_gsi = tauv*tausi/(tauv+tausi*Dt_factor_gsi)*(1.-a/(a+Gl)*taum/(taum+tauw));
	    var_ge = 2*std_ge*std_ge*scale_ge;
	    var_gi = 2*std_gi*std_gi*scale_gi;
	    var_gsi = 2*std_gsi*std_gsi*scale_gsi;}
	else{
	    double scale_ge = tauv*taue/(tauv+taue*Dt_factor_ge);
	    double scale_gi = tauv*taui/(tauv+taui*Dt_factor_gi);
	    double scale_gsi = tauv*tausi/(tauv+tausi*Dt_factor_gsi);
	    var_ge = 2*std_ge*std_ge*scale_ge;
	    var_gi = 2*std_gi*std_gi*scale_gi;
	    var_gsi = 2*std_gsi*std_gsi*scale_gsi;}
	
    }
    // calculating firing rate
    double temp_fr = 0.0;
    if (!run_instant_jump){
	    double flux = jVec_l[N-1][1]+jVec_exp[N-1][1]+jVec_ad_ge[N-1][1]*mean_ge+jVec_ad_gi[N-1][1]*mean_gi+jVec_ad_gsi[N-1][1]*mean_gsi-mean_w/C+mN/C;
	    if (flux > 0)
	        temp_fr += flux*den[N-1][1];
	    if (! run_mf_approx){
	        double s = jVec_diff_ge[N-1][1]*var_ge+jVec_diff_gi[N-1][1]*var_gi+jVec_diff_gsi[N-1][1]*var_gsi+jVec_diff_n[N-1][1]*stdN*stdN;
            temp_fr += -s/h*den[N-1][0] - s/h*den[N-1][1];
	    }
    }

    int * fr_porb = new int[N];
    if (run_instant_jump){
	for (int i = 0; i < N; i++){
	    if (i<N-1){
		Matrix2 temp;
		temp = kMat_l[i]+kMat_exp[i]+mean_w*kMat_w[i];
		K[2*i][2*i] += temp[0];
		K[2*i][2*i+1] += temp[2];
		K[2*i+1][2*i] += temp[1];
		K[2*i+1][2*i+1] += temp[3];


		double flux = jVec_l[i][1]+jVec_exp[i][1]-mean_w/C;
		if (flux>0){
		    K[2*i+1][2*i+1] -= flux;
		    K[2*i+2][2*i+1] += flux;
		}
		else{
		    K[2*i+1][2*i+2] -= flux;
		    K[2*i+2][2*i+2] += flux;	
		}
		
		double Vj = Vlb+i*h+h/2;
		K[2*i][2*i] -= exS*h/2.;
		K[2*i+1][2*i+1] -= exS*h/2.;
		double vjump = (Ve-Vj)*Ae;
		int numj = int(floor(vjump/h));
		if (numj>temp_ku) temp_ku = numj;
		if (i+numj >= N){
		    fr_porb[i] = 1;
		    if (! (t_ref>0)){
			K[2*iVr][2*i] += exS*h/2.;
			K[2*iVr+1][2*i+1] += exS*h/2.;
		    }
		}
		else{
		    fr_porb[i] = 0;
		    K[2*(i+numj)][2*i] += exS*h/2;
		    K[2*(i+numj)+1][2*i+1] += exS*h/2;
		}
		
		K[2*i][2*i] -= inS*h/2.;
		K[2*i+1][2*i+1] -= inS*h/2.;
		vjump = (Vi-Vj)*Ai;
		numj = int(floor(vjump/h));
		if (numj>temp_ku) temp_ku = numj;
		if (-numj > temp_kl) temp_kl = -numj;
		if (i+numj <0){
		    K[0][2*i] += inS*h/2.;
		    K[1][2*i+1] += inS*h/2.;
		}
		else{
		    K[2*(i+numj)][2*i] += inS*h/2.;
		    K[2*(i+numj)+1][2*i+1] += inS*h/2.;
		}

		K[2*i][2*i] -= sinS*h/2.;
		K[2*i+1][2*i+1] -= sinS*h/2.;
		vjump = (Vsi-Vj)*Asi;
		numj = int(floor(vjump/h));
		if (numj>temp_ku) temp_ku = numj;
		if (-numj > temp_kl) temp_kl = -numj;
		if (i+numj <0){
		    K[0][2*i] += sinS*h/2.;
		    K[1][2*i+1] += sinS*h/2.;
		}
		else{
		    K[2*(i+numj)][2*i] += sinS*h/2.;
		    K[2*(i+numj)+1][2*i+1] += sinS*h/2.;
		}		
		
	    }
	    else{
		Matrix2 temp;
		temp = kMat_l[i]+kMat_exp[i]+mean_w*kMat_w[i];
		K[2*i][2*i] += temp[0];
		K[2*i][2*i+1] += temp[2];
		K[2*i+1][2*i] += temp[1];
		K[2*i+1][2*i+1] += temp[3];


		double flux = jVec_l[i][1]+jVec_exp[i][1]-mean_w/C;
		if (flux>0){
		    K[2*i+1][2*i+1] -= flux;
		    if (! (t_ref>0))
			K[2*iVr][2*i+1] += flux;
		}
		
		double Vj = Vlb+i*h+h/2;
		K[2*i][2*i] -= exS*h/2.;
		K[2*i+1][2*i+1] -= exS*h/2.;
		double vjump = (Ve-Vj)*Ae;
		int numj = int(floor(vjump/h));
		if (numj>temp_ku) temp_ku = numj;
		if (i+numj >= N){
		    fr_porb[i] = 1;
		    if (! (t_ref>0)){
			K[2*iVr][2*i] += exS*h/2.;
			K[2*iVr+1][2*i+1] += exS*h/2.;
		    }
		}
		else{
		    fr_porb[i] = 0;
		    K[2*(i+numj)][2*i] += exS*h/2;
		    K[2*(i+numj)+1][2*i+1] += exS*h/2;
		}
		
		K[2*i][2*i] -= inS*h/2.;
		K[2*i+1][2*i+1] -= inS*h/2.;
		vjump = (Vi-Vj)*Ai;
		numj = int(floor(vjump/h));
		if (numj>temp_ku) temp_ku = numj;
		if (-numj > temp_kl) temp_kl = -numj;
		if (i+numj <0){
		    K[0][2*i] += inS*h/2.;
		    K[1][2*i+1] += inS*h/2.;
		}
		else{
		    K[2*(i+numj)][2*i] += inS*h/2.;
		    K[2*(i+numj)+1][2*i+1] += inS*h/2.;
		}		

		K[2*i][2*i] -= sinS*h/2.;
		K[2*i+1][2*i+1] -= sinS*h/2.;
		vjump = (Vsi-Vj)*Asi;
		numj = int(floor(vjump/h));
		if (numj>temp_ku) temp_ku = numj;
		if (-numj > temp_kl) temp_kl = -numj;
		if (i+numj <0){
		    K[0][2*i] += sinS*h/2.;
		    K[1][2*i+1] += sinS*h/2.;
		}
		else{
		    K[2*(i+numj)][2*i] += sinS*h/2.;
		    K[2*(i+numj)+1][2*i+1] += sinS*h/2.;
		}		
		
	    }
	}
	
    }
    else{
	for (int i = 0; i < N; i++){
	    if (i<N-1){
		Matrix2 temp;
		temp = kMat_l[i]+kMat_exp[i]+mean_w*kMat_w[i]+mean_ge*kMat_ad_ge[i]+mean_gi*kMat_ad_gi[i]+mean_gsi*kMat_ad_gsi[i]-mN*kMat_n[i];
		K[2*i][2*i] += temp[0];
		K[2*i][2*i+1] += temp[2];
		K[2*i+1][2*i] += temp[1];
		K[2*i+1][2*i+1] += temp[3];
		if (! run_mf_approx){
		    temp = (var_ge*kMat_diff_ge[i]+var_gi*kMat_diff_gi[i]+var_gsi*kMat_diff_gsi[i]+stdN*stdN*kMat_diff_n[i]) * temp_K;
		    K[2*i][2*i] += temp[0];
		    K[2*i][2*i+1] += temp[2];
		    K[2*i+1][2*i] += temp[1];
		    K[2*i+1][2*i+1] += temp[3];

		    K[2*i][2*i] -= var_ge*2/h*kMat_diff_ge[i][0]+var_gi*2/h*kMat_diff_gi[i][0]+var_gsi*2/h*kMat_diff_gsi[i][0]+stdN*stdN*2/h*kMat_diff_n[i][0];
		    K[2*i][2*i+2] +=  var_ge*2/h*kMat_diff_ge[i][2]+var_gi*2/h*kMat_diff_gi[i][2]+var_gsi*2/h*kMat_diff_gsi[i][2]+stdN*stdN*2/h*kMat_diff_n[i][2];
		    K[2*i+1][2*i] -=  var_ge*2/h*kMat_diff_ge[i][1]+var_gi*2/h*kMat_diff_gi[i][1]+var_gsi*2/h*kMat_diff_gsi[i][1]+stdN*stdN*2/h*kMat_diff_n[i][1];
		    K[2*i+1][2*i+2] +=  var_ge*2/h*kMat_diff_ge[i][3]+var_gi*2/h*kMat_diff_gi[i][3]+var_gsi*2/h*kMat_diff_gsi[i][3]+stdN*stdN*2/h*kMat_diff_n[i][3];	    
		    double s = jVec_diff_ge[i][1]*var_ge+jVec_diff_gi[i][1]*var_gi+jVec_diff_gsi[i][1]*var_gsi+jVec_diff_n[i][1]*stdN*stdN;
		    K[2*i+1][2*i] -= -s/h; 
		    K[2*i+1][2*i+1] -= -s/h; 
		    K[2*i+1][2*i+2] -= 2*s/h; 
		    K[2*i+2][2*i] += -s/h; 
		    K[2*i+2][2*i+1] += -s/h; 
		    K[2*i+2][2*i+2] += 2*s/h;
		    
		}
		double flux = jVec_l[i][1]+jVec_exp[i][1]+jVec_ad_ge[i][1]*mean_ge+jVec_ad_gi[i][1]*mean_gi+jVec_ad_gsi[i][1]*mean_gsi-mean_w/C+mN/C;
		if (flux>0){
		    K[2*i+1][2*i+1] -= flux;
		    K[2*i+2][2*i+1] += flux;
		}
		else{
		    K[2*i+1][2*i+2] -= flux;
		    K[2*i+2][2*i+2] += flux;	
		}
		
		
	    }
	    else{
		Matrix2 temp;
		temp = kMat_l[i]+kMat_exp[i]+mean_w*kMat_w[i]+mean_ge*kMat_ad_ge[i]+mean_gi*kMat_ad_gi[i]+mean_gsi*kMat_ad_gsi[i]-mN*kMat_n[i];
		K[2*i][2*i] += temp[0];
		K[2*i][2*i+1] += temp[2];
		K[2*i+1][2*i] += temp[1];
		K[2*i+1][2*i+1] += temp[3];

		if (! run_mf_approx){
		    temp = (var_ge*kMat_diff_ge[i]+var_gi*kMat_diff_gi[i]+var_gsi*kMat_diff_gsi[i]+stdN*stdN*kMat_diff_n[i]) * temp_K;
		    K[2*i][2*i] += temp[0];
		    K[2*i][2*i+1] += temp[2];
		    K[2*i+1][2*i] += temp[1];
		    K[2*i+1][2*i+1] += temp[3];

		    K[2*i][2*i] -= var_ge*2/h*kMat_diff_ge[i][0]+var_gi*2/h*kMat_diff_gi[i][0]+var_gsi*2/h*kMat_diff_gsi[i][0]+stdN*stdN*2/h*kMat_diff_n[i][0];
		    K[2*i+1][2*i] -=  var_ge*2/h*kMat_diff_ge[i][1]+var_gi*2/h*kMat_diff_gi[i][1]+var_gsi*2/h*kMat_diff_gsi[i][1]+stdN*stdN*2/h*kMat_diff_n[i][1];
		    double s = jVec_diff_ge[i][1]*var_ge+jVec_diff_gi[i][1]*var_gi+jVec_diff_gsi[i][1]*var_gsi+jVec_diff_n[i][1]*stdN*stdN;
		    K[2*i+1][2*i] -= -s/h; 
		    K[2*i+1][2*i+1] -= -s/h;
		    frVar_ref = -s/h;
		    if (! (t_ref>0)){
			K[2*iVr][2*i] += -s/h; 
			K[2*iVr][2*i+1] += -s/h;
		    }
		}
		double flux = jVec_l[i][1]+jVec_exp[i][1]+jVec_ad_ge[i][1]*mean_ge+jVec_ad_gi[i][1]*mean_gi+jVec_ad_gsi[i][1]*mean_gsi-mean_w/C+mN/C;
		if (flux>0){
		    K[2*i+1][2*i+1] -= flux;
		    if (! (t_ref>0))
			K[2*iVr][2*i+1] += flux;
		}		
	    }
	}
    }

    int flag = 1;
    int * p = new int[2*N];
    int ok;
    int size = 2*N;
    
    int ku;
    int kl;
    if (run_instant_jump){
	if (temp_ku*2>(N-iVr)*2)
	    ku = temp_ku*2+2;
	else
	    ku = (N-iVr)*2+2;
	kl = temp_kl*2+5;
    }
    else{
	ku =  (N-iVr)*2+2;
	kl = 5;
    }
    int lda = 2*kl+ku+1;
    int ldb = 2*N;

    double * A = new double[lda*2*N];
    for(int i = 0; i < lda*2*N; i++) A[i] = 0.0;
    double * B = new double[2*N];

    for(int j = 1; j < 2*N+1; j++){
        for (int i = int(fmax(1,j-ku)); i <= int(fmin(2*N,j+kl)); i++){
        if (i-1==j-1)
            A[kl+ku+i-j+(j-1)*lda]=1. - dt*K[i-1][j-1]*2./h;
        else
            A[kl+ku+i-j+(j-1)*lda]= - dt*K[i-1][j-1]*2./h;
        
        }
    }
    for(int i = 0; i < N; i++){
        B[2*i] = den[i][0];
        B[2*i+1] = den[i][1];
    }

    if (t_ref>0)
	B[2*iVr] += dt*refFlux*2/h;

    dgbsv_(&size, &kl, &ku, &flag, A, &lda, p, B, &ldb, &ok);

    for(int i = 0; i < N; i++){
        den[i][0] = B[2*i];
        den[i][1] = B[2*i+1];
    }
    
    slope_limiter(den);

    double sumD = 0.0;
    for(int i = 0; i < N; i++){
        sumD += den[i].ave()*h;
    }
    if (sumD>1.0){
    for(int i = 0; i < N; i++){
        den[i][0] /= sumD;
        den[i][1] /= sumD;
    }
    }
    
    if (run_instant_jump){
        temp_fr = 0.0;
	    double flux = jVec_l[N-1][1]+jVec_exp[N-1][1]-mean_w/C;
	    if (flux > 0)
	        temp_fr += flux*den[N-1][1];
	    for (int i = 0; i < N; i++){
	        if (fr_porb[i]>0)
		    temp_fr += exS*(den[i][0]+den[i][1])*h/2.;
	    }
    }
    frFlux_ref = temp_fr;
    for(int i = 0; i < 2*N; ++i){
        delete [] K[i];
    
    }
    delete [] K;
    delete [] A;
    delete [] B;
    delete [] p;
    delete [] fr_porb;
    
    return temp_fr;
}

void pdmAEIF::slope_limiter(Vector2 * density)
{
	double m1,m2,m3,d1,d2;
	for (size_t i = 0; i < N ; i++){
		if (i == 0){
			m1 = (density[i][0]+density[i][1])/2.0;
			m2 = (density[i][0]+density[i][1])/2.0;
			m3 = (density[i+1][0]+density[i+1][1])/2.0;
		}
		else if (i == N-1){
			m1 = (density[i-1][0]+density[i-1][1])/2.0;
			m2 = (density[i][0]+density[i][1])/2.0;
			m3 = 0.;
		}
		else{
			m1 = (density[i-1][0]+density[i-1][1])/2.0;
			m2 = (density[i][0]+density[i][1])/2.0;
			m3 = (density[i+1][0]+density[i+1][1])/2.0;
		}

		d1 = m2 - density[i][0];
		d2 = density[i][1] - m2;
			
		if ((d1 > 0) && (m3-m2 > 0) && (m2-m1 > 0))
			density[i][0] = m2 - fmin(fmin(d1,m3-m2), m2-m1);
		else if ((d1 < 0) && (m3-m2 < 0) && (m2-m1 < 0))
			density[i][0] = m2 + fmin(fmin(-d1,m2-m3),m1-m2);
		else
			density[i][0] = m2;

		if ((d2 > 0) && (m3-m2 > 0) && (m2-m1 > 0))
			density[i][1] = m2 + fmin(fmin(d2,m3-m2), m2-m1);
		else if ((d2 < 0) && (m3-m2 < 0) && (m2-m1 < 0))
			density[i][1] = m2 - fmin(fmin(-d2,m2-m3),m1-m2);
		else
			density[i][1] = m2;
		
		if (density[i][0] < 0 || density[i][1]<0){
		  if (m2>=0){
		    density[i][0] = m2;
		    density[i][1] = m2;
		  }
		  else{
		    if (density[i][0] < 0) density[i][0]=0.0;
		    if (density[i][1] < 0) density[i][1]=0.0;
		  }
		}
		
	}

}
void pdmAEIF::set_lnum(size_t n)
{
	loc = new Vector2[n];
	N = n;
}
void pdmAEIF::set_dnum(size_t n)
{
	den = new Vector2[n];
	N = n;
}

void pdmAEIF::make_buffer(double ts, double inter_delay, double intra_delay)
{
    size_t nBuffer;
    
    if (inter_delay == 0.0){
	frBuffer_inter = new double[1];
	frBuffer_inter[0] = 0.0;
	n_frBuffer_inter = 1;
    }
    else{
	nBuffer = int(inter_delay/ts);
	frBuffer_inter = new double[nBuffer];
	for (size_t i=0; i < nBuffer; i++)
	    frBuffer_inter[i] = 0.0;
	n_frBuffer_inter = nBuffer;
    }
    i_frBuffer_inter = 0;
    if (intra_delay == 0.0){
	frBuffer_intra = new double[1];
	frBuffer_intra[0] = 0.0;
	n_frBuffer_intra = 1;
    }
    else{
	nBuffer = int(intra_delay/ts);
	frBuffer_intra = new double[nBuffer];
	for (size_t i=0; i < nBuffer; i++)
	    frBuffer_intra[i] = 0.0;
	n_frBuffer_intra = nBuffer;
    }
    i_frBuffer_intra = 0;
    
    if (t_ref == 0.0){
	frFluxBuffer_ref = new double[1];
	frFluxBuffer_ref[0] = 0.0;
	frVarBuffer_ref = new double[1];
	frVarBuffer_ref[0] = 0.0;
	n_frBuffer_ref = 0;
    }
    else{
	nBuffer = int(t_ref/ts);
	frFluxBuffer_ref = new double[nBuffer];
	frVarBuffer_ref = new double[nBuffer];
	for (size_t i=0; i < nBuffer; i++){
	    frFluxBuffer_ref[i] = 0.0;
	    frVarBuffer_ref[i] = 0.0;
	}
	n_frBuffer_ref = nBuffer;
    }
    i_frBuffer_ref = 0;
}

void pdmAEIF::init_states()
{
    filtered_state.firingrate = fr;
    filtered_state.meanV = get_meanV();
    filtered_state.stdV = get_stdV();
    filtered_state.meanW = mean_w;
    filtered_state.meanGe = mean_ge;
    filtered_state.meanGi = mean_gi;
    filtered_state.meanGsi = mean_gsi;
    filtered_state.stdGe = std_ge;
    filtered_state.stdGi = std_gi;
    filtered_state.stdGsi = std_gsi;
    filtered_state.extExciInput = 0.0;
    filtered_state.extInhiInput = 0.0;
    filtered_state.extSlowInhiInput = 0.0;
    //store_state();

}

void pdmAEIF::w_dash(double dt)
{
    double mv = get_meanV();
    double nw = (mean_w+a*(mv-Vl)/tauw*dt+b*fr*dt)/(1+dt/tauw);
    mean_w = nw;
}

void pdmAEIF::ge_dash(double dt, double exS,bool run_white_approx)
{
    Matrix2 A, invA;
    Vector2 B, X;
    if (run_white_approx){
        mean_ge = jGe*exS;
        std_ge = sqrt(mean_ge*jGe);
    }
    else{
    A = Matrix2(1.+dt/taue,-2.*exS*jGe*dt,0.0,1.+2*dt/taue);
    B = Vector2(mean_ge,mean_ge*mean_ge+std_ge*std_ge)+dt*Vector2(exS*jGe,exS*jGe*jGe);
    invA = A.invert();
    X = invA*B;
        
    mean_ge = X[0];
    std_ge = sqrt(X[1]-X[0]*X[0]);
    }
}

void pdmAEIF::gi_dash(double dt, double inS,bool run_white_approx)
{
    Matrix2 A, invA;
    Vector2 B, X;
    if (run_white_approx){
        mean_gi = jGi*inS;
        std_gi = sqrt(mean_gi*jGi);
    }
    else{
    A = Matrix2(1.+dt/taui,-2.*inS*jGi*dt,0.0,1.+2*dt/taui);
    B = Vector2(mean_gi,mean_gi*mean_gi+std_gi*std_gi)+dt*Vector2(inS*jGi,inS*jGi*jGi);
    invA = A.invert();
    X = invA*B;
    
    mean_gi = X[0];
    std_gi = sqrt(X[1]-X[0]*X[0]);
    }

}

void pdmAEIF::gsi_dash(double dt, double sinS,bool run_white_approx)
{
    Matrix2 A, invA;
    Vector2 B, X;
    
    if (run_white_approx){
        mean_gsi = jGsi*sinS;
        std_gsi = sqrt(mean_gsi*jGsi);
    }
    else{
    A = Matrix2(1.+dt/tausi,-2.*sinS*jGsi*dt,0.0,1.+2*dt/tausi);
    B = Vector2(mean_gsi,mean_gsi*mean_gsi+std_gsi*std_gsi)+dt*Vector2(sinS*jGsi,sinS*jGsi*jGsi);
    invA = A.invert();
    X = invA*B;
    
    mean_gsi = X[0];
    std_gsi = sqrt(X[1]-X[0]*X[0]);
    }

}
void pdmAEIF::dash(double dt, double exS, double inS, double sinS, bool run_instant_jump, bool run_white_approx, bool run_mf_approx, bool run_var_mat, double extExciIn,double extInhiIn, double extSlowInhiIn)
{
    double frFlux = 0, frVar = 0;
    fr = step(dt, exS, inS, sinS, frFluxBuffer_ref[i_frBuffer_ref], frVarBuffer_ref[i_frBuffer_ref], frFlux, frVar, run_instant_jump, run_white_approx,run_mf_approx,run_var_mat);
    w_dash(dt);
    ge_dash(dt,exS,run_white_approx);
    gi_dash(dt,inS,run_white_approx);
    gsi_dash(dt,sinS,run_white_approx);
    //std::cout << fr << "\t" << get_meanV() << "\t"<<mean_w << "\t"<< mean_ge << "\t" << std_ge << "\t"<< mean_gi << "\t" << std_gi << std::endl;
    //std::cout << fr << std::endl;
    filtered_state.firingrate = fr;
    filtered_state.meanV = get_meanV();
    filtered_state.stdV = get_stdV();
    filtered_state.meanW = mean_w;
    filtered_state.meanGe = mean_ge;
    filtered_state.meanGi = mean_gi;
    filtered_state.meanGsi = mean_gsi;
    filtered_state.stdGe = std_ge;
    filtered_state.stdGi = std_gi;
    filtered_state.stdGsi = std_gsi;
    filtered_state.extExciInput = extExciIn;
    filtered_state.extInhiInput = extInhiIn;
    filtered_state.extSlowInhiInput = extSlowInhiIn;

    if (t_ref > 0.0){
	if (i_frBuffer_ref +1 == n_frBuffer_ref){
	    i_frBuffer_ref = 0;
	    frFluxBuffer_ref[n_frBuffer_ref-1] = frFlux;
	    frVarBuffer_ref[n_frBuffer_ref-1] = frVar;
	}
	else{
	    frFluxBuffer_ref[i_frBuffer_ref] = frFlux;
	    frVarBuffer_ref[i_frBuffer_ref] = frVar;
	    i_frBuffer_ref += 1;
	}
    }
    if (i_frBuffer_inter +1 == n_frBuffer_inter){
	i_frBuffer_inter = 0;
	frBuffer_inter[n_frBuffer_inter-1] = fr;
    }
    else{
	frBuffer_inter[i_frBuffer_inter] = fr;
	i_frBuffer_inter += 1;
    }
    if (i_frBuffer_intra +1 == n_frBuffer_intra){
	i_frBuffer_intra = 0;
	frBuffer_intra[n_frBuffer_intra-1] = fr;
    }
    else{
	frBuffer_intra[i_frBuffer_intra] = fr;
	i_frBuffer_intra += 1;
    }
    
}

void pdmAEIF::store_state()
{
    //state new_state;
    //new_state.firingrate = fr;
    //new_state.meanV = get_meanV();
    hist_states.push_back(filtered_state);
}

double pdmAEIF::get_meanV()
{
    double dv = (Vc-Vlb)/N;
    double sumD = 0, mv = 0;
    
    for (size_t i =0; i <N;i ++){
	sumD += den[i].ave() * dv;
	mv += den[i].ave() * dv * loc[i].ave();
    }
    /*
    if (sumD < 1.0)
        mv += (1-sumD) * Vr;
    else
        mv /= sumD;
    */
    mv /= sumD;
    return mv;
}

double pdmAEIF::get_stdV()
{
    double dv = (Vc-Vlb)/N;
    double sumD = 0, mv = 0,mv2 = 0;
    
    for (size_t i =0; i <N;i ++){
	sumD += den[i].ave() * dv;
	mv += den[i].ave() * dv * loc[i].ave();
    mv2 += den[i].ave() * dv * loc[i].ave()* loc[i].ave();
    }
    /*
    if (sumD < 1.0){
        mv += (1-sumD) * Vr;
        mv2 += (1-sumD) * Vr*Vr;
    }else{
        mv /= sumD;
        mv2 /= sumD;
    }
    */
    mv /= sumD;
    mv2 /= sumD;

    return sqrt(mv2-mv*mv);
}

void pdmAEIF::h5_save(hid_t group)
{
	char const * s = id.c_str();
	hid_t g = H5Gcreate(group, s, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	prm::h5_save(g, param);
	
	// storage of cells
	 prm::h5_save(g, cloc, "location");
	 prm::h5_save(g, cden, "density");
	 // extract data
	 state * hist_states_ptr = new state[hist_states.size()];
	for(size_t i = 0; i<hist_states.size();i++){
		hist_states_ptr[i].firingrate =hist_states.at(i).firingrate;
		hist_states_ptr[i].meanV =hist_states.at(i).meanV;
		hist_states_ptr[i].stdV =hist_states.at(i).stdV;
		hist_states_ptr[i].meanW =hist_states.at(i).meanW;
		hist_states_ptr[i].meanGe =hist_states.at(i).meanGe;
		hist_states_ptr[i].meanGi =hist_states.at(i).meanGi;
		hist_states_ptr[i].meanGsi =hist_states.at(i).meanGsi;
		hist_states_ptr[i].stdGe =hist_states.at(i).stdGe;
		hist_states_ptr[i].stdGi =hist_states.at(i).stdGi;
		hist_states_ptr[i].stdGsi =hist_states.at(i).stdGsi;
		hist_states_ptr[i].extExciInput =hist_states.at(i).extExciInput;
		hist_states_ptr[i].extInhiInput =hist_states.at(i).extInhiInput;
		hist_states_ptr[i].extSlowInhiInput =hist_states.at(i).extSlowInhiInput;
	}

	hsize_t dimf[] = {hist_states.size()};
	 
	hid_t stateS = H5Screate_simple(1, dimf, NULL);
	hid_t stateT = H5Tcreate (H5T_COMPOUND, sizeof(state));
	H5Tinsert(stateT, "firing_rate", HOFFSET(state, firingrate), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_V", HOFFSET(state, meanV), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_V", HOFFSET(state, stdV), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_w", HOFFSET(state, meanW), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_ge", HOFFSET(state, meanGe), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_gi", HOFFSET(state, meanGi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_gsi", HOFFSET(state, meanGsi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_ge", HOFFSET(state, stdGe), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_gi", HOFFSET(state, stdGi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_gsi", HOFFSET(state, stdGsi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "ext_exciIn", HOFFSET(state, extExciInput), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "ext_inhiIn", HOFFSET(state, extInhiInput), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "ext_slowinhiIn", HOFFSET(state, extSlowInhiInput), H5T_NATIVE_DOUBLE);
	hid_t stateD = H5Dcreate2(g, "hist_states", stateT, stateS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status = H5Dwrite(stateD, stateT, stateS, H5Dget_space(stateD),H5P_DEFAULT, hist_states_ptr);
	
	
	H5Tclose(stateT);	 
	H5Sclose(stateS);
	H5Dclose(stateD);

	 h5_save_doubleArr(g, frFluxBuffer_ref, n_frBuffer_ref, "frFluxBuffer_ref");
	 h5_save_doubleArr(g, frVarBuffer_ref, n_frBuffer_ref, "frVarBuffer_ref");
	 h5_save_doubleArr(g, frBuffer_intra, n_frBuffer_intra, "frBuffer_intra");
	 h5_save_doubleArr(g, frBuffer_inter, n_frBuffer_inter, "frBuffer_inter");
	 
	 h5_save_intScalar(g,&n_frBuffer_ref,"n_frBuffer_ref");
	 h5_save_intScalar(g,&i_frBuffer_ref,"i_frBuffer_ref");
	 h5_save_intScalar(g,&n_frBuffer_intra,"n_frBuffer_intra");
	 h5_save_intScalar(g,&i_frBuffer_intra,"i_frBuffer_intra");
	 h5_save_intScalar(g,&n_frBuffer_inter,"n_frBuffer_inter");
	 h5_save_intScalar(g,&i_frBuffer_inter,"i_frBuffer_inter");

	 // filter buffer
	 double ** filt_buffer = new double *[7];
	 size_t num_taps = all_filters.at(0)->get_num_taps();
	 for (size_t i = 0; i < 7; i++){
	     filt_buffer[i] = new double [all_filters.at(i)->get_num_taps()];
	     all_filters.at(i)->get_m(filt_buffer[i]);
	 }
	 h5_save_doubleMat(g, filt_buffer, 7, num_taps, "filter_buffers");

	H5Gclose(g);
	 delete [] hist_states_ptr;
	 for (size_t i = 0; i < 7; i++) delete [] filt_buffer[i];
	 delete [] filt_buffer;
	 

}
void pdmAEIF::h5_load(hid_t group)
{
	char const * s = id.c_str();
	hid_t g = H5Gopen(group, s, H5P_DEFAULT);
	prm::h5_load(g, param);
	
	 prm::h5_load(g, cloc, "location");
	 prm::h5_load(g, cden, "density");
	 // extract data
	 
	hid_t stateT = H5Tcreate (H5T_COMPOUND, sizeof(state));
	H5Tinsert(stateT, "firing_rate", HOFFSET(state, firingrate), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_V", HOFFSET(state, meanV), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_V", HOFFSET(state, stdV), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_w", HOFFSET(state, meanW), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_ge", HOFFSET(state, meanGe), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_gi", HOFFSET(state, meanGi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "mean_gsi", HOFFSET(state, meanGsi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_ge", HOFFSET(state, stdGe), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_gi", HOFFSET(state, stdGi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "std_gsi", HOFFSET(state, stdGsi), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "ext_exciIn", HOFFSET(state, extExciInput), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "ext_inhiIn", HOFFSET(state, extInhiInput), H5T_NATIVE_DOUBLE);
	H5Tinsert(stateT, "ext_slowinhiIn", HOFFSET(state, extSlowInhiInput), H5T_NATIVE_DOUBLE);

	 hid_t stateD = H5Dopen(g, "hist_states", H5P_DEFAULT);
	 hid_t stateS = H5Dget_space(stateD);
	hsize_t dims[1];
	H5Sget_simple_extent_dims(stateS, dims, NULL); 
	 
	 state * hist_states_ptr = new state[dims[0]];
	 herr_t status = H5Dread (stateD, stateT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hist_states_ptr);
	  
	for(size_t i = 0; i<dims[0];i++){
		hist_states.push_back(hist_states_ptr[i]);
	}	
	
	H5Tclose(stateT);	 
	H5Sclose(stateS);
	H5Dclose(stateD);

	H5Gclose(g);
	 delete [] hist_states_ptr;	
}

void pdmAEIF::h5_fr_den_load(hid_t group)
{
	char const * s = id.c_str();
	hid_t g = H5Gopen(group, s, H5P_DEFAULT);
	 prm::h5_load(g, cden, "density");

	h5_load_doubleArr(g, frFluxBuffer_ref, n_frBuffer_ref, "frFluxBuffer_ref");
	h5_load_doubleArr(g, frVarBuffer_ref, n_frBuffer_ref, "frVarBuffer_ref");
	 h5_load_doubleArr(g, frBuffer_intra, n_frBuffer_intra, "frBuffer_intra");
	 h5_load_doubleArr(g, frBuffer_inter, n_frBuffer_inter, "frBuffer_inter");
	 
	 h5_load_intScalar(g,&i_frBuffer_ref,"i_frBuffer_ref");
	 h5_load_intScalar(g,&i_frBuffer_intra,"i_frBuffer_intra");
	 h5_load_intScalar(g,&i_frBuffer_inter,"i_frBuffer_inter");
	 
	 // filter buffer
	 double ** filt_buffer = new double *[7];
	 size_t num_taps = all_filters.at(0)->get_num_taps();
	 for (size_t i = 0; i < 7; i++){
	     filt_buffer[i] = new double [all_filters.at(i)->get_num_taps()];
	 }
	 h5_load_doubleMat(g, filt_buffer, 7, num_taps, "filter_buffers");

	 for (size_t i = 0; i < 7; i++){
	     all_filters.at(i)->set_m(filt_buffer[i]);
	 }
	 for (size_t i = 0; i < 7; i++) delete [] filt_buffer[i];
	 delete [] filt_buffer;
	 
	 H5Gclose(g);
}

pdmAEIF::pdmAEIF(std::string const & id):
    InitModule(id)
{
	using namespace prm::tag;
	param.add_var("fr", fr = 0.0)
		<< Name("fr")
		<< Desc("mean firing rate (spikes/msec)")
		<< Save();
	param.add_var("meanW", mean_w = 0.0)
		<< Name("meanW")
		<< Desc("mean adaptation current (uA)")
		<< Save();
	param.add_var("meanGe", mean_ge = 0.0)
		<< Name("meanGe")
		<< Desc("mean excitatory synaptic conductance (mS)")
		<< Save();
	param.add_var("meanGi", mean_gi = 0.0)
		<< Name("meanGi")
		<< Desc("mean inhibitory synaptic conductance (mS)")
		<< Save();
	param.add_var("meanGsi", mean_gsi = 0.0)
		<< Name("meanGsi")
		<< Desc("mean slow inhibitory synaptic conductance (mS)")
		<< Save();
	param.add_var("stdGe", std_ge = 0.0)
		<< Name("stdGe")
		<< Desc("std excitatory synaptic conductance")
		<< Save();
	param.add_var("stdGi", std_gi = 0.0)
		<< Name("stdGi")
		<< Desc("std inhibitory synaptic conductance")
		<< Save();
	param.add_var("stdGsi", std_gsi = 0.0)
		<< Name("stdGsi")
		<< Desc("std slow inhibitory synaptic conductance")
		<< Save();
	param.add_var("init_V_loc", init_V_loc =-55)
		<< Name("init_V_loc")
		<< Desc("initial mean voltage (mV)")
		<< Save()
		<< CmdLine();
	param.add_var("init_V_scale", init_V_scale = 0.1)
		<< Name("init_V_scale")
		<< Desc("initial std voltage (mV)")
		<< Save()
		<< CmdLine();

	cloc.set_sizer(new prm::SzrFun<pdmAEIF>(* this, & pdmAEIF::set_lnum, & pdmAEIF::get_lnum));
	cden.set_sizer(new prm::SzrFun<pdmAEIF>(* this, & pdmAEIF::set_dnum, & pdmAEIF::get_dnum));

	prm::Struct * ss = new prm::Struct;
	ss->add_member("V1", & Vector2::x)
		<< Name("V1")
		<< Desc("left position");
	ss->add_member("V2", & Vector2::y)
		<< Name("V2")
		<< Desc("right position");
	ss->set_resolver(* this, & pdmAEIF::res_cloc);
	cloc.add_struct(ss);
	
	ss = new prm::Struct;
	ss->add_member("rho1", & Vector2::x)
		<< Name("rho1")
		<< Desc("left density");
	ss->add_member("rho2", & Vector2::y)
		<< Name("rho2")
		<< Desc("right density");
	ss->set_resolver(* this, & pdmAEIF::res_cden);
	cden.add_struct(ss);
	
    
}

pdmAEIF::~pdmAEIF()
{
    delete [] frFluxBuffer_ref;
    delete [] frVarBuffer_ref;
    delete [] frBuffer_intra;
    delete [] frBuffer_inter;
 }

