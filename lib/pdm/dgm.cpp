#include <cmath>
#include "dgm.h"


void LGM::make_mesh()
{

    set_loc();
    cal_MassStiffMatrix();
    set_initDen();
     allocMesh = true;
    return;
}

LGM::LGM()
{
    	using namespace prm::tag;
	param.add_var("N", N = 600)
		<< Name("N")
		<< Desc("The number of meshes in V-direction")
		<< Range(0, 10000)
		<< Save()
		<< CmdLine();
	allocMesh = false;
}
LGM::~LGM()
{
    if (allocMesh){
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
	
    }
}
