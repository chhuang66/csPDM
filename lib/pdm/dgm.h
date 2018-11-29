#ifndef DGM_H
#define DGM_H
#include <vector>
#include "Matrices.h"
#include "Vectors.h"
#include <sim/with_param.hh>
//必須加virtual destructor



class LGM:
	virtual protected WithParam{
	protected:
	    size_t N;				// the number of meshes on V

	    Vector2 * loc;			// nodes positions
	    Vector2 * den;			// nodal values
	    Matrix2 mMat;			// mass Matrix
	    Matrix2 * kMat_l;		// stiffness Matrix for leaky
	    Matrix2 * kMat_exp;		// stiffness Matrix for exponential
	    Matrix2 * kMat_w;		// stiffness Matrix for adaptive
	    Matrix2 * kMat_n;		// stiffness Matrix for mean current noise
	    Matrix2 * kMat_ad_ge;	// stiffness Matrix for excitatory advection
	    Matrix2 * kMat_ad_gi;	// stiffness Matrix for inhitory advection
	    Matrix2 * kMat_ad_gsi;	// stiffness Matrix for slow inhitory advection
	    Matrix2 * kMat_diff_ge;	// stiffness Matrix for excitatory diffusion
	    Matrix2 * kMat_diff_gi;	// stiffness Matrix for inhitory diffusion
	    Matrix2 * kMat_diff_gsi;	// stiffness Matrix for slow inhitory diffusion
	    Matrix2 * kMat_diff_n;	// stiffness Matrix for current noise

	    Vector2 * jVec_l;			// flux vector for leaky
	    Vector2 * jVec_exp;		// flux vector for exponential
	    Vector2 * jVec_w;		// flux vector for adaptive
	    Vector2 * jVec_n;		// flux vector for mean current noise
	    Vector2 * jVec_ad_ge;	// flux vector for excitatory advection
	    Vector2 * jVec_ad_gi;		// flux vector for inhitory advection
	    Vector2 * jVec_ad_gsi;	// flux vector for slow inhitory advection
	    Vector2 * jVec_diff_ge;	// flux vector for excitatory diffusion
	    Vector2 * jVec_diff_gi;	// flux vector for inhitory diffusion
	    Vector2 * jVec_diff_gsi;	// flux vector for slow inhitory diffusion
	    Vector2 * jVec_diff_n;		// flux vector for current noise

	    size_t iVr, iVl;
	    bool allocMesh;		// allocate memories for meshing or not
	public:
	    LGM();
	    virtual ~LGM();
	    void make_mesh();
	    
	    virtual void cal_MassStiffMatrix() = 0;
	    virtual void set_initDen() = 0;
	    virtual void set_loc() = 0;
	    size_t get_N() const {return N;}
	    void set_N(size_t n){N = n;}
	    
	    Vector2 * get_den() {return den;}
	    Vector2 * get_loc() {return loc;}
	    
};
#endif
