#ifndef MY_IK_SIMULATOR_H
#define MY_IK_SIMULATOR_H

#include <GLModel/GLModel.h>
#include <shared/defs.h>
#include <util/util.h>
#include "animTcl.h"
#include "BaseSimulator.h"
#include "hermite.h"
#include "util/util.h"
#include <string>



struct Human { 
	double angle[7]; 
	int list[9];
};

// Inverse Kinematics Simulator

class IKSimulator : public BaseSimulator 
{
public:

	IKSimulator( const std::string& name, Hermite* target );
	~IKSimulator();

	int step(double time);
	void initialize_pos();


	void draw_human();
	void draw_head();
	void draw_hand();
	void draw_lowerarm();
	void draw_arm();
	void draw_leg();
	void draw_lowerleg();
	void load_human();
	void getColumnJacobian(int i);
	void compute_Jacobian();
	int init(double time);

	int command(int argc, char **argv) ;

	void getShoulderMatrix(double c[4][4],double t1, double t2, double t3);
	void getElbowMatrix(double c[4][4],double t4, double t5);
	void getWristMatrix(double c[4][4],double t6, double t7);
	void MatrixMult4x4(double c[4][4], double a[4][4], double b[4][4]);
	void inverse3x3Matrix(double inv[3][3], double A[3][3]);
	void MatrixMult4x1(double c[4] ,double a[4][4] , double b[4] );
	void display(GLenum mode);

protected:

	Vector m_pos0; // initial position
	Vector m_vel0; // initial velocity
	
	double m_J[3][7];
	double m_J_plus[7][3];
	double m_J_trans[7][3];
	double m_dist_from_board;

	double delta_theta[7];
	double ** m_temp_c;

	Hermite* m_hermite;

	Human m_guy;

	double m_timestep;
	VectorObj m_p;	//point position world frame
	VectorObj m_init;
	VectorObj m_prev_p;
	VectorObj m_delta_x;

	int m_n; //counter for hermite curves

	double sh[4][4];
	double wr[4][4];
	double el[4][4];
	double elXwr[4][4];
	double sh_xXsh_yXsh_z[4][4];
	double sh_yXsh_z[4][4];
	double shXsh_tran[4][4];
	double el_yXwr[4][4];
	double el_xXel_yXwr[4][4];
	
	double hand_tranXwr_tran[4][4];
	double wr_zXhand[4][4];
	double wr_yXwr_zXhand[4][4];
	double total[4][4];
	double delta_f_thetai[4];
};


#endif