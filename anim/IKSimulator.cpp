#include "IKSimulator.h"

const char HEAD = 0;
const char BODY = 1;
const char ARM_L = 2;
const char ARM_U = 3;
const char LEG_L = 4; 
const char LEG_U = 5;
const char HAND = 6;
const char FOOT = 7;
const double ARM_U_LEN = 2.5;
const double ARM_L_LEN = 2.2;
const double SHOULDER_ROOT_LEN_X = 0.5;
const double SHOULDER_ROOT_LEN_Y = 2.5;
const double HAND_LEN = 0.75;

using namespace std;

IKSimulator::IKSimulator( const std::string& name, Hermite* target ):
	BaseSimulator( name ),
	m_hermite( target )
{
	for(int i = 0; i< 7; i++)
		m_guy.angle[i] = 0.0;
	
	m_n = 0;

	load_human();
	
	m_dist_from_board = 3.0;
	m_timestep = 0.0;

	for(int i = 0; i< 3; i++)
		for(int j = 0; j<7; j++)
			m_J[i][j]=0;
	for(int i = 0; i<4; i++)
		for(int j =0; j<4;j++)
			total[i][j] = 0.0;
	total[3][3] = 1.0;

	//initialize previous point as initial position
	m_prev_p.assign(SHOULDER_ROOT_LEN_X+ARM_U_LEN+ARM_L_LEN+HAND_LEN,SHOULDER_ROOT_LEN_Y-1,m_dist_from_board);

	//INITIALIZE HERMITE
	m_hermite = new Hermite( "hermite" );
	m_hermite->loadFromFile2D( "data/spline_hello.txt" );
	VectorObj temp_p;
	m_hermite->getPoint(temp_p,0.0);
	m_init = temp_p;


}	// IKSimulator

IKSimulator::~IKSimulator()
{
}	// IKSimulator::~SampleGravitySimulator

int IKSimulator::init(double time) 
	{ 

	return 0;
}


int IKSimulator::step(double time)
{
	
	
	
	//increment time step
	m_timestep = m_timestep+0.001;

	//when end is reached
	if(m_timestep >= 0.999){
		m_timestep = 0.0;
		m_prev_p = m_p;
		m_p = m_init;

		//animTcl::OutputMessage("INIT");
		return 0;
	}
	m_hermite->getPoint(m_p, m_timestep);
	m_p = m_p*0.5;
	VectorObj delta_p = m_p-m_prev_p ;
	
	//compute J
	compute_Jacobian();
	
	for(int k = 0; k < 7; k++){
		
		delta_theta[k] = 0.0;
		for(int j = 0; j < 3; j++)
		{
			
				delta_theta[k] += m_J_plus[k][j]*delta_p[j];
				//animTcl::OutputMessage("m_delta_x[%d] %lf", j, m_delta_x[j]);
		}
		
		m_guy.angle[k] += delta_theta[k];
		//animTcl::OutputMessage("angle%d %lf", k, m_guy.angle[k]);
	}

	VectorObj checker = m_prev_p;
	for(int k = 0; k < 3; k++){
		delta_p[k] = 0.0;
		for(int j = 0; j < 7; j++)
		{
				delta_p[k] += m_J[k][j]*delta_theta[j];
				//animTcl::OutputMessage("m_delta_x[%d] %lf", j, m_delta_x[j]);
		}
		
		checker[k] += delta_p[k];

	}
	//save state
	m_prev_p = m_p;
	return 0;

}// IKSimulator::step

int IKSimulator::command(int argc, char **argv){
	if(argc<1){
		animTcl::OutputMessage("system %s: wrong number of params.", m_name) ;
		return TCL_ERROR;
	}
	else if (strcmp(argv[0], "read") == 0){
		if(argc == 2){
			delete(m_hermite);
			m_n++;
			m_hermite = new Hermite("hermite"+ m_n);
			
			m_hermite->loadFromFile2D(argv[1]);
		//	bool success = GlobalResourceManager::use()->addSystem(m_hermite, true );
		//	assert( success );
			return TCL_OK;
		}
	}
	return TCL_ERROR;
}	//IKSimulator::command

void IKSimulator::getShoulderMatrix(double c[4][4], double t1, double t2, double t3){

	double mx [4][4] = {{1,0,0,0}, 
						{0,cos(t1), -sin(t1), 0}, 
						{0,sin(t1), cos(t1), 0}, 
						{0,0,0,1}};
	double my [4][4] = {{cos(t2),0,sin(t2),0}, 
						{0,1, 0, 0}, 
						{-sin(t2),0, cos(t2), 0}, 
						{0,0,0,1}};
	double mz [4][4] = {{cos(t3),-sin(t3),0,0}, {sin(t3),cos(t3), 0, 0}, {0,0, 1, 0}, {0,0,0,1}};
	double intermediate[4][4] = {0};
	intermediate[3][3] = 1.0;
	double intermediate2[4][4] = {0};
	intermediate2[3][3] = 1.0;
	MatrixMult4x4(intermediate, mx, my);
	MatrixMult4x4(intermediate2, intermediate,mz);

	//Translation
	double sh_tran[4][4]  = {{1,0,0,SHOULDER_ROOT_LEN_X},{0,1,0,SHOULDER_ROOT_LEN_Y},{0,0,1,0},{0,0,0,1}};
	MatrixMult4x4(c,sh_tran, intermediate2 );
}
void IKSimulator::getElbowMatrix(double c[4][4], double t4, double t5){

	double mx [4][4] = {{1,0,			0,	   0}, 
						{0,	cos(t4), -sin(t4), 0}, 
						{0, sin(t4), cos(t4),  0}, 
						{0,		0,		0,	   1}};
	double my [4][4] = {{cos(t5),0,sin(t5),0}, 
						{0,1, 0, 0}, 
						{-sin(t5),0, cos(t5), 0}, 
						{0,0,0,1}};
	double intermediate[4][4] ={0};
	intermediate[3][3] = 1.0;
	MatrixMult4x4(intermediate,mx, my);

	//Translation
	double el_tran[4][4]  = {{1,0,0,ARM_U_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	MatrixMult4x4(c,el_tran,intermediate);
}
void IKSimulator::getWristMatrix(double c[4][4], double t6, double t7){
	
	
	double my [4][4] = {{cos(t6),0,sin(t6),0}, 
						{0,1, 0, 0}, 
						{-sin(t6),0, cos(t6), 0}, 
						{0,0,0,1}};
	double mz [4][4] = {{cos(t7),-sin(t7),0,0}, 
						{sin(t7),cos(t7), 0, 0}, 
						{0, 0, 1, 0}, 
						{0, 0, 0, 1}};
	double intermediate[4][4] = {0};
	intermediate[3][3] = 1.0;
	MatrixMult4x4(intermediate,my, mz);
	//Translation
	double wr_tran[4][4]  = {{1,0,0,ARM_L_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	
	double d[4][4];
	MatrixMult4x4(d, wr_tran,intermediate );	
	//Apply Hand Translation (for convenience)
	double hand_tran[4][4]  = {{1,0,0,HAND_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	MatrixMult4x4(c, d, hand_tran);
}
inline void IKSimulator::MatrixMult4x4(double c[4][4], double a[4][4] , double b[4][4] ){

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			c[i][j] = 0.0;
			for(int k = 0; k < 4; k++)
			{
				c[i][j]+= a[i][k]*(b[k][j]);
			}
			
		}
	}
}
inline void IKSimulator::MatrixMult4x1(double c[4], double a[4][4] , double b[4] ){
	
	
	for(int i = 0; i< 4; i++){
		c[i] = 0.0;
		for(int j = 0; j < 4; j++){
			c[i] += a[i][j]*b[j];
		}
	}
	
}

void IKSimulator::getColumnJacobian(int i){
	
	double p_end_effector[4] = {0,0,0,1};
	if(i == 1){
		
		getWristMatrix(wr, m_guy.angle[5], m_guy.angle[6]);	//theta 6, theta 7
		getElbowMatrix(el, m_guy.angle[3], m_guy.angle[4]); //theta 4, theta 5
	
		MatrixMult4x4(elXwr, el, wr);
		double sh_x_theta1[4][4] = {{0,0,0,0},
									{0,-sin(m_guy.angle[0]), -cos(m_guy.angle[0]),0},
									{0, cos(m_guy.angle[0]), -sin(m_guy.angle[0]),0},
									{0,0,0,0}};
		double my [4][4] = {{cos(m_guy.angle[1]),0,sin(m_guy.angle[1]),0}, {0,1, 0, 0}, {-sin(m_guy.angle[1]),0, cos(m_guy.angle[1]), 0}, {0,0,0,1}};
		double mz [4][4] = {{cos(m_guy.angle[2]),-sin(m_guy.angle[2]),0,0}, {sin(m_guy.angle[2]),cos(m_guy.angle[2]), 0, 0}, {0,0, 1, 0}, {0,0,0,1}};
		double sh_tran[4][4]  = {{1,0,0,SHOULDER_ROOT_LEN_X},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		MatrixMult4x4(sh_yXsh_z, my , mz);
		MatrixMult4x4(sh_xXsh_yXsh_z,sh_x_theta1,sh_yXsh_z);
		MatrixMult4x4(sh,sh_tran, sh_xXsh_yXsh_z);//Translation
		MatrixMult4x4(total,sh, elXwr); //T_root = I
	
	}
	else if(i==2){
		getWristMatrix(wr, m_guy.angle[5], m_guy.angle[6]);	//theta 6, theta 7
		getElbowMatrix(el, m_guy.angle[3], m_guy.angle[4]); //theta 4, theta 5
		MatrixMult4x4(elXwr,el,wr);
		
		double mx [4][4] = {{1,0,0,0}, {0,cos(m_guy.angle[0]), -sin(m_guy.angle[0]), 0}, {0,sin(m_guy.angle[0]), cos(m_guy.angle[0]), 0}, {0,0,0,1}};
		double sh_y_theta2[4][4] = {{-sin(m_guy.angle[1]),0,cos(m_guy.angle[1]),0},
									{0, 0, 0, 0},
									{-cos(m_guy.angle[1]), 0, -sin(m_guy.angle[1]),0},
									{0, 0, 0, 0}};
		double mz [4][4] = {{cos(m_guy.angle[2]),-sin(m_guy.angle[2]),0,0}, {sin(m_guy.angle[2]),cos(m_guy.angle[2]), 0, 0}, {0,0, 1, 0}, {0,0,0,1}};
		double sh_tran[4][4]  = {{1,0,0,SHOULDER_ROOT_LEN_X},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		
		MatrixMult4x4(sh_yXsh_z, sh_y_theta2 , mz);
		MatrixMult4x4(sh_xXsh_yXsh_z,mx,sh_yXsh_z);
		MatrixMult4x4(sh,sh_tran, sh_xXsh_yXsh_z);//Translation
		MatrixMult4x4(total,sh, elXwr); //T_root = I
	}
	else if(i==3)
	{
		getWristMatrix(wr, m_guy.angle[5], m_guy.angle[6]);	//theta 6, theta 7
		getElbowMatrix(el, m_guy.angle[3], m_guy.angle[4]); //theta 4, theta 5
		MatrixMult4x4(elXwr, el, wr);
		
		double mx [4][4] = {{1,0,0,0}, {0,cos(m_guy.angle[0]), -sin(m_guy.angle[0]), 0}, {0,sin(m_guy.angle[0]), cos(m_guy.angle[0]), 0}, {0,0,0,1}};
		double my [4][4] = {{cos(m_guy.angle[1]),0,sin(m_guy.angle[1]),0}, 
							{0,1, 0, 0}, 
							{-sin(m_guy.angle[1]),0, cos(m_guy.angle[1]), 0}, 
							{0,0,0,1}};
		double sh_z_theta3[4][4] = {{-sin(m_guy.angle[2]),-cos(m_guy.angle[2]),0,0},
									{cos(m_guy.angle[2]),-sin(m_guy.angle[2]), 0,0},
									{0, 0, 0, 0},
									{0, 0, 0, 0}};

		double sh_tran[4][4]  = {{1,0,0,SHOULDER_ROOT_LEN_X},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	
		MatrixMult4x4(sh_yXsh_z, my , sh_z_theta3);
		MatrixMult4x4(sh_xXsh_yXsh_z,mx,sh_yXsh_z);
		MatrixMult4x4(sh,sh_tran, sh_xXsh_yXsh_z);//Translation
		MatrixMult4x4(total,sh, elXwr); //T_root = I
		//Translation
		

	}
	else if(i==4)
	{
		getWristMatrix(wr, m_guy.angle[5], m_guy.angle[6]);	//theta 6, theta 7
		getShoulderMatrix(sh,m_guy.angle[0], m_guy.angle[1], m_guy.angle[2]); //theta 1,2,3
		double el_x_theta4[4][4] = {{0,0,0,0},
									{0,-sin(m_guy.angle[3]), -cos(m_guy.angle[3]),0},
									{0, cos(m_guy.angle[3]), -sin(m_guy.angle[3]),0},
									{0,0,0,0}};
		double my[4][4] = {{cos(m_guy.angle[4]),0,sin(m_guy.angle[4]),0}, 
						{0,1, 0, 0}, 
						{-sin(m_guy.angle[4]),0, cos(m_guy.angle[4]), 0}, 
						{0,0,0,1}};
		double el_tran[4][4]  = {{1,0,0,ARM_U_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		
		MatrixMult4x4(el_yXwr, my, wr );
		MatrixMult4x4(el_xXel_yXwr, el_x_theta4, el_yXwr);
		MatrixMult4x4(elXwr, el_tran,el_xXel_yXwr);
		//Translation
	
		MatrixMult4x4(total, sh, elXwr);
		
		
	}
	else if(i==5)
	{
		getWristMatrix(wr,m_guy.angle[5], m_guy.angle[6]);	//theta 6, theta 7
		getShoulderMatrix(sh, m_guy.angle[0], m_guy.angle[1], m_guy.angle[2]); //theta 1,2,3
		double el_y_theta5[4][4] ={{-sin(m_guy.angle[4]),0,cos(m_guy.angle[4]),0},
									{0,0, 0,0},
									{-cos(m_guy.angle[4]), 0, -sin(m_guy.angle[4]),0},
									{0,0,0,0}};
		double mx[4][4] =  {{1,0,			0,	   0}, 
						{0,	cos(m_guy.angle[3]), -sin(m_guy.angle[3]), 0}, 
						{0, sin(m_guy.angle[3]), cos(m_guy.angle[3]),  0}, 
						{0,		0,		0,	   1}};
		double el_tran[4][4]  = {{1,0,0,ARM_U_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		MatrixMult4x4(el_yXwr, el_y_theta5, wr );
		MatrixMult4x4(el_xXel_yXwr, mx, el_yXwr);
		MatrixMult4x4(elXwr, el_tran,el_xXel_yXwr);
		//Translation
	
		MatrixMult4x4(total, sh, elXwr);
		
	
	}
	else if(i==6)
	{
		getShoulderMatrix(sh, m_guy.angle[0], m_guy.angle[1], m_guy.angle[2]); //theta 1,2,3
		getElbowMatrix(el, m_guy.angle[3], m_guy.angle[4]); //theta 4, theta 5
		
		double wr_y_theta6[4][4] = {{-sin(m_guy.angle[5]),0,cos(m_guy.angle[5]),0},
									{0,0, 0,0},
									{-cos(m_guy.angle[5]), 0, -sin(m_guy.angle[5]),0},
									{0,0,0,0}};
		double mz[4][4] = {{cos(m_guy.angle[6]),-sin(m_guy.angle[6]),0,0}, 
							{sin(m_guy.angle[6]),cos(m_guy.angle[6]), 0, 0}, 
							{0,0, 1, 0}, 
							{0,0,0,1}};
		double hand_tran[4][4] = {{1,0,0,HAND_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		double wr_tran[4][4] = {{1,0,0,ARM_L_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	
	
		MatrixMult4x4(wr_zXhand,mz, hand_tran);
		MatrixMult4x4(wr_yXwr_zXhand, wr_y_theta6, wr_zXhand);
		MatrixMult4x4(wr, wr_tran,wr_yXwr_zXhand);
		MatrixMult4x4(elXwr, el, wr);
		MatrixMult4x4(total, sh, elXwr);
	}
	else if(i==7)
	{
		getShoulderMatrix(sh, m_guy.angle[0], m_guy.angle[1], m_guy.angle[2]); //theta 1,2,3
		getElbowMatrix(el, m_guy.angle[3], m_guy.angle[4]); //theta 4, theta 5
		double wr_z_theta7[4][4] = {{-sin(m_guy.angle[6]),-cos(m_guy.angle[6]),0,0},
									{cos(m_guy.angle[6]),-sin(m_guy.angle[6]), 0,0},
									{0, 0, 0, 0},
									{0, 0, 0, 0}};
		double my[4][4] =  {{cos(m_guy.angle[5]),0,sin(m_guy.angle[5]),0}, 
							{0,1, 0, 0}, 
							{-sin(m_guy.angle[5]),0, cos(m_guy.angle[5]), 0}, 
							{0,0,0,1}};
		double hand_tran[4][4] = {{1,0,0,HAND_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		double wr_tran[4][4] = {{1,0,0,ARM_L_LEN},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	
		MatrixMult4x4(wr_zXhand,wr_z_theta7, hand_tran);
		MatrixMult4x4(wr_yXwr_zXhand, my, wr_zXhand);
		MatrixMult4x4(wr, wr_tran,wr_yXwr_zXhand);
		MatrixMult4x4(elXwr, el, wr);
		MatrixMult4x4(total, sh, elXwr);
	}

	MatrixMult4x1(delta_f_thetai, total, p_end_effector);
	
	for(int j = 0; j< 3; j++)
			m_J[j][i-1] = delta_f_thetai[j];
}
void IKSimulator::inverse3x3Matrix(double inv[3][3], double A[3][3]){
	double det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
				-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
				+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]); 
	
	if(det!=0){
		inv[0][0] = (A[1][1]*A[2][2]-A[1][2]*A[2][1])/det;
		inv[0][1] = (A[0][2]*A[2][1]-A[0][1]*A[2][2])/det;
		inv[0][2] = (A[0][1]*A[1][2]-A[0][2]*A[1][1])/det;
		inv[1][0] = (A[1][2]*A[2][0]-A[1][0]*A[2][2])/det;
		inv[1][1] = (A[0][0]*A[2][2]-A[0][2]*A[2][0])/det;
		inv[1][2] = (A[0][2]*A[1][0]-A[0][0]*A[1][2])/det;
		inv[2][0] = (A[1][0]*A[2][1]-A[1][1]*A[2][0])/det;
		inv[2][1] = (A[0][1]*A[2][0]-A[0][0]*A[2][1])/det;
		inv[2][2] = (A[0][0]*A[1][1]-A[0][1]*A[1][0])/det;
	}
	else{
		inv[0][0] = 1.0;
		inv[1][1] = 1.0;
		inv[2][2] = 1.0;
		inv[0][1] = 0.0;
		inv[0][2] = 0.0;
		inv[1][0] = 0.0;
		inv[1][2] = 0.0;
		inv[2][0] = 0.0;
		inv[2][1] = 0.0;

	}
	

}
void IKSimulator::compute_Jacobian(){
	//Get columns for Jacobian
	for(int i = 1; i<=7; i++)
		getColumnJacobian(i);
	//Get transpose of Jacobian
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 7; j++){
			m_J_trans[j][i] = m_J[i][j]; 
		//	animTcl::OutputMessage("[%d][%d] %lf",i,j, m_J[i][j]);
		}
	}
	//J * J_t
	double JJ_t[3][3];
	for(int k = 0; k < 3; k++){
			for(int i = 0; i< 3; i++)
			{
				JJ_t[k][i] = 0.0;
				for(int j = 0; j < 7; j++)
				{
					JJ_t[k][i] += m_J[k][j] *m_J_trans[j][i];	
				//	animTcl::OutputMessage("m_j[%d][%d] %lf", k,j,m_J[k][j]);
				}	
			}
	}
	//Invert JJ_t
	double JJ_t_inv[3][3];
	inverse3x3Matrix(JJ_t_inv, JJ_t);
	//Computer J+
	for(int k = 0; k < 7; k++){
			for(int i = 0; i< 3; i++)
			{
				m_J_plus[k][i] = 0.0;
				for(int j = 0; j < 3; j++)
				{
					m_J_plus[k][i] += m_J_trans[k][j]*JJ_t_inv[j][i];
					
				}
			}
	}

}
void IKSimulator::load_human(){
	int temp, a;
	temp = glGenLists(8);
	for (a=0; a<8; a++)
		m_guy.list[a] = temp + a;
	//Head
	glPushMatrix();
	glNewList(m_guy.list[HEAD], GL_COMPILE);
    glTranslatef(0,0.5,0);
    glScalef(1,1.5,1);
    glutSolidSphere(0.8,10,10);
	glEndList();
	glPopMatrix();
	//Body
	glPushMatrix();
	glNewList(m_guy.list[BODY], GL_COMPILE);
    glScalef(0.5,1.5,0.5);
    glutSolidSphere(2,10,10);
	glEndList();
	glPopMatrix();
	//Lower Arm
	glPushMatrix();
	glNewList(m_guy.list[ARM_L], GL_COMPILE);
	glTranslatef(ARM_L_LEN/2,0,0);
	glPushMatrix();
    glScalef(2,0.75,0.75);
    glutSolidSphere(.6,10,10);
	glPopMatrix();
	glEndList();
	glPopMatrix();
	//Upper Arm
	glPushMatrix();
	glNewList(m_guy.list[ARM_U], GL_COMPILE);
	glTranslatef(ARM_U_LEN/2,0,0);
	glPushMatrix();
    glScalef(2.0,0.75,.75);
    glutSolidSphere(.65,10,10);
	glPopMatrix();
	glEndList();
	glPopMatrix();
	//Lower Leg
	glPushMatrix();
	glNewList(m_guy.list[LEG_L], GL_COMPILE);
    glTranslatef(0,-3.75,0);
    glScalef(1,1.5,1);
    glutSolidSphere(.6,10,10);
	glEndList();
	glPopMatrix();
	//Upper Leg
	glPushMatrix();
	glNewList(m_guy.list[LEG_U], GL_COMPILE);
    glTranslatef(0,-4,0);
    glScalef(1,2,1);
    glutSolidSphere(.7,10,10);
	glEndList();
	glPopMatrix();
	//Hand
	glPushMatrix();
	glNewList(m_guy.list[HAND], GL_COMPILE);
    glTranslatef(HAND_LEN/2,0,0);
	glPushMatrix();
    glScalef(2,1,1);
    glutSolidSphere(.25,10,10);
	glPopMatrix();
	glEndList();
	glPopMatrix();
	//Foot
	glPushMatrix();
	glNewList(m_guy.list[FOOT], GL_COMPILE);
    glTranslatef(0,-.25,0);
    glScalef(1,.5,2);
    glutSolidSphere(.5,10,10);
	glEndList();
	glPopMatrix();
}

void IKSimulator::draw_head()
{
  glCallList(m_guy.list[HEAD]);
}

void IKSimulator::draw_hand()
{
  glPushMatrix();
    glCallList(m_guy.list[HAND]);
  glPopMatrix();
}  

void IKSimulator::draw_lowerarm()
{
	//glPushMatrix();
	
	glCallList(m_guy.list[ARM_L]);
	
	glTranslated(HAND_LEN,0,0);
	glRotated(m_guy.angle[5]*180/PI, 0.0,1.0,0.0);
	glRotated(m_guy.angle[6]*180/PI, 0.0, 0.0, 1.0);
	
	glCallList(m_guy.list[HAND]);
	
	
	
//	glPopMatrix();
	
}  

void IKSimulator::draw_arm()
{   
	glCallList(m_guy.list[ARM_U]);
	glTranslated(ARM_L_LEN/2,0, 0);
	glRotated(m_guy.angle[3]*180/PI, 1.0,0.0,0.0);
	glRotated(m_guy.angle[4]*180/PI,0.0,1.0,0.0);  
	
	
	draw_lowerarm();

}
void IKSimulator::draw_lowerleg()
{
	glPushMatrix();
		glCallList(m_guy.list[LEG_L]);
	glPopMatrix();
}
void IKSimulator::draw_leg()
{
	glPushMatrix();
		glCallList(m_guy.list[LEG_U]);
	glPopMatrix();

	glPushMatrix();
		glTranslated(0,-2.5,0);
		draw_lowerleg();
	glPopMatrix();
}
void IKSimulator::draw_human(){
	glPushMatrix();
    glCallList(m_guy.list[BODY]);
	glPopMatrix();

	glPushMatrix();
    glTranslatef(0,3.5,.25);
    draw_head();
	glPopMatrix();
	
	//MAIN IK

	
	//Right Arm
	
	glPushMatrix();
	glTranslatef(SHOULDER_ROOT_LEN_X,SHOULDER_ROOT_LEN_Y,0);
	
	
	glRotated(m_guy.angle[0]*180/PI, 1.0, 0.0, 0.0);
	
	glRotated(m_guy.angle[1]*180/PI, 0.0, 1.0, 0.0);
	glRotated(m_guy.angle[2]*180/PI, 0.0, 0.0, 1.0);
	draw_arm();
	
	glPopMatrix();


	//Left Arm
	glPushMatrix();
	glTranslated(-SHOULDER_ROOT_LEN_X,SHOULDER_ROOT_LEN_Y,0);
	glRotated(180,0.0,0.0,1.0);

	glPushMatrix();
			
			glPushMatrix();
				glCallList(m_guy.list[ARM_U]);
			glPopMatrix();
			glPushMatrix();
				glTranslatef(ARM_U_LEN,0, 0);
				glPushMatrix();
					glCallList(m_guy.list[ARM_L]);
				glPopMatrix();
				glPushMatrix();
				glTranslatef(ARM_L_LEN,0,0);	
					glPushMatrix();
					glCallList(m_guy.list[HAND]);
					glPopMatrix();
				glPopMatrix();
			glPopMatrix();
	glPopMatrix();
	glPopMatrix();

	//Legs
	glPushMatrix();
	glTranslatef(-0.6,0,0);
	draw_leg();
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.6, 0,0);
	draw_leg();
	glPopMatrix();
	
}
void IKSimulator::display(GLenum mode)
{

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glColor3f(1.0, 1.0, 1.0);
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix() ;
	//Draw Floor

	glColor3f(0.80392,0.3607843,0.3607843);
	glBegin(GL_POLYGON);
	glVertex3f(8,-8,10);
	glVertex3f(8,-8,-10);
	glVertex3f(-8,-8,-10);
	glVertex3f(-8,-8,10);
	glEnd();
	
	//Blackboard
	glColor3f(0.0,0.0,0.0);
	glBegin(GL_POLYGON);
	glVertex3f(5,-3,0);
	glVertex3f(5,3,0);
	glVertex3f(-5,3,0);
	glVertex3f(-5,-3,0);
	glEnd();

	//Wall
	glColor3f(0.545098,0.270588,0.0745098);
	glBegin(GL_POLYGON);
	glVertex3f(6,-8,0);
	glVertex3f(6,4,0);
	glVertex3f(-6,4,0);
	glVertex3f(-6,-8,0);
	glEnd();

	//Side Walls
	glBegin(GL_POLYGON);
	glVertex3f(6, -8,0);
	glVertex3f(6, -8, -5);
	glVertex3f(-6, -8, -5);
	glVertex3f(-6, -8, 0);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3f(6,-8,0);
	glVertex3f(6, 4, 0);
	glVertex3f(6, 4, -5);
	glVertex3f(6, -8, -5);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3f(-6,-8,0);
	glVertex3f(-6, 4, 0);
	glVertex3f(-6, 4, -5);
	glVertex3f(-6, -8, -5);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3f(6, -8,-5);
	glVertex3f(6, 4, -5);
	glVertex3f(-6, 4, -5);
	glVertex3f(-6, -8, -5);
	glEnd();

	glPopMatrix();
	
	glEnable(GL_LIGHTING);
	glPushMatrix();
	//Draw human
	glTranslatef(0,-1,m_dist_from_board);
	
	draw_human();
	glPopMatrix();
	

	//Draw Points
	glPushMatrix();
	glScalef(0.5,0.5,0.5);
	glColor3f(1.0,1.0,1.0);
	glDisable(GL_LIGHTING);
	glBegin(GL_POINTS);
	VectorObj temp_p;
	for(double t = 0.0; t< 1.0; t+= 0.001){
		temp_p = m_hermite->getIntermediatePoint(t);
		
		glVertex3f(temp_p[0], temp_p[1], temp_p[2]);
	}
	glEnd();
	glPopMatrix();
}
