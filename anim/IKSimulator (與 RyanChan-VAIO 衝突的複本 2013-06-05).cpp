#include "IKSimulator.h"

const char HEAD = 0;
const char BODY = 1;
const char ARM_L = 2;
const char ARM_U = 3;
const char LEG_L = 4; 
const char LEG_U = 5;
const char HAND = 6;
const char FOOT = 7;
const double ARM_U_LEN = 2.4;
const double ARM_L_LEN = 2.0;
const double SHOULDER_ROOT_LEN_X = 0.5;
const double SHOULDER_ROOT_LEN_Y = 2.5;
const double HAND_LEN = 0.75;

IKSimulator::IKSimulator( const std::string& name, Hermite* target ):
	BaseSimulator( name ),
	m_hermite( target )
{
	for(int i = 0; i< 7; i++)
		m_guy.angle[i] = 0.0;
	m_guy.angle[0] = 0;
	m_guy.angle[1] = 0;
	m_guy.angle[2] = -90;

	load_human();
	
	m_dist_from_board = 3.5;
	m_timestep = 0.0;

	for(int i = 0; i< 3; i++)
		for(int j = 0; j<7; j++)
			m_J[i][j]=0;
	for(int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
			m_J_plus[i][j] = 0;

	
	
}	// IKSimulator

IKSimulator::~IKSimulator()
{
}	// IKSimulator::~SampleGravitySimulator

int IKSimulator::step(double time)
{
	
	//Calculate point on curve (World-frame)
	if(m_timestep<1){
		m_hermite->getPoint(m_p,m_timestep);
		m_timestep+=0.001;
	}
	else m_timestep =0;
	m_p.assign(m_p.x()*0.5, m_p.y()*0.5, m_p.z()*0.5);
	return 0;

}// IKSimulator::step

int IKSimulator::command(int argc, char **argv){
	if(argc<1){
		animTcl::OutputMessage("system %s: wrong number of params.", m_name) ;
		return TCL_ERROR;
	}
	else if (strcmp(argv[0], "read") == 0){
		if(argc == 2){

			m_hermite->loadFromFile2D(argv[1]);
			bool success = GlobalResourceManager::use()->addSystem(m_hermite, true );
			assert( success );
			return TCL_OK;
		}
	}
	return TCL_ERROR;
}	//IKSimulator::command

double** IKSimulator::getShoulderMatrix(double t1, double t2, double t3){
	double** c = new double*[4];
	for(int i =0; i< 4; i++)
		c[i] = new double[4];
	double mx [4][4] = {{1,0,0,0}, {0,cos(t1), -sin(t1), 0}, {0,sin(t1), cos(t1), 0}, {0,0,0,1}};
	double my [4][4] = {{cos(t2),0,-sin(t2),0}, {0,1, 0, 0}, {sin(t2),0, cos(t2), 0}, {0,0,0,1}};
	double mz [4][4] = {{cos(t3),-sin(t3),0,0}, {sin(t3),cos(t3), 0, 0}, {0,0, 1, 0}, {0,0,0,1}};
	double intermediate[4][4];
	compRotMat4(intermediate, my, mx);
	
	//Translation
	double sh_tran[4][4]  = {{1,0,0,SHOULDER_ROOT_LEN_X},{0,1,0,SHOULDER_ROOT_LEN_Y},{0,0,1,0},{0,0,0,1}};
	c = MatrixMult4x4(intermediate,  sh_tran);
	
	return c;
}
double** IKSimulator::getElbowMatrix(double t4, double t5){
	double** c = new double*[4];
	for(int i =0; i< 4; i++)
		c[i] = new double[4];
	double mx [4][4] = {{1,0,			0,	   0}, 
						{0,	cos(t4), -sin(t4), 0}, 
						{0, sin(t4), cos(t4),  0}, 
						{0,		0,		0,	   1}};
	double my [4][4] = {{cos(t5),0,-sin(t5),0}, 
						{0,1, 0, 0}, 
						{sin(t5),0, cos(t5), 0}, 
						{0,0,0,1}};
	double intermediate[4][4];
	compRotMat4(intermediate,my, mx);

	//Translation
	double el_tran[4][4]  = {{1,0,0,0},{0,1,0,ARM_U_LEN},{0,0,1,0},{0,0,0,1}};
	c = MatrixMult4x4(intermediate,  el_tran);
	return c;
}
double** IKSimulator::getWristMatrix(double t6, double t7){
	double** c = new double*[4];
	for(int i =0; i< 4; i++)
		c[i] = new double[4];
	double my [4][4] = {{cos(t6),0,-sin(t6),0}, {0,1, 0, 0}, {sin(t6),0, cos(t6), 0}, {0,0,0,1}};
	double mz [4][4] = {{cos(t7),-sin(t7),0,0}, {sin(t7),cos(t7), 0, 0}, {0,0, 1, 0}, {0,0,0,1}};
	double intermediate[4][4];
	compRotMat4(intermediate,mz, my);

	//Translation
	double wr_tran[4][4]  = {{1,0,0,0},{0,1,0,ARM_L_LEN},{0,0,1,0},{0,0,0,1}};
	c = MatrixMult4x4(intermediate,  wr_tran);
	
	//Apply Hand Translation (for convenience)
	double hand_tran[4][4]  = {{1,0,0,0},{0,1,0,HAND_LEN},{0,0,1,0},{0,0,0,1}};
	c = MatrixMult4x4(hand_tran,(double(*)[4])c);
	return c;
}
double ** IKSimulator::MatrixMult4x4(double a[4][4] , double b[4][4] ){
	double** c = new double*[4];
	for(int i =0; i< 4; i++)
		c[i] = new double[4];
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			double sum = 0;
			for(int k = 0; k < 4; k++)
			{
				sum+= a[i][k]*b[k][j];
			}
			c[i][j] = sum;
		}

	}
	return c;
}
double** IKSimulator::getShoulder_X(double t1){
		double** c = new double*[4];
		for(int i =0; i< 4; i++)
			c[i] = new double[4];
		setIdentMat((double*)c, 4);
		double temp[4][4] = {{1,0,0,0}, {0,cos(t1), -sin(t1), 0}, {0, sin(t1), cos(t1), 0}, {0,0,0,1}};
		return c;
	}
void IKSimulator::getColumnJacobian(int i){
	double temp_c[3];
	switch(i){
	case 1:
		double ** wr = getWristMatrix(m_guy.angle[5], m_guy.angle[6]);	//theta 6, theta 7
		double ** el = getElbowMatrix(m_guy.angle[3], m_guy.angle[4]); //theta 4, theta 5
		double **wrXel = MatrixMult4x4((double(*)[4])wr, (double(*)[4])el);
		double sh_x_theta1[4][4] = {{0,0,0,0},
									{0,-sin(m_guy.angle[0]), -cos(m_guy.angle[0]),0},
									{0, cos(m_guy.angle[0]), -sin(m_guy.angle[0]),0},
									{0,0,0,0}};
		double **wrXelX
		break;
	case 2:
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	}
	m_J[0][i] = temp_c[0];
	m_J[1][i] = temp_c[1];
	m_J[2][i] = temp_c[2];
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
    glTranslatef(0,-1.0,0);
    glScalef(0.75,2,0.75);
    glutSolidSphere(.6,10,10);
	glEndList();
	glPopMatrix();
	//Upper Arm
	glPushMatrix();
	glNewList(m_guy.list[ARM_U], GL_COMPILE);
    glTranslatef(0,-1.2,0);
    glScalef(.75,2,.75);
    glutSolidSphere(.65,10,10);
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
    glTranslatef(0,-HAND_LEN,0);
    glScalef(1,2,1);
    glutSolidSphere(.25,10,10);
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
  glPushMatrix();
    glCallList(m_guy.list[ARM_L]);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(0,-ARM_L_LEN,0);
	glRotatef(m_guy.angle[6], 0, 0, 1);
    glRotatef(m_guy.angle[5], 0, 1, 0);
    draw_hand();
  glPopMatrix();
}  

void IKSimulator::draw_arm()
{
  glPushMatrix();
    glCallList(m_guy.list[ARM_U]);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(0,-ARM_U_LEN, 0);
    glRotatef(m_guy.angle[4], 0, 1, 0);
    glRotatef(m_guy.angle[3], 1, 0, 0);
    draw_lowerarm();
  glPopMatrix();
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
		glTranslatef(0,-2.5,0);
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

	//Left Arm

	glPushMatrix();
	glTranslatef(-SHOULDER_ROOT_LEN_X,SHOULDER_ROOT_LEN_Y,0);
		glRotatef(m_guy.angle[2], 0, 0, 1);
		glRotatef(m_guy.angle[1], 0, 1, 0);
		glRotatef(m_guy.angle[0], 1, 0, 0);
    draw_arm();
	glPopMatrix();

	//Right Arm
	glPushMatrix();
	glRotatef(180, 0, 1, 0);
	glPushMatrix();
		
			glTranslatef(-SHOULDER_ROOT_LEN_X,SHOULDER_ROOT_LEN_Y,0);
			glRotatef(m_guy.angle[2], 0, 0, 1);
			glRotatef(m_guy.angle[1], 0, 1, 0);
			glRotatef(m_guy.angle[0], 1, 0, 0);
			glPushMatrix();
				glCallList(m_guy.list[ARM_U]);
			glPopMatrix();
			glPushMatrix();
				glTranslatef(0,-ARM_U_LEN, 0);
				glPushMatrix();
					glCallList(m_guy.list[ARM_L]);
				glPopMatrix();
				glPushMatrix();
					glTranslatef(0,-ARM_L_LEN,0);	
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
	glVertex3f(10,-8,10);
	glVertex3f(10,-8,-10);
	glVertex3f(-10,-8,-10);
	glVertex3f(-10,-8,10);
	glEnd();
	
	//Blackboard
	glColor3f(0.0,0.0,0.0);
	glBegin(GL_POLYGON);
	glVertex3f(8,-4,0);
	glVertex3f(8,5,0);
	glVertex3f(-8,5,0);
	glVertex3f(-8,-4,0);
	glEnd();

	//Wall
	glColor3f(0.545098,0.270588,0.0745098);
	glBegin(GL_POLYGON);
	glVertex3f(10,-5,0);
	glVertex3f(10,6,0);
	glVertex3f(-10,6,0);
	glVertex3f(-10,-5,0);
	glEnd();
	glPopMatrix();
	glEnable(GL_LIGHTING);
	glPushMatrix();
	
	//Draw human
	glTranslatef(0,-1,m_dist_from_board);
	draw_human();
	glPopMatrix();
	

	
	//test path
	/*
	VectorObj temp_p;
	glPushMatrix();
	glTranslatef(m_p.x(), m_p.y(),m_p.z());
	glutSolidSphere(1.0,10,10);
	glPopMatrix();
	glPopAttrib();
	*/
}
