#include <bits/stdc++.h>
#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif 

#define pi (2*acos(0.0))

#define cubSphr 1


using namespace std;

const int INF = 2999999999;

const double PI = acos(-1.0);

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
int tsk;

struct Point
{
	double x,y,z;
	Point() {}
	Point(double x, double y, double z) : x(x), y(y), z(z) {}
	Point(const Point &p) : x(p.x), y(p.y), z(p.z)	{}
	
};


int dirX[]={1,-1,-1,1,1,-1,-1,1};
int dirY[]={1,1,-1,-1,1,1,-1,-1};
int dirZ[]={1,1,1,1,-1,-1,-1,-1};


Point cmPos,cmUp,cmLk,cmRt;
double cubeLen,sphereRad;
void drwCubSphr(double a,double r);

void drawAxes()
{
    glPushMatrix();{
	if(drawaxes==1)
	{
		glBegin(GL_LINES);{
            glColor3f(1.0, 0, 0);
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

            glColor3f(0, 1.0, 0);
			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

            glColor3f(0, 0, 1.0);
			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
    }glPopMatrix();
}

void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		double a=30;
		glBegin(GL_LINES);{
			for(i=-a;i<=a;i++){
				//lines parallel to Y-axis
				glVertex3f(i*10, -10*a, 0);
				glVertex3f(i*10,  10*a, 0);

				//lines parallel to X-axis
				glVertex3f(-10*a, i*10, 0);
				glVertex3f( 10*a, i*10, 0);
			}
		}glEnd();
	}
}

void drwSqr(double a)
{
    glColor3f(1.0,0.0,0.0);
    a/=2;
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}

void drwSphrQrtr(double radius,int slices,int stacks,int quarter) //quarter is 0-based
{
    glPushMatrix();{
    glRotatef(90*(quarter%4),0,0,1);
    int rev=1;
    if(quarter>=4)  rev=-1;
	struct Point pnts[stacks+2][slices+2];
	int i,j;
	double h,r;
	for(i=0;i<=stacks;i++)
	{
	    double ang=((double)i/(double)stacks)*(pi/2);
		h=radius*sin(ang);
		r=radius*cos(ang);
		for(j=0;j<=slices;j++)
		{
		    double ang2=((double)j/(double)slices)*pi/2;
			pnts[i][j].x=r*cos(ang2);
			pnts[i][j].y=r*sin(ang2);
			pnts[i][j].z=h*rev;
		}
	}
	//draw quads using generated pnts
	for(i=0;i<stacks;i++)
	{

        double shade;
        if(i<stacks/2)  shade=((double)i/(double)stacks);
        else  shade=(1-(double)i/(double)stacks);

        glColor3f(shade,shade,shade);

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
				glVertex3f(pnts[i][j].x,pnts[i][j].y,pnts[i][j].z);
				glVertex3f(pnts[i][j+1].x,pnts[i][j+1].y,pnts[i][j+1].z);
				glVertex3f(pnts[i+1][j+1].x,pnts[i+1][j+1].y,pnts[i+1][j+1].z);
				glVertex3f(pnts[i+1][j].x,pnts[i+1][j].y,pnts[i+1][j].z);

			}glEnd();
		}
	}
	}glPopMatrix();
}

void drwSphr(double radius,int slices,int stacks)
{
    for(int i=0;i<8;i++)    drwSphrQrtr(radius,slices,stacks,i);
}

void drwRctngl(double width,double height)
{
    //glColor3f(1.0,0.0,0.0);
    glPushMatrix();{
    double x=width/2,y=height/2;
//glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);{
        glVertex3f(x,y,0);
        glVertex3f(x,-y,0);
        glVertex3f(-x,-y,0);
        glVertex3f(-x,y,0);
    }glEnd();
    }glPopMatrix();
}


void drwQrtrCylndr(double height,double radius,int sgmnt,int quarter)   //half cylinder is on positive z-axis, another half is on the negative.
{
    glPushMatrix();{
    glRotatef(90*quarter,0,0,1);
    Point pnts[sgmnt+2];
    int i;
    height/=2;
    for(i=0;i<=sgmnt;i++) {
        double ang=((double)i/(double)sgmnt)*pi/2;
        pnts[i].x=radius*cos(ang);
        pnts[i].y=radius*sin(ang);
        pnts[i].z=height;
    }
    double shade;
    for(i=0;i<sgmnt;i++) {
               
            glBegin(GL_TRIANGLES);{
            glVertex3f(0,0,height);
            glVertex3f(pnts[i].x,pnts[i].y,pnts[i].z);
            glVertex3f(pnts[i+1].x,pnts[i+1].y,pnts[i+1].z);
        }glEnd();

        glBegin(GL_QUADS);{
            glVertex3f(pnts[i].x,pnts[i].y,pnts[i].z);
            glVertex3f(pnts[i+1].x,pnts[i+1].y,pnts[i+1].z);
            glVertex3f(pnts[i+1].x,pnts[i+1].y,-pnts[i+1].z);
            glVertex3f(pnts[i].x,pnts[i].y,-pnts[i].z);
        }glEnd();
    }

    }glPopMatrix();
}

void drwCubSphr(double a,double r)
{
    int i;
    int sgmnt=20;
    glPushMatrix(); {
        a=a/2;
        a=a-r;
      
        /****QUARTER-SPHERES****/
        for(i=0;i<8;i++)    {

            glPushMatrix();{
                glTranslatef(a*dirX[i],a*dirY[i],a*dirZ[i]);
                drwSphrQrtr(r,sgmnt,sgmnt,i);
            }glPopMatrix();
        }

        /***quarter cylinders***/


        for(i=0;i<=2;i++)    {
            if(i==1)   glRotatef(90,1,0,0);
            if(i==2)   glRotatef(90,0,1,0);
            for(int j=0;j<4;j++)    {
                glPushMatrix();{
                    glTranslatef(a*dirX[j],a*dirY[j],0);
                    drwQrtrCylndr(a*2,r,sgmnt,j);
                }glPopMatrix();
            }
            if(i==2)   glRotatef(-90,0,1,0);
            if(i)   glRotatef(-90,1,0,0);
        }
        a+=r;

        /***side walls***/
        glPushMatrix(); {
            glColor3f(1,1,1);

            double sqLen=2*(a-r);
            for(i=0;i<4;i++)    {
                glPushMatrix();{
                    glRotatef(90*i,0,0,1);
                    glTranslatef(a,0,0);
                    glRotatef(90,0,1,0);

                    drwSqr(sqLen);
                }glPopMatrix();
            }
            glTranslatef(0,0,a);
            drwSqr(sqLen);
            glTranslatef(0,0,-2*a);
            drwSqr(sqLen);
        }glPopMatrix();

    }
    glPopMatrix();
}

Point crsPrdct(const Point &a,const Point &b)
{
    Point rt;
    rt.x=a.y*b.z-a.z*b.y;
    rt.y=-(a.x*b.z-a.z*b.x);
    rt.z=a.x*b.y-a.y*b.x;
    return rt;
}

/***Keyboard and mouse listeners***/
void keyboardListener(unsigned char key, int x,int y)
{

    if(tsk==cubSphr)  {
        Point l=cmLk,r=cmRt,u=cmUp;
        double ang=.05;// 3 degree

        switch(key) {
            case '1':   
                ang*=-1;
                cmRt.x=r.x*cos(ang)+l.x*sin(ang);
                cmRt.y=r.y*cos(ang)+l.y*sin(ang);
                cmRt.z=r.z*cos(ang)+l.z*sin(ang);

                cmLk=crsPrdct(cmUp,cmRt);

                break;
            case '2':
                cmRt.x=r.x*cos(ang)+l.x*sin(ang);
                cmRt.y=r.y*cos(ang)+l.y*sin(ang);
                cmRt.z=r.z*cos(ang)+l.z*sin(ang);

                cmLk=crsPrdct(cmUp,cmRt);
                break;
            case '3':       
                cmLk.x=l.x*cos(ang)+u.x*sin(ang);
                cmLk.y=l.y*cos(ang)+u.y*sin(ang);
                cmLk.z=l.z*cos(ang)+u.z*sin(ang);

                cmUp=crsPrdct(cmRt,cmLk);
                break;
            case '4':       
                ang*=-1;
                cmLk.x=l.x*cos(ang)+u.x*sin(ang);
                cmLk.y=l.y*cos(ang)+u.y*sin(ang);
                cmLk.z=l.z*cos(ang)+u.z*sin(ang);

                cmUp=crsPrdct(cmRt,cmLk);
                break;

            case '5':       
                cmUp.x=u.x*cos(ang)+r.x*sin(ang);
                cmUp.y=u.y*cos(ang)+r.y*sin(ang);
                cmUp.z=u.z*cos(ang)+r.z*sin(ang);

                cmRt=crsPrdct(cmLk,cmUp);

                break;
            case '6':
                ang*=-1;
                cmUp.x=u.x*cos(ang)+r.x*sin(ang);
                cmUp.y=u.y*cos(ang)+r.y*sin(ang);
                cmUp.z=u.z*cos(ang)+r.z*sin(ang);

                cmRt=crsPrdct(cmLk,cmUp);
                break;
        }

    }

}

void specialKeyListener(int key, int x,int y)
{
    if(tsk==cubSphr)   {
        switch(key) {
            case GLUT_KEY_F3:
                if(sphereRad<cubeLen/2) sphereRad++;
                break;
            case GLUT_KEY_F4:
                if(sphereRad)   sphereRad--;
                break;
            case GLUT_KEY_PAGE_DOWN:		//down arrow key
                cmPos.x -= 3*cmUp.x;
                cmPos.y -= 3*cmUp.y;
                cmPos.z -= 3*cmUp.z;
                break;
            case GLUT_KEY_PAGE_UP:		// up arrow key
                cmPos.x += 3*cmUp.x;
                cmPos.y += 3*cmUp.y;
                cmPos.z += 3*cmUp.z;
                break;

            case GLUT_KEY_RIGHT:
                cmPos.x+=3*cmRt.x;
                cmPos.y+=3*cmRt.y;
                cmPos.z+=3*cmRt.z;
                break;
            case GLUT_KEY_LEFT:
                cmPos.x-=3*cmRt.x;
                cmPos.y-=3*cmRt.y;
                cmPos.z-=3*cmRt.z;
                break;

            case GLUT_KEY_UP:
                cmPos.x+=3*cmLk.x;
                cmPos.y+=3*cmLk.y;
                cmPos.z+=3*cmLk.z;
                break;
            case GLUT_KEY_DOWN:
                cmPos.x-=3*cmLk.x;
                cmPos.y-=3*cmLk.y;
                cmPos.z-=3*cmLk.z;
                break;
        }

    }
}

void mouseListener(int button, int state, int x, int y)
{	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			cameraHeight--;
            break;

	
		default:
			break;
	}
}


/***main display function***/
void display()
{

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	
	if(tsk==cubSphr)   {
          gluLookAt(cmPos.x,cmPos.y,cmPos.z,cmPos.x+5*cmLk.x,cmPos.y+5*cmLk.y,cmPos.z+5*cmLk.z,cmUp.x,cmUp.y,cmUp.z);
	}

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	glPushMatrix();{
	//add objects

	drawAxes();
	drawGrid();

    tsk=cubSphr;
    switch(tsk)    {
        case cubSphr:
            drwCubSphr(cubeLen,sphereRad);
            break;
    }

    }glPopMatrix();
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate()
{
	
	glutPostRedisplay();
}

void init()
{
	
	tsk=cubSphr;

	/** Camera initialization **/
        drawgrid=0;
	drawaxes=1;
	cameraHeight=80;
	cameraAngle=pi/4;

    cmPos=Point(-90,-90,10);
    cmLk=Point(1/sqrt(3.0),1/sqrt(3.0),0);
    cmUp=Point(0,0,1);
    cmRt=Point(1/sqrt(3.0),-1/sqrt(3.0),0);
	

    /** Cube-Sphere initialization **/
    cubeLen=50;
    sphereRad=10;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,1,1,1000.0);
	
}




int main(int argc, char **argv)
{
	glutInit(&argc,argv);
	glutInitWindowSize(800,800);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("***Cube-Sphere Transformation***");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
