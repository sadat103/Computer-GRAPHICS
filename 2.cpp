#include <bits/stdc++.h>
#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif 

#define pi (2*acos(0.0))

using namespace std;

#define WHEEL 4

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


double whlAng,whlRd,whlDst;
Point whlPos;
void drwWhl(Point pos,double ang,double dist,double rad);

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


void drwCrcl(double radius,int sgmnts)
{
    int i;
    struct Point points[100];
    glColor3f(0.0,0.5,0.0);
    //generate points
    for(i=0;i<=sgmnts;i++)
    {
        points[i].x=radius*cos(((double)i/(double)sgmnts)*2*pi);
        points[i].y=radius*sin(((double)i/(double)sgmnts)*2*pi);
    }
    //draw sgmnts using generated points
    for(i=0;i<sgmnts;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
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

void drwSqr(double a)
{
    glColor3f(1.0,0.0,0.0);
    a=a/2;
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}

void drwWhl(Point pos,double ang,double dist,double rad)
{
    int sgmnts=40;
    int i,j,k,l;
    double prmtr=2*pi*rad;
    double rtAng = 360*dist/prmtr;
    glPushMatrix();
    glTranslatef(0,0,rad);
    glTranslatef(pos.x,pos.y,pos.z);
    glRotatef(ang,0,0,1);
    glRotatef(rtAng,0,1,0); //rotation Y-axis
    double shd;

    Point point[101];

    for(i=0;i<=sgmnts;i++) {
        double crclAng=2*pi*(double)i/(double)sgmnts;   //radian
        point[i].x=rad*cos(crclAng);
        point[i].y=0;
        point[i].z=rad*sin(crclAng);
        
    }
    i=0;
    double sqLen=prmtr/sgmnts;

    for(i=0;i<sgmnts;i++) {
        if(i<sgmnts/2)    shd=2*((double)i/(double)sgmnts);
        else    shd=2*(1-(double)i/(double)sgmnts);
        double crclAng=360.0*(double)i/(double)sgmnts;    //degree
        glPushMatrix(); {
            glColor3f(shd,shd,shd);
                glRotatef(crclAng,0,1,0);
                glTranslatef(rad,0,0);
                glRotatef(90,0,1,0);
                drwRctngl(sqLen,4*sqLen);

        }
        glPopMatrix();
    }
    glColor3f(1.0,1.0,0.0);
    drwRctngl(2*whlRd,sqLen);
    glColor3f(1.0,0.5,0.0);
    glRotatef(90,0,1,0);
    drwRctngl(2*whlRd,sqLen);

    glPopMatrix();
}
/***Keyboard and mouse listeners***/
void keyboardListener(unsigned char key, int x,int y)
{
    if(tsk==WHEEL)    {
        double whlRdian = whlAng*pi/180.0;
        switch(key){
            case 'w':
                whlPos.x+=3*cos(whlRdian);   //3 steps
                whlPos.y+=3*sin(whlRdian);   //3steps
                whlDst+=3;
                break;
            case 's':
                whlPos.x-=3*cos(whlRdian);   //3steps
                whlPos.y-=3*sin(whlRdian);   //3 steps
                whlDst-=3;
                break;
            case 'a':
                whlAng+=3;
                break;
            case 'd':
                whlAng-=3;
                break;
            default:
                break;
        }
    }
}

void specialKeyListener(int key, int x,int y)
{
    switch(key){
        case GLUT_KEY_DOWN:		
            cameraHeight -= 5.0;
            break;
        case GLUT_KEY_UP:		
            cameraHeight += 5.0;
            break;

        case GLUT_KEY_RIGHT:
            cameraAngle += 0.05;
            break;
        case GLUT_KEY_LEFT:
            cameraAngle -= 0.05;
            break;

    

        default:
            break;
    }
}

void mouseListener(int button, int state, int x, int y)
{	
	switch(button){
		case GLUT_LEFT_BUTTON:
			
			break;

		case GLUT_RIGHT_BUTTON:
			
			break;

		case GLUT_MIDDLE_BUTTON:
			
			break;

        case 3:             //scroll up
            cameraHeight--;
            break;
        case 4:             //scroll down
            cameraHeight++;
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
	glClearColor(0,0,0,0);	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_MODELVIEW);

	
	glLoadIdentity();
        gluLookAt(-150*cos(cameraAngle), -150*sin(cameraAngle), cameraHeight,	0,0,0,	0,0,1);


	glMatrixMode(GL_MODELVIEW);

	
	glPushMatrix();{
	

	drawAxes();
	drawGrid();

    tsk=WHEEL;
    switch(tsk)    {
        case WHEEL:
            drawgrid=1;
            drawaxes=0;
            drwWhl(whlPos,whlAng,whlDst,whlRd);
            break;
    }

    }glPopMatrix();
	
	glutSwapBuffers();
}

void animate()
{
	
	glutPostRedisplay();
}

void init()
{
      tsk=WHEEL;
      
      whlRd=50;
      whlAng=0;
      whlDst=0;
      whlPos.x=0;
      whlPos.y=0;
      whlPos.z=0;
	
       drawgrid=0;
       drawaxes=1;
      cameraHeight=100;
      cameraAngle=pi/4;

       glClearColor(0,0,0,0);

	
	glMatrixMode(GL_PROJECTION);

	
	glLoadIdentity();

	
	gluPerspective(80,1,1,1000.0);
	
}

int main(int argc, char **argv)
{
	glutInit(&argc,argv);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(10,0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	

	glutCreateWindow("**Wheel Simulation**");

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

