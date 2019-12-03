/*
Assignment: OpenGL Drawing-1 (Cube-Sphere Transformation)
Author: Amlan Saha
*/
#include <bits/stdc++.h>
#include "bitmap_image.hpp"

#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif // windows

#define pi (2*acos(0.0))

using namespace std;

typedef long long LL;
typedef unsigned long long LLU;
typedef vector<int> VI;
typedef vector<long long> VLL;
typedef map<int, int> MAPII;
typedef map<string,int> MAPSI;
typedef pair<int, int> PII;
typedef pair<LL, LL> PLL;
typedef pair<double, double> PDD;

#define FOR(i,a,b) for(i=a;i<=b;i++)
#define ROF(i,a,b) for(i=a;i>=b;i--)
#define FR(i,n)    for(i=0;i<n;i++)
#define RF(i,n) for(i=n;i>0;i--)
#define CLR(a) memset ( a, 0, sizeof ( a ) )
#define RESET(a) memset ( a, -1, sizeof ( a ) )
#define PB(a)    push_back ( a )
#define XX first
#define YY second

//#define X v[0]
//#define Y v[1]
//#define Z v[2]

#define CAMERA_ANGLE 1
#define CUBE_SPHERE 1
#define TORUS 2
#define HAND 3
#define WHEEL 4
#define SPHERICON 5

const int INF = 2000000009;
const int Max = 1000007;
const double PI = acos(-1.0);
const double EPS = 1e-10;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
int TASK;

//struct Point
//{
//	double x,y,z;
//	Point() {}
//	Point(double x, double y, double z) : x(x), y(y), z(z) {}
//	Point(const Point &p) : x(p.x), y(p.y), z(p.z)	{}
//	void print(){
//        cout<< "(" << x << "," << y << ","<<z<<")"<<endl;
//	}
//};

struct Vector
{
    double X,Y,Z;
    int dimension;
    Vector(){}
    Vector(double ii,double jj,double kk)   {
        X=ii;
        Y=jj;
        Z=kk;
        dimension=3;
    }
    Vector(const Vector &vv)   {
//        for(int i=0;i<dimension;i++)    v[i]=vv.v[i];
        X = vv.X;
        Y = vv.Y;
        Z = vv.Z;
    }
    bool operator == (const Vector &vv) const  {
        if(abs(X-vv.X)>EPS)   return false;
        if(abs(Y-vv.Y)>EPS)   return false;
        if(abs(Z-vv.Z)>EPS)   return false;

        return true;
    }
    Vector operator *(const double &x)   const {
        Vector ret;
//        for(int i=0;i<dimension;i++)    ret.v[i]=v[i]*x;
        ret.X=X*x;
        ret.Y=Y*x;
        ret.Z=Z*x;
//        ret.v[dimension]=1;
        return ret;
    }
    Vector operator *(const Vector& x)   const { //cross product
        Vector ret;
        ret.X=Y*x.Z-Z*x.Y;
        ret.Y=-(X*x.Z-Z*x.X);
        ret.Z=X*x.Y-Y*x.X;
//        ret.v[3]=1;
//        for(int i=0;i<dimension;i++)    ret.v[i]=v[i]*x;
        return ret;
    }
    Vector operator +(const Vector &v)  const {
        Vector ret;
        ret.X = X+v.X;
        ret.Y = Y+v.Y;
        ret.Z = Z+v.Z;
        return ret;
    }
    Vector operator -(const Vector &v)  const {
        Vector ret;
        ret.X = X-v.X;
        ret.Y = Y-v.Y;
        ret.Z = Z-v.Z;
        return ret;
    }
    double absVal() {
//        if(v[dimension]>1+EPS)  {
//            for(int i=0;i<dimension;i++)    v[i]/=v[dimension];
//            v[dimension]=1;
//        }
        return sqrt(X*X+Y*Y+Z*Z);
    }
    void normalize()    {
        double val=this->absVal();
        if(val<EPS) return;
        X/=val;
        Y/=val;
        Z/=val;
    }
//    void print(ostream &ff = cout,int precision=7)    {
    void print()    {
//        const ostream &fout=*ff;
        int precision = 7;
        cout<<fixed<<setprecision(precision)<<X<<" "<<Y<<" "<<Z<<" "<<endl;
//        fprintf(fp,"%.7lf %.7lf %.7lf\n",X,Y,Z);
    }
};

typedef Vector Point;
//#define Point Vector


void keyboardListener(unsigned char key, int x,int y);
void specialKeyListener(int key, int x,int y);
void mouseListener(int button, int state, int x, int y);
void animate();
void init();

void drawAxes();
void drawGrid();
void drawSquare(double a);
void drawRectangle(double width,double height);
void drawCircle(double radius,int segments);
void drawSphereQuarter(double radius,int slices,int stacks,int quarter);
void drawSphere(double radius,int slices,int stacks);
void drawQuarterCylinder(double height,double radius,int segments,int quarter);
void printImage();
void rayTracer();

double dot(const Vector &a, const Vector &b)
{
    double ret = a.X*b.X+a.Y*b.Y+a.Z*b.Z;
    return ret;
}

inline double matrixDeterminant(const Vector &col1, const Vector &col2, const Vector &col3)
{
    return dot(col1, (col2*col3));
}

inline double triangleArea(const Vector &col1, const Vector &col2, const Vector &col3)
{
    Vector a = col1-col3;
    Vector b = col2-col3;
    double area = .5*sqrt(dot(a*b, a*b));
    return area;
}

struct Ray
{
    Point p0;
    Vector dir;
    Ray()  {}
    Ray(Point p0, Point p1)    {
        dir = p1+(p0*-1);
        dir.normalize();
    }
//    Ray(const Ray& cl)    {
//        p0=cl.p0;
//        dir=cl.dir;
//    }
};

struct Plane
{
    double a,b,c,d;
    Plane(double a=0,double b=0,double c=0,double d=0): a(a),b(b),c(c),d(d){}
};

struct Color
{
    double r,g,b;
    Color() {}
    Color(double r,double g,double b):   r(r),g(g),b(b){}
    Color(const Color &c):   r(c.r),g(c.g),b(c.b){}
    Color operator *(const double &x)   const {
        Color ret;
//        for(int i=0;i<dimension;i++)    ret.v[i]=v[i]*x;
        ret.r=r*x;
        ret.g=g*x;
        ret.b=b*x;
//        ret.v[dimension]=1;
        return ret;
    }
    Color operator +(const Color &v)  const {
        Color ret;
        ret.r = r+v.r;
        ret.g = g+v.g;
        ret.b = b+v.b;
        return ret;
    }
    Color operator -(const Color &v)  const {
        Color ret;
        ret.r = r-v.r;
        ret.g = g-v.g;
        ret.b = b-v.b;
        return ret;
    }
    Color dotProduct(const Color &c)    {
        Color ret(r*c.r,g*c.g,b*c.b);
        return ret;
    }
    bool operator == (const Color &c)   const   {
        return abs(r-c.r)<EPS && abs(g-c.g)<EPS && abs(b-c.b)<EPS;
    }
    void print()    {
        cout<<"("<<r<<","<<g<<","<<b<<")"<<endl;
    }
};

struct Triangle
{
    Point p[3];
    Color c;
    double amb, diffuse, spec, refl;
    double shine;

    Triangle()  {}
    Triangle(Point p0,Point p1,Point p2)    {
        p[0]=p0;
        p[1]=p1;
        p[2]=p2;
    }
    Triangle(Point p0,Point p1,Point p2,Color colr)    {
        p[0]=p0;
        p[1]=p1;
        p[2]=p2;
        c=colr;
    }
    Triangle(const Triangle& t) {
        p[0]=t.p[0];
        p[1]=t.p[1];
        p[2]=t.p[2];
        c=t.c;
        amb = t.amb;
        diffuse = t.diffuse;
        spec = t.spec;
        shine = t.shine;
        refl = t.refl;
    }

    bool isOnTriangle(const Point &pp)   {
        double tot = 0;
        for(int i = 0; i < 3; i++)  tot+= triangleArea(p[i],p[(i+1)%3],pp);
        double ar = triangleArea(p[0],p[1],p[2]);
        if(abs(tot-ar)<EPS) return true;
        return false;
    }

    double rayIntersection(Ray ri, Vector &normal) {       //returns the parametric value of t
                                                            //where (ri.p0+ri.dir*t) is the point of intersection
                                                            //returns negative if no intersection
                                                            //if valid intersection, the direction of the normal is stored
        Vector ab = p[1]-p[0];
        Vector ac = p[2]-p[0];
        ri.dir.normalize();

        double det = matrixDeterminant(ri.dir, ab,ac);
        if(abs(det)<EPS)    {   //co-plane rays
            ///TO BE IMPLEMENTED
            return -INF;
        }

        Vector cons = ri.p0-p[0];
        double numerator = matrixDeterminant(ab,ac,cons);
        double denom = matrixDeterminant(ab,ac,ri.dir*-1);
        if(abs(denom)<EPS)  {       //no intersection
            return -INF;
        }
        double t = numerator/denom;   //cramers rule

        bool onTrngl = isOnTriangle(ri.p0+ri.dir*t);
        if(onTrngl==false)  return -INF;
        normal = ab*ac;
        normal.normalize();
        return t;
    }

    void draw() {
        glColor3f(c.r,c.g,c.b);
        glBegin(GL_TRIANGLES);{
            glVertex3f(p[0].X,p[0].Y,p[0].Z);
            glVertex3f(p[1].X,p[1].Y,p[1].Z);
            glVertex3f(p[2].X,p[2].Y,p[2].Z);
        }glEnd();

    }

    void read(istream &fin = cin)  {
        for(int i = 0; i<3; i++)    cin>>p[i].X>>p[i].Y>>p[i].Z;
    }
};

struct Sphere
{
    Point center;
    double rad;
    Color c;
    double amb, diffuse, spec, refl;
    double shine;

    Sphere()  {}
    Sphere(Point cen, double radius, Color colr=Color(0,0,0))    {
        center = cen;
        rad = radius;
        c=colr;
    }
    Sphere(const Sphere& t) {
        center = t.center;
        rad = t.rad;
        c=t.c;
        amb = t.amb;
        diffuse = t.diffuse;
        spec = t.spec;
        shine = t.shine;
        refl = t.refl;
    }

    double rayIntersection(Ray ri, Vector &normal) {       //returns the parametric value of t
                                                            //where (ri.p0+ri.dir*t) is the point of intersection
                                                            //returns negative if no intersection
                                                            //if valid intersection, the direction of the normal is stored
        ri.dir.normalize();
        double b = 2*dot(ri.dir, (ri.p0-center));
        double k = dot(ri.p0-center, ri.p0-center)-rad*rad;
        double a = dot(ri.dir,ri.dir);
        double tmp = b*b-4*a*k;
        if(tmp<0)   return -INF;    //no intersection
        tmp = sqrt(tmp);
        double t1 = (-b+tmp)/(2*a);
        double t2 = (-b-tmp)/(2*a);
//        cout<<"TTTT: "<<t1<<" "<<t2<<" rad: "<<rad<<" center: ";
//        center.print();
        if(t1<0)    {
            t1 = INF;
            swap(t1,t2);
        }
        if(t1<0)    {
            return -INF;    //no intersection
        }
        double t = min(t1,t2);
        normal = (ri.p0+(ri.dir)*t)-center;
        normal.normalize();
        return t;
    }

    void draw() {
        glPushMatrix(); {
            glColor3f(c.r,c.g,c.b);
            glTranslatef(center.X, center.Y,center.Z);
            drawSphere(rad, 20, 20);
        }
        glPopMatrix();
    }

    void read()  {
        cin>>center.X>>center.Y>>center.Z;
//        center.print();
        cin>>rad;
//        cout<<rad<<endl;
        cin>>c.r>>c.g>>c.b;
//        cout<<"COLLLL: "<<c.r<<" "<<c.g<<" "<<c.b<<endl;
        cin>>amb>>diffuse>>spec>>refl;
        cin>>shine;
//        cout<<"shine: "<<shine<<endl;
    }
};

struct Square
{
    Point p[4];
    Color c;
    double amb, diffuse, spec,refl;
    double shine;

    Square()  {}
    Square(Point lowerLeft, double length, Color c)   {
        p[0] = lowerLeft;
        p[1] = p[0]+Point(length,0,0);
        p[2] = p[1]+Point(0,length,0);
        p[3] = p[0]+Point(0,length,0);
        this->c = c;
    }
    Square(Point p0,Point p1,Point p2,Point p3,Color colr)    {
        p[0]=p0;
        p[1]=p1;
        p[2]=p2;
        p[3]=p3;
        c=colr;
    }
    Square(const Square& t) {
        p[0]=t.p[0];
        p[1]=t.p[1];
        p[2]=t.p[2];
        p[3]=t.p[3];
        c=t.c;
        amb = t.amb;
        diffuse = t.diffuse;
        spec = t.spec;
        refl = t.refl;
        shine = t.shine;
    }

    double rayIntersection(const Ray &ri, Vector &normal)   {
        Triangle tr(p[0],p[1],p[2]);
        Triangle tr2(p[0],p[3],p[2]);
        double t = tr.rayIntersection(ri,normal);
        if(t>=0)    {
            return t;
        }
        return tr2.rayIntersection(ri,normal);
    }

    void draw() {
        glColor3f(c.r,c.g,c.b);
        glBegin(GL_QUADS);{
            glVertex3f(p[0].X,p[0].Y,p[0].Z);
            glVertex3f(p[1].X,p[1].Y,p[1].Z);
            glVertex3f(p[2].X,p[2].Y,p[2].Z);
            glVertex3f(p[3].X,p[3].Y,p[3].Z);
        }glEnd();

    }
};

struct Pyramid
{
    Triangle t[4];
    Square sq;
    double width, height;
    Color c;

    double amb, diffuse, spec, refl;
    double shine;

    Pyramid()   {}
    Pyramid(Point lowLeft, double width, double height, Color colr) {
        Point base[4];
        this->width = width;
        this->height = height;
        base[0] = lowLeft;
        base[1] = lowLeft+Point(width,0,0);
        base[2] = lowLeft+Point(width,width,0);
        base[3] = lowLeft+Point(0,width,0);
        Point top = lowLeft+Point(width/2,width/2,height);
        sq = Square(base[0],base[1],base[2],base[3], colr);
        for(int i = 0; i < 4; i++)  {
            t[i] = Triangle(top,base[i],base[(i+1)%4],colr);
        }
        c = colr;
    }
    Pyramid (const Pyramid &t)    {
        for(int i = 0; i <4; i++)   this->t[i] = t.t[i];
        this->sq=t.sq;
        c=t.c;
        amb = t.amb;
        diffuse = t.diffuse;
        spec = t.spec;
        refl = t.refl;
        shine = t.shine;
    }


    double rayIntersection(const Ray &rri, Vector &normal)   {
        double t = INF;
        Ray ri = rri;
        ri.dir.normalize();

        for(int i = 0; i < 4; i++)  {
            Vector nm;
            double tmp = this->t[i].rayIntersection(ri,nm);
            if(tmp<0)   continue;
            if(tmp<t)   {
                t = tmp;
                normal = nm;
            }
        }
        Vector nm;
        double tmp = sq.rayIntersection(ri,nm);
        if(tmp>=0 && tmp<t)    {
            t = tmp;
            normal = nm;
        }
        if(t >INF-(EPS+EPS))    t = -INF;
        return t;
    }

    void draw() {
        glColor3f(c.r,c.g,c.b);
        for(int i = 0; i < 4; i++)  {
            t[i].draw();
        }
        sq.draw();
    }
};

struct CheckBoard
{
    Point center;
    double rad;
    Color c;
    double amb, diffuse, spec, refl;
    double shine;
    double height, width;
    double boardHeight, boardWidth;

    CheckBoard()  {}
    CheckBoard(double h, double w)    {
        height = h;
        width = w;
        boardHeight = 20;
        boardWidth = 20;
    }
    CheckBoard(const CheckBoard& t) {
        height=t.height;
        width=t.width;
        amb = t.amb;
        diffuse = t.diffuse;
        spec = t.spec;
        shine = t.shine;
        refl = t.refl;
        boardHeight = t.boardHeight;
        boardWidth = t.boardWidth;
    }

    bool onBoard(Point p)   {
        if(abs(p.Z)>EPS)    return false;
        if(-width/2<=p.X && p.X<=width/2)    {
            if(-height/2<=p.Y && p.Y<=height/2) return true;
        }
//        cout<<"PlanePoint: ";
//        p.print();
        return false;
    }

    Color getColor(Point p) {
        if(!onBoard(p)) return Color(0,0,0);
        Color ret;
        Color white(1,1,1),black(0,0,0);
        p=p+Point(width/2,height/2,0);
        int row = p.X/boardWidth, colm = p.Y/boardHeight;
        if((row+colm)%2==0) return white;
        return black;
    }
    void draw() {
        int row = 0;
        for(double xx = -width/2; xx<= width/2; xx+= boardWidth, row++)   {
            int colm = 0;
            for(double yy = -height/2; yy<=height/2; yy+= boardHeight, colm++)   {
                Color colr(1,1,1);
                if((row+colm)%2)    colr = Color(0,0,0);
                Square s(Point(xx,yy,0),boardHeight, colr);
                s.draw();
            }
        }
    }

    double rayIntersection(const Ray &ri, Vector &normal)   {
//        cout<< "kkasdf ray: ";{
//            Ray rr = ri;
//            rr.dir.print();
//        }
        if(abs(ri.dir.Z)<EPS)  return -INF;

        double t = -ri.p0.Z/ri.dir.Z;
        Point pt = ri.p0+(ri.dir*t);
        bool valid = onBoard(pt);
//        cout<<"LLLL: "<<valid<<" "<<t<<endl;
//        Ray rr = ri;
//        if(ri.dir.Z<-.97)   {
////            rr.p0.print();
////            rr.dir.print();
//            cout<<"PPQTTT: "<<t<<endl;
//            pt.print();
////            if(ri.dir.Z<0)  {
////            }
//        }
        if(!valid)  return -INF;
//        cout<<"VALID: "<<t<<"...";
//        rr.dir.print();
        normal = Vector(0,0,1);
        return t;
    }
};


int dirX[]={1,-1,-1,1,1,-1,-1,1};
int dirY[]={1,1,-1,-1,1,1,-1,-1};
int dirZ[]={1,1,1,1,-1,-1,-1,-1};

/***Prototypes for Cube-Sphere Transformation***/
Point cameraPos,cameraUp,cameraLook,cameraRight;
double cubeLen,sphereRad;
double screenHeight, screenWidth;
double recursionLevel;

void drawCubeSphere(double a,double r);

Point crossProduct(const Point &a,const Point &b)
{
    Point ret;
    ret.X=a.Y*b.Z-a.Z*b.Y;
    ret.Y=-(a.X*b.Z-a.Z*b.X);
    ret.Z=a.X*b.Y-a.Y*b.X;
    return ret;
}

Ray calculateReflection(const Ray &ri, const Vector &normal, const Point p0)
{
    Vector nm = normal;
    nm.normalize();
    Ray ret;
    ret.dir = ri.dir-(nm*(2*dot(ri.dir, nm)));
    ret.dir.normalize();
    ret.p0 = p0;
    return ret;
}

vector<Sphere>spheres;
vector<Pyramid>pyramids;
vector<Point>lights;
CheckBoard checkBoard;


Color imageMap[2002][2002];
Color sourcePower;

#define CHECKBOARD 0
#define PYRAMID 1
#define SPHERE 2
#define EYE 3

///the direction of the lightRay SHOULD NOT be normalized
bool hasObstacle(Ray lightRay)
{
    double tt = lightRay.dir.absVal();
    lightRay.dir.normalize();
    for(int i = 0; i < spheres.size(); i++) {
        Vector tmp;
        double t = spheres[i].rayIntersection(lightRay, tmp);
        if(t<0) {
            if(tt-EPS<t)    return true;
            continue;
        }
        else    {
            if(t<tt+EPS)    return true;
            continue;
        }
    }
    for(int i = 0; i < pyramids.size(); i++) {
        Vector tmp;
        double t = pyramids[i].rayIntersection(lightRay, tmp);
        if(t<0) {
            if(tt-EPS<t)    return true;
            continue;
        }
        else    {
            if(t<tt+EPS)    return true;
            continue;
        }
    }

    Vector tmp;
    double t = checkBoard.rayIntersection(lightRay, tmp);
    if(t<0) {
        if(tt-EPS<t)    return true;
    }
    else    {
        if(t<tt+EPS)    return true;
    }

    return false;
}

struct ObjectID
{
    int type;
    int id;
    ObjectID(int tp=EYE, int i=0) {type=tp; id=i;}
    ObjectID(const ObjectID &ob) {type=ob.type; id=ob.id;}
    bool operator == (const ObjectID &ob) const {
        return type==ob.type && id == ob.id;
    }
};

Color rayCast(Ray ri, int level, ObjectID from)
{
    if(level==0)    {
        return Color(0,0,0);
    }

    double spPt = INF, pyPt = INF;
    int sphr = -1, prmd=-1;
    Vector spNormal,pyNormal;
    double amb, spec, diffuse, reflCoeff,shine;

//    cout<<"******"<<endl;
    for(int i = 0; i < spheres.size(); i++) {
        if(from.type==SPHERE && from.id==i) continue;
        Vector nm;
        double tmp = spheres[i].rayIntersection(ri,nm);
        if(tmp<-1 || tmp>spPt)  continue;
//        cout<<"sp: "<<i<<" "<<spheres.size()<<" "<<tmp<<endl;
        sphr = i;
        spPt = tmp;
        spNormal = nm;
    }
    if(spPt>INF-EPS)    sphr= -1;
////    cout<<"______"<<endl;
    for(int i = 0; i < pyramids.size(); i++) {
        if(from.type==PYRAMID && from.id==i) continue;
        Vector nm;
        double tmp = pyramids[i].rayIntersection(ri,nm);
        if(tmp<-1 || tmp>=pyPt)  continue;
        prmd = i;
        pyPt = tmp;
        pyNormal = nm;
    }
    if(pyPt>INF-EPS)    prmd = -1;
    Vector nm;
    double chkPt = checkBoard.rayIntersection(ri,nm);
    if(from.type==CHECKBOARD)   {
        chkPt = -INF;
    }

    int selected;
    ObjectID newObj;

    Vector normal;
    Color colr;
    double minDist;

//    if(chkPt>=0)   {
//        cout<<"CHKLKAJDF: "<<chkPt<<" "<<spPt<<" "<<pyPt<<endl;
//    }
    if(sphr<0 && prmd<0)    {   //no intersection
        if(chkPt<0 || chkPt>INF-EPS) {
            return Color(0,0,0);
        }
//        cout<<"CHECKERBOARDAAAA : "<<chkPt<<endl;
        selected = CHECKBOARD;
        normal = Vector(0,0,1);
        minDist = chkPt;
    }

    else if(sphr<0 || pyPt<spPt)   { //a pyramid is closer
        selected = PYRAMID;
        minDist = pyPt;
        if(chkPt>=0 && chkPt<pyPt)    {
            selected = CHECKBOARD;
            normal = Vector(0,0,1);
            minDist = chkPt;
        }
    }
    else if(prmd<0 || spPt<=pyPt)    {   //a sphere is closer
        selected = SPHERE;
        minDist = spPt;
        if(chkPt>=0 && chkPt<spPt)  {
            selected = CHECKBOARD;
            normal = Vector(0,0,1);
            minDist = chkPt;
        }
    }
//    if(selected==PYRAMID) cout<<"selected: "<<selected<<" "<<sphr<<" "<< prmd<<" "<<pyPt<<" "<<chkPt<<endl;//" "<<amb<<" "<<diffuse<<endl;
//    cout<<"selected: "<<selected<<" "<<sphr<<" "<< prmd<<" "<<pyPt<<" "<<chkPt<<endl;//" "<<amb<<" "<<diffuse<<endl;
//    if(selected==CHECKBOARD)    return Color(0,0,0);
    Point incidentPoint = ri.p0+(ri.dir*minDist);

    if(selected == CHECKBOARD)   {
        colr = checkBoard.getColor(ri.p0+(ri.dir*minDist));
        normal = Vector(0,0,1);
        amb = .4;
        spec = .15, diffuse=.2,reflCoeff=.25;
        shine = 4;
        newObj = ObjectID(CHECKBOARD,0);
    }
    else if(selected == PYRAMID)    {
        colr = pyramids[prmd].c;
        amb = pyramids[prmd].amb;
        spec = pyramids[prmd].spec;
        diffuse = pyramids[prmd].diffuse;
        reflCoeff = pyramids[prmd].refl;
        shine = pyramids[prmd].shine;
        normal = pyNormal;
        newObj = ObjectID(PYRAMID, prmd);
    }
    else if(selected == SPHERE) {
        colr = spheres[sphr].c;
        amb = spheres[sphr].amb;
        spec = spheres[sphr].spec;
        diffuse = spheres[sphr].diffuse;
        reflCoeff = spheres[sphr].refl;
        shine= spheres[sphr].shine;
        normal = spNormal;
        newObj = ObjectID(SPHERE, sphr);
//        cout<<"PPP: "<<amb<<endl;
    }

//    if(selected!=SPHERE) cout<<"AMBB: "<<selected<<" "<<sphr<<" "<<pyPt<<" "<<chkPt<<" "<<amb<<" "<<diffuse<<endl;
    Ray rflct = calculateReflection(ri,normal,incidentPoint);

    Color ret(0,0,0);

    Color ambientLight = sourcePower.dotProduct(colr);
    Color diffuseLight = sourcePower.dotProduct(colr);
    Color specLight = sourcePower.dotProduct(colr);

    ret = (ambientLight*amb);

//    cout<<"INCIDENT POINTS: ";
//    incidentPoint.print();
//    rflct.dir.print();

    for(int i = 0; i < lights.size(); i++)  {
        Vector lt = incidentPoint-lights[i];
        Ray lightRay(lights[i],lt);

//        ///SHADOW
//        if(hasObstacle(lightRay))    continue;

        Ray lightRflct = calculateReflection(lightRay,normal,incidentPoint);
//        lt.print();
//        normal.print();

        double cosTheta = max(0.0,dot(lightRflct.dir,normal)/(lightRflct.dir.absVal()*normal.absVal()));
        double cosPhi = max(0.0,dot(ri.dir,lightRflct.dir)/(ri.dir.absVal()*normal.absVal()));
//        cout<<cosTheta<<" "<<cosPhi<<endl;
        Color df = diffuseLight*diffuse*cosTheta;
        Color spc = specLight*spec*pow(cosPhi,shine);
        ret = ret+df+spc;
//        ret = ret + (diffuseLight*diffuse*cosTheta) + (sourcePower*spec*pow(cosPhi,shine));
//        ret = ret+crnt;
    }

    Color lght = rayCast(rflct,level-1, newObj);
//    cout<<"refl col: ";lght.print();
    ret = ret+ (lght*reflCoeff);
//    cout<<"OIUOWEItuklansklgh asdv ";
//    ret.print();

    ret.r = min(ret.r,1.0);
    ret.g = min(ret.g,1.0);
    ret.b = min(ret.b,1.0);
    return ret;
}

void rayTracer()
{
    Ray src;
    src.p0 = cameraPos;
    sourcePower = Color(1,1,1);
//    screenHeight = 10;
//    screenWidth = 10;
    for(int row = 0; row<screenHeight; row++)   {
        for(int colm = 0; colm<screenWidth; colm++) {
            double ht = row-screenHeight/2;
            ht/=screenHeight/2;
            double wd = colm-screenWidth/2;
            wd/=screenWidth/2;
//            cout<<ht<<" "<<wd<<endl;
            Vector dir = (cameraUp*ht) + (cameraRight*wd);
            dir = cameraLook+dir;
//            dir.print();
            src.dir = dir;
//            cout<<"src dir:  ";
//            src.dir.print();
            imageMap[row][colm] = rayCast(src,recursionLevel, ObjectID(EYE,0));
//            if(std::isnan(imageMap[row][colm].r) || imageMap[row][colm]==Color(0,0,0)) continue;
//            cout<<"row: "<<row<<" column: "<<colm<<" ==> ";
//            imageMap[row][colm].print();
        }
    }
    cout<<"Image is ready: "<<endl;
}

/***main display function***/
void display()
{

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	if(TASK==CAMERA_ANGLE)   {
        gluLookAt(cameraPos.X,cameraPos.Y,cameraPos.Z,
                  cameraPos.X+10*cameraLook.X,cameraPos.Y+10*cameraLook.Y,cameraPos.Z+10*cameraLook.Z,
                  cameraUp.X,cameraUp.Y,cameraUp.Z);
	}

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);

	/****************************
	/ Add your objects from here
	****************************/
	glPushMatrix();{
	//add objects

	drawAxes();
	drawGrid();

//    TASK=CUBE_SPHERE;
//    switch(TASK)    {
//        case CUBE_SPHERE:
////            drawCubeSphere(cubeLen,sphereRad);
//            break;
//    }

//    drawCheckerBoard(height,width);

    for(int i = 0; i < spheres.size(); i++) {
        spheres[i].draw();
    }
    for(int i = 0; i < pyramids.size(); i++) {
        pyramids[i].draw();
    }
    for(int i = 0; i < lights.size(); i++) {
        Sphere sp(lights[i], 1,Color(0,1,0));
//        cout<<"light drawing: ";
//        sp.center.print();
        sp.draw();
    }
    checkBoard.draw();



//    Vector pos=Point(-100,-100,50);
//    Vector Look=Point(1/sqrt(2.0),1/sqrt(2.0),0);
//    Vector Up=Point(0,0,1);
//    Vector Right=Point(1/sqrt(2.0),-1/sqrt(2.0),0);
//
//    glPushMatrix();{
//        glTranslatef(pos.X,pos.Y,pos.Z);
//        glColor3f(.9,.3,.6);
//        drawSphere(4,20,20);
//    }glPopMatrix();
//
//    int height=100,width=100;
//    for(int row = 0; row<height; row++)   {
//        for(int colm = 0; colm<width; colm++) {
//            double ht = row-height/2;
//            ht/=height/2;
//            double wd = colm-width/2;
//            wd/=width/2;
////            cout<<ht<<" "<<wd<<endl;
//            Vector dir = (Up*ht) + (Right*wd);
//
//            dir = Look+dir;
//            dir = dir*300;
//            dir = dir+pos;
//            glColor3f(.5,.6,.7);
//            glBegin(GL_LINES);   {
//                glVertex3f(pos.X,pos.Y,pos.Z);
//                glVertex3f(dir.X,dir.Y,dir.Z);
//            }glEnd();
//        }
//    }

    }glPopMatrix();
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate()
{
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init()
{
	//codes for initialization
	TASK=CUBE_SPHERE;

	/** Camera initialization **/
    cameraPos=Point(100,100,50);
    cameraLook=Point(-1/sqrt(2.0),-1/sqrt(2.0),0);
    cameraUp=Point(0,0,1);
    cameraRight=Point(-1/sqrt(2.0),1/sqrt(2.0),0);
	drawgrid=0;
	drawaxes=1;
	cameraHeight=80;
	cameraAngle=pi/4;

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
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void readInput()
{
    cin>>recursionLevel;
//    cout<<recursionLevel<<endl;
    cin>>screenWidth;
    screenHeight=screenWidth;
    checkBoard = CheckBoard(screenHeight,screenWidth);
    int totalObj;
    cin>>totalObj;

    string type;
//    cout<<totalObj<<endl;
    for(int i = 0; i < totalObj;i++)    {
        cin>>type;
//        cout<<i<<type<<endl;
        if(type=="sphere")  {
            Sphere sp;
            sp.read();
            spheres.push_back(sp);
        }
        else if(type == "pyramid")  {
            Point lp;
            double w, h;
            Color c;
            double am,df,spc, refl,shn;
            cin>>lp.X>>lp.Y>>lp.Z;
            cin>>w>>h;
            cin>>c.r>>c.g>>c.b;
            cin>>am>>df>>spc>>refl;
            cin>>shn;
//            cout<<"py: "<<shn<<endl;
            Pyramid py(lp,w,h,c);
            py.amb = am;
            py.diffuse = df;
            py.spec = spc;
            py.refl = refl;
            py.shine = shn;
            pyramids.push_back(py);
        }
    }
    int lightCnt;
    cin>>lightCnt;
    lights.clear();
    for(int i = 0; i < lightCnt; i++)   {
        Point lg;
        cin>>lg.X>>lg.Y>>lg.Z;
        lights.push_back(lg);
//        cout<<"light: "<<endl;
//        lg.print();
    }
//    cout<<"END OF INPUT"<<endl;
}

void printImage()
{
    bitmap_image image(screenWidth,screenHeight);  // Creating an image
    for(int i=0;i<screenWidth;i++)  {
        for(int j=0;j<screenHeight;j++) {
            Color &c = imageMap[i][j];
            image.set_pixel(j,screenWidth-i-1,c.r*255,c.g*255,c.b*255);  // Setting the color of the pixels
        }
    }
    image.save_image("out.bmp");    // Saving the image in a file
}

int main(int argc, char **argv)
{
//    ifstream fin("description.txt");

    freopen("description.txt","r",stdin);
    freopen("testing.txt","w",stdout);
    readInput();

	glutInit(&argc,argv);
	glutInitWindowSize(screenHeight,screenWidth);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray-Tracer");

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

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
    a/=2;
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}

void drawCircle(double radius,int segments)
{
    int i;
    Point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].X=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].Y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].X,points[i].Y,0);
			glVertex3f(points[i+1].X,points[i+1].Y,0);
        }
        glEnd();
    }
}

void drawSphereQuarter(double radius,int slices,int stacks,int quarter) //quarter is 0-based
{
    glPushMatrix();{
//    quarter--;
    glRotatef(90*(quarter%4),0,0,1);
    int rev=1;
    if(quarter>=4)  rev=-1;
	Point points[stacks+2][slices+2];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
	    double ang=((double)i/(double)stacks)*(pi/2);
		h=radius*sin(ang);
		r=radius*cos(ang);
		for(j=0;j<=slices;j++)
		{
		    double ang2=((double)j/(double)slices)*pi/2;
			points[i][j].X=r*cos(ang2);
			points[i][j].Y=r*sin(ang2);
			points[i][j].Z=h*rev;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
//        double shade;
//        if(i<stacks/2)  shade=((double)i/(double)stacks);
//        else  shade=(1-(double)i/(double)stacks);
////        shade=1-shade;
//        glColor3f(shade,shade,shade);

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].X,points[i][j].Y,points[i][j].Z);
				glVertex3f(points[i][j+1].X,points[i][j+1].Y,points[i][j+1].Z);
				glVertex3f(points[i+1][j+1].X,points[i+1][j+1].Y,points[i+1][j+1].Z);
				glVertex3f(points[i+1][j].X,points[i+1][j].Y,points[i+1][j].Z);
//                //lower hemisphere
//                glVertex3f(points[i][j].X,points[i][j].Y,-points[i][j].Z);
//				glVertex3f(points[i][j+1].X,points[i][j+1].Y,-points[i][j+1].Z);
//				glVertex3f(points[i+1][j+1].X,points[i+1][j+1].Y,-points[i+1][j+1].Z);
//				glVertex3f(points[i+1][j].X,points[i+1][j].Y,-points[i+1][j].Z);
			}glEnd();
		}
	}
	}glPopMatrix();
}

void drawSphere(double radius,int slices,int stacks)
{
    for(int i=0;i<8;i++)    drawSphereQuarter(radius,slices,stacks,i);
}

void drawQuarterCylinder(double height,double radius,int segments,int quarter)   //half cylinder is on positive z-axis, another half is on the negative.
{
    glPushMatrix();{
    glRotatef(90*quarter,0,0,1);
    Point points[segments+2];
    int i;
    height/=2;
    for(i=0;i<=segments;i++) {
        double ang=((double)i/(double)segments)*pi/2;
        points[i].X=radius*cos(ang);
        points[i].Y=radius*sin(ang);
        points[i].Z=height;
    }
    double shade;
    for(i=0;i<segments;i++) {
//        double shade;
//        if(i<segments/2)  shade=((double)i/(double)segments);
//        else  shade=(1-(double)i/(double)segments);
////        shade=1-shade;
//        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);{
            glVertex3f(0,0,height);
            glVertex3f(points[i].X,points[i].Y,points[i].Z);
            glVertex3f(points[i+1].X,points[i+1].Y,points[i+1].Z);
        }glEnd();

        glBegin(GL_QUADS);{
            glVertex3f(points[i].X,points[i].Y,points[i].Z);
            glVertex3f(points[i+1].X,points[i+1].Y,points[i+1].Z);
            glVertex3f(points[i+1].X,points[i+1].Y,-points[i+1].Z);
            glVertex3f(points[i].X,points[i].Y,-points[i].Z);
        }glEnd();
    }

    }glPopMatrix();
}

void drawCubeSphere(double a,double r)
{
    int i;
    int segments=20;
    glPushMatrix(); {
        a/=2;
        a-=r;
//        int i=0;
        /****DRAWING CORNER QUARTER-SPHERES****/
        for(i=0;i<8;i++)    {
//            if(i%2==0) continue;
            glPushMatrix();{
                glTranslatef(a*dirX[i],a*dirY[i],a*dirZ[i]);
                drawSphereQuarter(r,segments,segments,i);
            }glPopMatrix();
        }

        /***Drawing corner quarter cylinders***/


        for(i=0;i<=2;i++)    {
            if(i==1)   glRotatef(90,1,0,0);
            if(i==2)   glRotatef(90,0,1,0);
            for(int j=0;j<4;j++)    {
                glPushMatrix();{
                    glTranslatef(a*dirX[j],a*dirY[j],0);
                    drawQuarterCylinder(a*2,r,segments,j);
                }glPopMatrix();
            }
            if(i==2)   glRotatef(-90,0,1,0);
            if(i)   glRotatef(-90,1,0,0);
        }
        a+=r;

        /***Drawing side walls***/
        glPushMatrix(); {
            glColor3f(1,1,1);
//            drawSquare(a-r);
            double sqLen=2*(a-r);
            for(i=0;i<4;i++)    {
                glPushMatrix();{
                    glRotatef(90*i,0,0,1);
                    glTranslatef(a,0,0);
                    glRotatef(90,0,1,0);

                    drawSquare(sqLen);
                }glPopMatrix();
            }
            glTranslatef(0,0,a);
            drawSquare(sqLen);
            glTranslatef(0,0,-2*a);
            drawSquare(sqLen);
//            cout<<"a-r: "<<a-r<<endl;
        }glPopMatrix();

    }
    glPopMatrix();
}

/***Keyboard and mouse listeners***/
void keyboardListener(unsigned char key, int x,int y)
{
//    cout<<key<<" "<<x<<" "<<y<<endl;
//    if(TASK==CAMERA_ANGLE)  {
    Vector l=cameraLook;
    Vector r=cameraRight;
    Vector u=cameraUp;
    double ang=.05;

    switch(key) {
        case '1':   //rotate right. rotate r and l WRT u;
            ang*=-1;
            cameraRight.X=r.X*cos(ang)+l.X*sin(ang);
            cameraRight.Y=r.Y*cos(ang)+l.Y*sin(ang);
            cameraRight.Z=r.Z*cos(ang)+l.Z*sin(ang);

            cameraLook=crossProduct(cameraUp,cameraRight);
//                cameraUp.print();
            break;
        case '2':
            cameraRight.X=r.X*cos(ang)+l.X*sin(ang);
            cameraRight.Y=r.Y*cos(ang)+l.Y*sin(ang);
            cameraRight.Z=r.Z*cos(ang)+l.Z*sin(ang);

            cameraLook=crossProduct(cameraUp,cameraRight);
            break;
        case '3':       //rotate up. rotate l and u vectors WRT r.
            cameraLook.X=l.X*cos(ang)+u.X*sin(ang);
            cameraLook.Y=l.Y*cos(ang)+u.Y*sin(ang);
            cameraLook.Z=l.Z*cos(ang)+u.Z*sin(ang);

            cameraUp=crossProduct(cameraRight,cameraLook);
            break;
        case '4':       //rotate down. rotate l and u vectors WRT r.
            ang*=-1;
            cameraLook.X=l.X*cos(ang)+u.X*sin(ang);
            cameraLook.Y=l.Y*cos(ang)+u.Y*sin(ang);
            cameraLook.Z=l.Z*cos(ang)+u.Z*sin(ang);

            cameraUp=crossProduct(cameraRight,cameraLook);
            break;

        case '5':       //tilt camera clockwise.rotate r and u WRT l.
            cameraUp.X=u.X*cos(ang)+r.X*sin(ang);
            cameraUp.Y=u.Y*cos(ang)+r.Y*sin(ang);
            cameraUp.Z=u.Z*cos(ang)+r.Z*sin(ang);

            cameraRight=crossProduct(cameraLook,cameraUp);

            break;
        case '6':
            ang*=-1;
            cameraUp.X=u.X*cos(ang)+r.X*sin(ang);
            cameraUp.Y=u.Y*cos(ang)+r.Y*sin(ang);
            cameraUp.Z=u.Z*cos(ang)+r.Z*sin(ang);

            cameraRight=crossProduct(cameraLook,cameraUp);
            break;
        case '0':
            cout<<"PRINTING IMAGE: "<<endl;
            rayTracer();
            printImage();
            break;
    }
//    cout<<"KKEY: "<<key<<endl;
//        cout<<"Up: ";cameraUp.print();
//        cout<<"Right: ";cameraRight.print();
//        cout<<"Look: ";cameraLook.print();
//        cout<<"CamPos: ";cameraPos.print();
//    }

}

void specialKeyListener(int key, int x,int y)
{
//    if(TASK==CUBE_SPHERE)   {
    double mov=3;
    switch(key) {
        case GLUT_KEY_HOME:
            if(sphereRad<cubeLen/2) sphereRad++;
            break;
        case GLUT_KEY_END:
            if(sphereRad)   sphereRad--;
            break;
        case GLUT_KEY_PAGE_DOWN:		//down arrow key
            cameraPos.X -= mov*cameraUp.X;
            cameraPos.Y -= mov*cameraUp.Y;
            cameraPos.Z -= mov*cameraUp.Z;
            break;
        case GLUT_KEY_PAGE_UP:		// up arrow key
            cameraPos.X += mov*cameraUp.X;
            cameraPos.Y += mov*cameraUp.Y;
            cameraPos.Z += mov*cameraUp.Z;
            break;

        case GLUT_KEY_RIGHT:
            cout<<"RIGHT PRESSED"<<endl;
            cameraPos.X+=mov*cameraRight.X;
            cameraPos.Y+=mov*cameraRight.Y;
            cameraPos.Z+=mov*cameraRight.Z;
            break;
        case GLUT_KEY_LEFT:
            cameraPos.X-=mov*cameraRight.X;
            cameraPos.Y-=mov*cameraRight.Y;
            cameraPos.Z-=mov*cameraRight.Z;
            break;

        case GLUT_KEY_UP:
            cameraPos.X+=mov*cameraLook.X;
            cameraPos.Y+=mov*cameraLook.Y;
            cameraPos.Z+=mov*cameraLook.Z;
            break;
        case GLUT_KEY_DOWN:
            cameraPos.X-=mov*cameraLook.X;
            cameraPos.Y-=mov*cameraLook.Y;
            cameraPos.Z-=mov*cameraLook.Z;
            break;
    }
//    cout<<"SP KEY:  "<<key<<endl;
//        cout<<"CamPos: ";
//        cameraPos.print();
//    }
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
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
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
