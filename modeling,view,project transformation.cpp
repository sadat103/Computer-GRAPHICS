#include <bits/stdc++.h>

using namespace std;

#define CLR(a) memset ( a, 0, sizeof ( a ) )

const double PI = acos(-1.0);
const double BN = 1e-12;


inline double degreeToRadian(double ang)
{
    return ang*PI/180.0;
}

inline double radianToDegree(double ang)
{
    return ang*180/PI;
}

struct Vector
{
    double vec[4];
    int dimnsn;
    Vector(){CLR(vec);dimnsn=3;}
    Vector(double a,double b,double c)   
    {
        vec[0]=a;
        vec[1]=b;
        vec[2]=c;
        vec[3]=1;
        dimnsn=3;
    }
    Vector(Vector &vv)   
    {
        for(int i=0;i<=dimnsn;i++)   
	{ 
 		vec[i]=vv.vec[i];
	}
    }
    bool operator == (const Vector &vvc) const  
    {
        for(int i=0;i<=dimnsn;i++)   
        {
            if(abs(vec[i]-vvc.vec[i])>BN)   return false;
        }
        return true;
    }
    Vector operator *(const double &x)   const 
    {
        Vector a;
        for(int i=0;i<=dimnsn;i++)   
 	{
 	        a.vec[i]=vec[i]*x;
        }
        a.vec[dimnsn]=1;
        return a;
    }
    Vector operator *(const Vector& x)   const 
    { 
        Vector a;
        a.vec[0]=vec[1]*x.vec[2]-vec[2]*x.vec[1];
        a.vec[1]=-(vec[0]*x.vec[2]-vec[2]*x.vec[0]);
        a.vec[2]=vec[0]*x.vec[1]-vec[1]*x.vec[0];
        a.vec[3]=1;

        return a;
    }
    Vector operator +(const Vector &x)  const 
    {
        Vector a;
        for(int i=0;i<dimnsn;i++)    
        {
            a.vec[i]=vec[i]+x.vec[i];

        }
        a.vec[dimnsn]=1;
        return a;
    }
    double absoulateVal() 
    {
        if(vec[dimnsn]>1+BN)  
        {
            for(int i=0;i<dimnsn;i++)    
            {
                 vec[i]/=vec[dimnsn];
            }
            vec[dimnsn]=1;
        }
        return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    }
    void normlz()   
    {
        double val=this->absoulateVal();
        for(int i=0;i<dimnsn;i++)    vec[i]/=val;
    }
    void print(ostream &ff = cout,int prcsn=7)    
    {

        ff<<fixed<<setprecision(prcsn)<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<endl;

    }
};

double dotPrdct(Vector a,Vector b)
{
    double p=0;
    for(int i=0;i<a.dimnsn;i++)  p+=a.vec[i]*b.vec[i];
    return p;
}

struct Matrix
{
    int dimnsn;
    double array[4][4];
    Matrix()    {dimnsn=3;CLR(array);}
    Matrix(double a[4][4])  
    {
        for(int i =0;i<4;i++)   
          {
              for(int j=0;j<4;j++)    
                {  
                       array[i][j]=a[i][j];
                }
	  }
    }
    Matrix(const Matrix &mtx) 
    {
        dimnsn=mtx.dimnsn;
        for(int i=0;i<=dimnsn;i++)   
        {
            for(int j=0;j<=dimnsn;j++)
               {
                     array[i][j]=mtx.array[i][j];
               }
             
        }
    }
    void scale(double p)    
    {
        for(int i=0;i<=dimnsn;i++)   
        {
            for(int j=0;j<=dimnsn;j++)   
            {
                array[i][j]*=p;
            }
        }
    }
    Matrix operator * (const Matrix &mtx) const   
    {
        Matrix a;
        for(int i=0;i<=dimnsn;i++)   
        {
            for(int j=0;j<=dimnsn;j++)   
        {
                double tmp=0;
                for(int k=0;k<=dimnsn;k++)   
                {
                    tmp+=array[i][k]*mtx.array[k][j];
                }
                a.array[i][j]=tmp;
            }
        }
        return a;
    }
    Matrix operator + (const Matrix &m) const   
    {
        Matrix a;
        for(int i=0;i<=dimnsn;i++)   
        {
            for(int j=0;j<=dimnsn;j++)   
            {
                a.array[i][j]=m.array[i][j]+array[i][j];
            }
        }
        return a;
    }
    
    Vector operator * (const Vector &vct) const   
    {
        Vector a;
        for(int i=0;i<=dimnsn;i++)   
        {
            for(int j=0;j<1;j++)   
            {
                double tmp=0;
                for(int k=0;k<=dimnsn;k++)   
                {
                    tmp+=array[i][k]*vct.vec[k];
                }
                a.vec[i]=tmp;
            }
        }
        return a;
    }
  
};

Vector eye,cmLk,cmRt,cmUp;
double fovY,aspect,near,far;
stack<Matrix>mtStck;
stack<int>mtIndStck;
Matrix scaleMat;
Matrix idn;

void stage1()
{
    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");


    fin>>eye.vec[0]>>eye.vec[1]>>eye.vec[2];
    fin>>cmLk.vec[0]>>cmLk.vec[1]>>cmLk.vec[2];
    fin>>cmUp.vec[0]>>cmUp.vec[1]>>cmUp.vec[2];
    fin>>fovY>>aspect>>near>>far;

    /** stack with an idntity matrix**/
    while(mtStck.empty()==false)    
    {
        mtStck.pop();
    }
    while(!mtIndStck.empty())    mtIndStck.pop();
    string str;
    scaleMat=idn;
    mtStck.push(idn);
    mtIndStck.push(mtStck.size());

    while(fin>>str) {
        if(str=="triangle")	{
            Vector tring[3];
            for(int i=0;i<3;i++)    {
                fin>>tring[i].vec[0]>>tring[i].vec[1]>>tring[i].vec[2];
                tring[i].vec[3]=1;
                Matrix tmp=mtStck.top();
                Vector vt;
                vt = (tmp*tring[i]);
                vt=vt*(1.0/vt.vec[3]);
                vt.print(fout);
            }
            fout<<endl;
        }
        if(str=="translate")	{
            Matrix translt = idn;
            for(int i=0;i<idn.dimnsn;i++)   {
                fin>>translt.array[i][translt.dimnsn];
            }
            Matrix tmp=mtStck.top();
            translt = tmp*translt;
            int d=translt.dimnsn;
            translt.scale(1.0/translt.array[d][d]);
            mtStck.push(translt);
        }
        if(str=="scale")	{
            Matrix scl = idn;
            for(int i=0;i<scl.dimnsn;i++)   {
                fin>>scl.array[i][i];
            }
            Matrix tmp=mtStck.top();
            scl = tmp*scl;
            int d = scl.dimnsn;
            scl.scale(1.0/scl.array[d][d]);
            mtStck.push(scl);
        }
        if(str=="rotate")	{
            double ang;
            Vector a;
            Vector x(1,0,0),y(0,1,0),z(0,0,1);
            fin>>ang>>a.vec[0]>>a.vec[1]>>a.vec[2];
            a.vec[3]=1;
            a.normlz();
            ang=degreeToRadian(ang);
            //rodrigrez formula
            x = ((x*cos(ang)+((a*x)*sin(ang)))+((a*dotPrdct(a,x))*(1-cos(ang))));
            y = ((y*cos(ang)+((a*y)*sin(ang)))+((a*dotPrdct(a,y))*(1-cos(ang))));
            z = ((z*cos(ang)+((a*z)*sin(ang)))+((a*dotPrdct(a,z))*(1-cos(ang))));

            Matrix R;

            for(int i=0;i<3;i++)   
            {
                R.array[i][0]=x.vec[i];
                R.array[i][1]=y.vec[i];
                R.array[i][2]=z.vec[i];
            }
            R.array[3][3]=1;
            Matrix tmp=mtStck.top();
            tmp=tmp*R;
            tmp.scale(1.0/tmp.array[3][3]);
            mtStck.push(tmp);
        }
        if(str=="push")	 
        {
            mtIndStck.push(mtStck.size());
        }
        if(str=="pop")	
        {
            if(mtStck.size()==1)  continue;
            int e=mtIndStck.top();
            mtIndStck.pop();
            while(mtStck.size()>e)  
               {
                    mtStck.pop();
               }
        }
        if(str=="end")	
        {
            break;
        }
    }
    fin.close();
    fout.close();
}

void stage2()
{
    ifstream fin("stage1.txt");
    ofstream fout("stage2.txt");
    Vector e1,L,R,U;
    e1=eye*-1.0;
    L = cmLk+e1;
    L.normlz();
    R = L*cmUp;
    R.normlz();
    U = R*L;
    U.normlz();
    Matrix R1;
    Matrix T = idn;
    for(int i=0;i<3;i++)    {
        R1.array[0][i]=R.vec[i];
        R1.array[1][i]=U.vec[i];
        R1.array[2][i]=-L.vec[i];
        T.array[i][3]=-eye.vec[i];
    }
    R1.array[3][3]=1;
    Matrix V = R1*T;
    for(int cntv=1;;cntv++) {
        Vector p;
        if(!(fin>>p.vec[0]>>p.vec[1]>>p.vec[2]))  break;
        p.vec[3]=1;
        p = V*p;
        p.print(fout);
        if(cntv%3==0) fout<<endl;
    }
    fin.close();
    fout.close();
}

void stage3()
{
    ifstream fin("stage2.txt");
    ofstream fout("stage3.txt");
    Matrix P1;
    double fovX,t,r;
    fovX = fovY * aspect;
    t=near * tan(degreeToRadian(fovY/2.0));
    r=near * tan(degreeToRadian(fovX/2.0));
    CLR(P1.array);

    P1.array[0][0]=near/r;
    P1.array[1][1]=near/t;
    P1.array[2][2]=-(far+near)/(far-near);
    P1.array[2][3]=-(2*far*near)/(far-near);
    P1.array[3][2]=-1;
    int i,j;

    for(int cntv=1;;cntv++) 
   {
        Vector p;
        if(!(fin>>p.vec[0]>>p.vec[1]>>p.vec[2]))  break;
        p.vec[3]=1;
        p = P1*p;
        p = p*(1/p.vec[3]);
        p.print(fout,7);
        if(cntv%3==0) fout<<endl;
    }
    fin.close();
    fout.close();
}

int main ()
{
    int i, j, k;
    idn.dimnsn=3;
    for(int i=0;i<=idn.dimnsn;i++)  
    {
        idn.array[i][i]=1;
    }

    stage1();
    stage2();
    stage3();
    return 0;
}
