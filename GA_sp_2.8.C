#include <lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <matrixtypes.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>//libreria estandar modo c
#include <complex>
#include <functional>
#include <numeric>
#include <math.h>
#include <sstream>
#include <omp.h>

#include <functions.h>

using namespace ula;
using namespace std;
using namespace myf; //include the functions E(...); F_IT(...); RK4(...); M; Q; X;and teh complex I;

#define root2 1.41421356237
#define zero 0.00000001
#define PI 3.14159265359

int main( int argc, char *argv[] ){

  if( argc == 3 ) {
	int s=stoi( argv[1] );
	double th = atof( argv[2] );

 string dir  ("Data_2_8/");
  string name ("GA_D_"), ext(".dat");

  stringstream Ds;
  Ds << fixed << setprecision(0) << s;
  stringstream Dthi;
  Dthi << fixed << setprecision(3) << th;

fstream Mout, Qout, Xout, Eout, Sxout, Syout, MMout, SxSxout, SySyout, SpSpout, SpSmout;
	Mout.open(dir+name+"M_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	Qout.open(dir+name+"Q_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	Xout.open(dir+name+"X_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	Eout.open(dir+name+"E_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	Sxout.open(dir+name+"Sx_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	Syout.open(dir+name+"Sy_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	MMout.open(dir+name+"MM_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	SxSxout.open(dir+name+"SxSx_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	SySyout.open(dir+name+"SySy_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	SpSpout.open(dir+name+"SpSp_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);
	SpSmout.open(dir+name+"SpSm_th"+Dthi.str()+"_s"+Ds.str()+ext, ios::out);

int n=3 ; //s=30 ;// n -> internal satates , s -> sites in the lattice
double dt=0.0001;// time - step
double tf=100000000;//time iteration
double D=0; //th=-0.85*PI; //Interaction Parameters.
double DD=0.01; //Step variations of the interaction parameters.
double Ma=0, Tau=0, Xi=0;//Observables: Ma -> Magnetization ; Tau -> Quirality ; Xi -> Suceptibility.
Complex SX=0,SY=0.;
double SxSx=0., SySy=0.,MM=0.,SpSm=0., SpSp=0.;
RealVector En(2);
int countzero = 0; 

ComplexMatrix f(s,n); //Gutzwiller's Coefficient
ComplexMatrix prevf(s,n);

//----------------------------------------------------
//----------------------------------------------------
//Initial Condition
for(int si=0; si<s; si++){
  f(si,0) = 0.632455532;
  f(si,1) = 0.547722558;
  f(si,2) = 0.547722558;
}

//----------------------------------------------------
//----------------------------------------------------


//for(int si=0; si<s; si++)
//cout<<f(si,0)<<"	"<<f(si,1)<<"	"<<f(si,2)<<"	"<<endl;

 for(int Di=0; Di<=int(4/DD); Di++){
  D=(Di-200)*DD;
//----------------------------------------------------
//----------------------------------------------------
//---------------->RK4<-------------------------------
prevf=f;//Saves the previus result.

//Initial Condition again
for(int si=0; si<s; si++){
  f(si,0) = 0.632455532;
  f(si,1) = 0.547722558;
  f(si,2) = 0.547722558;
}
for(int t=0; t<tf; t++){
	RK4(dt,f,n,s,D,0.,th);
//----------------------------------------------------
//----------------------------------------------------

//norm
RealVector norm(s);

for(int si=0; si<s; si++){
	norm(si)=sqrt(real(conj(f(si,0))*f(si,0) + conj(f(si,1))*f(si,1) +conj(f(si,2))*f(si,2)));

// cout<<"norma1="<<norm(si)<<endl;

	
if(norm(si) != 1 ){
	f(si,0)=f(si,0)/norm(si);
	f(si,1)=f(si,1)/norm(si);
	f(si,2)=f(si,2)/norm(si);
	}

	norm(si)=real(conj(f(si,0))*f(si,0) + conj(f(si,1))*f(si,1) +conj(f(si,2))*f(si,2));
 //cout<<"norma2="<<norm(si)<<endl;
}
//----------------------------------------------------
//----------------------------------------------------
En(1) = En(0);
  En(0) = E( f, n, s, th, D ,0.);
if( abs( En(1) - En(0) )<zero )
  countzero += 1;
else
  countzero = 0;

//cout<<countzero<<"	"<<t<<"	"<<(abs( En(1) - En(0) )<zero)<<endl;

if( countzero==100*s*s )
  t=tf;
//----------------------------------------------------
//----------------------------------------------------
 }

RealVector norm(s);

//for(int si=0; si<s; si++){
// norm(si) = real(conj(f(si,0))*f(si,0) + conj(f(si,1))*f(si,1) +conj(f(si,2))*f(si,2));
//cout<<"norma final = "<<norm(si)<<endl;
//}


//for(int si=0; si<s; si++)
//cout<<f(si,0)<<"	"<<f(si,1)<<"	"<<f(si,2)<<"	"<<endl;

//cout<<"Interaction D="<<D<<endl;
Ma=M(f,n,s);
//cout<<"Magnetization: "<<Ma<<endl;
Tau=T(f,n,s);
//cout<<"Quirality: "<<Tau<<endl;
Xi=X(f,prevf,n,s, DD);
//cout<<"Suceptibility: "<<Xi<<endl;
En(0)=E(f, n, s, th, D, 0.);
//cout<<"Energy :"<<En(0)<<endl;
SX=Sx(f,n,s);
SY=Sy(f,n,s);
SxSx=SxbySx(f,n,s);
SySy=SybySy(f,n,s);
MM=MbyM(f,n,s);
SpSp=SpbySp(f,n,s);
SpSm=SpbySm(f,n,s);
Mout<<D<<"	"<<Ma<<endl;
Qout<<D<<"	"<<Tau<<endl;
Xout<<D<<"	"<<Xi<<endl;
Eout<<D<<"	"<<En(0)<<endl;
Sxout<<D<<"	"<<real(conj(SX)*SX)<<"	"<<atan2(imag(SX),real(SX))<<endl;
Syout<<D<<"	"<<real(conj(SY)*SY)<<"	"<<atan2(imag(SY),real(SY))<<endl;
SxSxout<<D<<"	"<<SxSx<<endl;
SySyout<<D<<"	"<<SySy<<endl;
SpSpout<<D<<"	"<<SpSp<<endl;
SpSmout<<D<<"	"<<SpSm<<endl;
MMout<<D<<"	"<<MM<<endl;
}//-----end-of-D-for------------


Mout.close();
Qout.close();
Xout.close();
Eout.close();
Sxout.close();
Syout.close();
SxSxout.close();
SySyout.close();
SpSpout.close();
SpSmout.close();
MMout.close();

   }
   else if( argc > 3 ) {
      printf("Too many arguments supplied.\n");
   }
   else {
      printf("two argument expected.\n");
   }
return 0;
}



