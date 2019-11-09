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
using namespace myf; //include the functions E(...); F_IT(...); RK4(...); and teh complex I;

#define root2 1.41421356237
#define zero 0.00000001
#define PI 3.14159265359


int main( int argc, char *argv[] ){
   if( argc == 3 ) {
	double tf = atof(argv[1]);
	double dt = atof(argv[2]);
	if(tf==1 && dt==1){
		dt=0.0001;// time - step
		tf=100000000;//time iteration
		}
      cout<<"The arguments supplieds are "<<endl;
	cout<<"Iteration number :"<<tf<<endl;
	cout<<"Time step :"<<dt<<endl;
 string dir  ("Data_1_5/Unbalanced/");
  string name ("GA_"), ext(".dat");

const int n=3 , s=50 ;// n -> internal satates , s -> sites in the lattice
const double D1=-1.5, D2=0., th=-0.9*PI;
double Ma=0, Tau=0;//Observables: Ma -> Magnetization ; Tau -> Quirality
Complex SX=0., SY=0.;
RealVector En(2);
int countzero = 0; 

En(0)=0;
En(1)=0;

ComplexMatrix f(s,n); //Gutzwiller's Coefficient
string stf= argv[1], sdt =argv[2];
fstream fout, Densout;

stringstream thstr;
  thstr << fixed << setprecision(3) <<th;
stringstream Dstr;
  Dstr << fixed << setprecision(2) <<D1;
stringstream D2str;
  D2str << fixed << setprecision(1) <<D2;
stringstream Ds;
  Ds << fixed << setprecision(0) << s;

 fout.open(dir+name+"Conv_D"+Dstr.str()+"_th"+thstr.str()+"_s"+Ds.str()+"_D2"+D2str.str()+ext, ios::out);
 Densout.open(dir+name+"Perfil_D"+Dstr.str()+"_th"+thstr.str()+"_s"+Ds.str()+"_D2"+D2str.str()+ext, ios::out);

//----------------------------------------------------
//----------------------------------------------------
//Initial Condition
for(int si=0; si<s; si++){
  f(si,0) = 0.632455532;//0.8944427191;//0.577350269
  f(si,1) = 0.547722558;
  f(si,2) = 0.547722558;//0.316227766;
}

//----------------------------------------------------
//----------------------------------------------------


for(int si=0; si<s; si++)
cout<<f(si,0)<<"	"<<f(si,1)<<"	"<<f(si,2)<<"	"<<endl;


//----------------------------------------------------
//----------------------------------------------------
//---------------->RK4<-------------------------------

for(int t=0; t<tf; t++){
	RK4(dt,f,n,s,D1,D2,th);
//----------------------------------------------------
//----------------------------------------------------


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

  En(1) = En(0);
  En(0) = E( f, n, s, th, D1, D2 );
  fout<<t<<"	"<<En(0)<<endl;
if( abs( En(1) - En(0) )<zero )
  countzero += 1;
else
  countzero = 0;

cout<<countzero<<"	"<<t<<"	"<<(abs( En(1) - En(0) )<zero)<<endl;

if( countzero==200*s*s )
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

cout<<"Interaction D2 = "<<D2<<endl;
cout<<"Interaction D1 = "<<D1<<endl;
cout<<"Interaction th = "<<th<<endl;
Ma=M(f,n,s);
cout<<"Magnetization: "<<Ma<<endl;
Tau=T(f,n,s);
cout<<"Quirality: "<<Tau<<endl;
En(0)=E(f, n, s, th, D1,D2);
cout<<"Energy :"<<En(0)<<endl;
SX=Sx(f,n,s);
SY=Sy(f,n,s);
cout<<"Sx: "<<real(conj(SX)*SX)<<"	"<<atan2(imag(SX),real(SX))<<endl;
cout<<"Sy: "<<real(conj(SY)*SY)<<"	"<<atan2(imag(SY),real(SY))<<endl;

	for(int si=0; si<s; si++)
	Densout<<si<<"	"<<real(f(si,0)*conj(f(si,0)))<<"	"<<real(f(si,1)*conj(f(si,1)))<<"	"<< real(f(si,2)*conj(f(si,2)))<<endl;

      cout<<"The arguments supplieds were "<<endl;
	cout<<"Iteration number :"<<tf<<endl;
	cout<<"Time step :"<<dt<<endl;
	cout<<"Sites :"<<s<<endl;
fout.close();
Densout.close();

   }
   else if( argc > 3 ) {
      printf("Too many arguments supplied.\n");
   }
   else {
      printf("Two argument expected.\n");
   }
return 0;
}
