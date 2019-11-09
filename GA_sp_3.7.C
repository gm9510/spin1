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

  if( argc == 2 ) {
	int s=stoi( argv[1] );
  stringstream Ds;
  Ds << fixed << setprecision(0) << s;

 string dir  ("Data_3_7/");
  string name ("GA_"), ext(".dat");

fstream fout;

const int Dsize=201, thsize=21;
const int n=3 ;// s=20 ;// n -> internal satates , s -> sites in the lattice
const double dt=0.0001;// time - step
const double tf=100000000;//time iteration
//double D, th; //Interaction Parameters.
const double DD=0.01, Dth=0.01; //Step variations of the interaction parameters.
//double Ma[Dsize*thsize], Tau[Dsize*thsize], Xi[Dsize*thsize];//Observables: Ma -> Magnetization ; Tau -> Quirality ; Xi -> Suceptibility.
double En0=0,En1=0;
int countzero = 0; 

double Parameterth[Dsize*thsize];
double ParameterD[Dsize*thsize];

//----------------------------------------------------
//----------------------------------------------------
for(int Di=0; Di<Dsize; Di++){
for(int thi=0; thi<thsize; thi++){
Parameterth[thi+thsize*Di]= (thi-90)*Dth*PI;
ParameterD[thi+thsize*Di]= (Di-100)*DD;
 }
}
//----------------------------------------------------
//----------------------------------------------------


#pragma omp parallel for private(fout,countzero,En0,En1) schedule(dynamic)
for(int i=0; i<Dsize*thsize; i++){
#pragma omp critical (cout)
{
cout<<"("<<Parameterth[i]<<","<<ParameterD[i]<<")"<<endl;
}
ComplexMatrix f(s,n); //Gutzwiller's Coefficient
ComplexMatrix prevf(s,n);
countzero=0;

stringstream thstr;
  thstr << fixed << setprecision(3) <<Parameterth[i];
stringstream Dstr;
  Dstr << fixed << setprecision(3) <<ParameterD[i];

fout.open(dir+name+"D"+Dstr.str()+"_th"+thstr.str()+"_s"+Ds.str()+ext, ios::out);

for(int si=0; si<s; si++){
  f(si,0) = 0.5773;
  f(si,1) = 0.5773;
  f(si,2) = 0.5773;
}
//----------------------------------------------------
//----------------------------------------------------
//---------------->RK4<-------------------------------
for(int t=0; t<tf; t++){
 RK4(dt,f,n,s,ParameterD[i],0.,Parameterth[i]);
//----------------------------------------------------
//----------------------------------------------------
//---------------------norm---------------------------
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
  En1 = En0;
  En0 = E( f, n, s, Parameterth[i], ParameterD[i] , 0.);

if( abs( En1 - En0 )<zero )
  countzero += 1;
else
  countzero = 0;

if( countzero==100*s*s )
  t=tf;

 }//-------------------End_Of_RK4--------------------------------

for(int si=0; si<s; si++){
fout<<si<<"	"<<f(si,0)<<"	"<<f(si,1)<<"	"<<f(si,2)<<endl;
}

fout.close();
}//--------------------End_Of_Parallel---------------------------


   }
   else if( argc > 2 ) {
      printf("Too many arguments supplied.\n");
   }
   else {
      printf("One argument expected.\n");
   }

return 0;
}
