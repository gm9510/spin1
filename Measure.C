#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>//libreria estandar modo c
#include <time.h> //numeros aleatorios
#include <complex>
#include <functional>
#include <numeric>
#include <math.h>
#include <sstream>
#include <omp.h>

#include <lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <matrixtypes.h>

#define PI 3.14159265359

using namespace ula;
using namespace std;

void Block_Reader(fstream& file, RealVector& Ra,RealVector& Rb){
        int j=0;
	char Buffer[128],A[16],B[16];
	const char a[2]="";

	  if(!file.is_open())  //Check if file is not open
	    perror("Error opening file");
  
	  else{
	
	  while(!file.eof()){

	    file.getline(Buffer,64); 
	    if(file.fail()||file.eof()){   //If failbit flag is up then break the while loop.
	      perror("Error while reading file");    
	      break;
	           }
	  else{
	   if(strcmp(Buffer,a)==0){ // If we get an empty line then break the while loop.
	    break;
	       }
	
	   else{
	    sscanf(Buffer,"%s	%s",A,B);   //scan the extracted line
	    Ra(j) = atof( A );
	    Rb(j) = atof( B );
	    j++;

	       }
	        }
	         }
		  }
     }


int main(){
fstream finput,fout;
 string dir  ("Data_3_8/");
  string name ("GA_D_"), ext(".dat");
fout.open(dir+"fullX"+ext, ios::out);
double th=0;
double Dth = 0.01;
int a=201;
RealVector X(a),Y(a);

for(int thi=5; thi<=15; thi++){
 th=(thi-90)*Dth*PI;
  stringstream Dthi;
  Dthi << fixed << setprecision(3) << th;
  finput.open(dir+name+"X_th"+Dthi.str()+"_s50"+ext, ios::in);


cout<<dir+name+"X_th"+Dthi.str()+"_s50"+ext<<"	"<<th<<endl;
Block_Reader(finput, X, Y);

for(int i=0; i<a; i++)
 fout<<th<<"	"<<X(i)<<"	"<<Y(i)<<endl;

fout<<endl;

finput.close();
}

fout.close();

return 0;
}
