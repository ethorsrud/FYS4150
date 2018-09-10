#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <armadillo>

using namespace arma;
using namespace std;
ofstream ofile;

int main(int argc, char *argv[]){

  int expn = atoi(argv[1]); //Reading exponent from commandline
  string file = "power-";
  file.append(argv[1]); file.append(".txt"); //Making filename power-x.txt where x = exponent

  int n = pow(10.0,expn); //10^expn
  double h = 1.0/n; //Stepsize

  //Array defining and memory allocation
  double *d = new double [n+1];
  double *b = new double [n+1];
  double *x = new double [n+1];
  double *u = new double [n+1];
  double *u_exact = new double [n+1];

  //initials
  u[0] = u[n] = 0.0;

  clock_t clock_start, clock_finish;

  for(int i = 0;i<=n;i++){
    d[i] = 2.0;
  }

  for(int i = 0;i<=n;i++) {
    x[i] = i*h;
    b[i] = h*h*100.0*exp(-10.0*x[i]); //defining b = h^2*f
    u_exact[i] = 1-(1-exp(-10))*x[i]-exp(-10*x[i]); //exact solution to u
  }

  clock_start = clock();

//Forward subst.
  for(int i = 2;i<n;i++) {
    d[i] = (i+1.0)/((double) i);
    b[i] = b[i] + (b[i-1]/d[i-1]);
  }

// Backward subst.
  u[n-1] = b[n-1]/d[n-1];
  for(int i = n-2;i>0;i--){
    u[i] = (b[i]+u[i+1])/d[i];
  }

  clock_finish = clock();
  double time_used = (clock_finish-clock_start)/(float)CLOCKS_PER_SEC;
  cout << "Tridiagonal special case calculations took: "<< time_used << " Seconds" << endl;
  cout << "Now writing to file " << file << " ..." << endl;
  
  //Write to file
  ofile.open(file);
  for(int i = 1; i<n;i++){
    // x - Tridiagonal solution - Exact solution - Relative error
    double rel_err_log = log10(abs((u[i]-u_exact[i])/u_exact[i]));
    ofile << setw(15) << setprecision(8) << x[i];
    ofile << setw(15) << setprecision(8) << u[i];
    ofile << setw(15) << setprecision(8) << u_exact[i];
    ofile << setw(15) << setprecision(8) << rel_err_log << endl;
  }

//Closing file and freeing memory
  ofile.close();
  delete [] x;
  delete [] b;
  delete [] d;
  delete [] u;
  delete [] u_exact;

}