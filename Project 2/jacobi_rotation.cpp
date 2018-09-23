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

//Function declaration
void maxdiag(mat A,int *r,int *c,int n);
void jacobi(mat& A, mat& D,int k,int l,int n);
void eigenvalues(mat A,int n,double eig[]);
void fillAD(mat& A, mat& D,int n,double d,double e,double rhomax,double rhomin,double h,double wr,int choice);
void unit_test(int r,int c);


int main(int argc, char *argv[]){
  clock_t clock_start, clock_finish;
  double wr; //omega_r
  int choice; //for deciding the different potentials in the diagonal elements
  cout<<" 1 for buckling beam\n 2 for one electron quantum dots\n 3 for two electron quantum dots\n";
  cin >> choice;
  if(choice==3){
    cout<<"omega_r = "<<endl;
    cin >> wr;
  }
  else{
    wr=0.0;
  }
  int n = atoi(argv[1]); //Iteration points
  double rhomax = 1; double rhomin = 0;
  if(choice > 1){
    rhomax = 5; rhomin = 0;
  }
  double h = (rhomax-rhomin)/n; //Stepsize
  int r; int c; //Row and column index for max off diagonal elements
  double d = 2.0/(h*h); //diagonal element
  double e = -1.0/(h*h); //element above and below diagonal
  double *eig = new double [n]; //eigenvalues
  mat D(n,n); //Diagonal matrix
  mat A(n,n); //Matrix A
  unit_test(r,c); //testing max_off_diag and orthogonality
  //Filling matrixes
  fillAD(A,D,n,d,e,rhomax,rhomin,h,wr,choice); //filling A and D
  //Jacobi rotations
  double epsilon = 1.0e-10;
  double max_off_diag;
  int max_iterations = pow(10.0,6);
  int iterations = 0;
  maxdiag(A,&r,&c,n);
  max_off_diag = fabs(A(r,c));
  clock_start = clock();
  while(max_off_diag>epsilon && iterations<max_iterations){
    jacobi(A,D,r,c,n);
    maxdiag(A,&r,&c,n);
    max_off_diag = fabs(A(r,c));
    iterations+=1;
    cout<<"Calculating: "<<iterations<<" iterations"<<"\r";
    cout.flush();
  }
  cout<<"\n"<<endl;
  clock_finish = clock();
  double time_used = (clock_finish-clock_start)/(float)CLOCKS_PER_SEC;
  cout << "Diagonalizing matrix took: "<< time_used << " Seconds" << endl;
  cout << "\n";
  eigenvalues(A,n,eig);
  cout << "lambda1 = "<<eig[0]<<" "<<"lambda2 = "<<eig[1]<<" "<<"lambda3 = "<<eig[2]<<endl;


}

void maxdiag(mat A,int *r,int *c,int n)
{
  double largest=0.0;
  for(int i = 0; i<n; i++)
  {
    for(int j = i+1;j<n;j++)
    {
      if(fabs(A(i,j))>=largest)
      {
        largest = fabs(A(i,j));
        *r = i;
        *c = j;
      }
    }
  }
}

void jacobi(mat& A, mat& D,int k,int l,int n)
{
  double s,c;
  double t,tau;
  if(A(k,l) == 0){
    c = 1.0;
    s = 0.0;
  }
  else{
    tau = (A(l,l)-A(k,k))/(2.0*A(k,l));
    if (tau>=0){
      t = 1.0/(tau+sqrt(1.0+tau*tau));
    }
    else{
      t = -1.0/(-tau+sqrt(1.0+tau*tau));
    }
    c = 1.0/sqrt(1+t*t);
    s = c*t;

  }
  double a_kk = A(k,k);
  double a_ll = A(l,l);
  A(k,k) = c*c*a_kk-2.0*c*s*A(k,l)+s*s*a_ll;
  A(l,l) = s*s*a_kk+2.0*c*s*A(k,l)+c*c*a_ll;
  A(k,l) = 0.0;
  A(l,k) = 0.0;
  double a_ik,a_il,d_ik,d_il;
  for (int i = 0; i < n; i++) {
  if (i != k && i != l){
    a_ik = A(i,k);
    a_il = A(i,l);
    A(i,k) = c*a_ik - s*a_il;
    A(k,i) = A(i,k);
    A(i,l) = c*a_il + s*a_ik;
    A(l,i) = A(i,l);
  }
  d_ik = D(i,k);
  d_il = D(i,l);
  D(i,k) = c*d_ik-s*d_il;
  D(i,l) = c*d_il+s*d_ik;
  }
}

//For setting an array with sorted eigenvalues
void eigenvalues(mat A,int n,double eig[]){
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(i==j){
        eig[i] = A(i,j);
      }
    }
  }
  sort(eig,eig+n);
}

void fillAD(mat& A, mat& D,int n,double d,double e,double rhomax,double rhomin,double h,double wr,int choice){
  double V;
  for (int i = 0;i<=n-1;i++){
    if(wr==0.0){
    V = (rhomin+(i+1)*h)*(rhomin+(i+1)*h);
    }
    else{
      V = wr*wr*(rhomin+(i+1)*h)*(rhomin+(i+1)*h)+1./(rhomin+(i+1)*h);
    }
    for (int j = 0;j<=n-1;j++){
      if(i==j && choice == 1){
        A(i,j) = d;
        D(i,j) = 1;
      }
      else if(i==j && choice > 1){
        A(i,j) = d+V;
        D(i,j) = 1;
      }
      if(abs(i-j)==1){
        A(i,j) = e;
      }
    }
  }
}

void unit_test(int r,int c){
  mat test = "1,2,3,4,5;2,3,4,10,6;3,4,5,6,7;4,5,6,7,8;5,6,7,8,9";
  mat test2 = "2,-2;1,1";
  mat test2_eig(2,2);
  maxdiag(test,&r,&c,5);
  if(test(r,c)==10){
    cout << "\033[1;32mUnit test for finding max off-diagonal element has passed\033[0m\n" << endl;
  }
  else{
    cout << "\033[1;31mUnit test for finding max off-diagonal element has failed\033[0m\n" << endl;
  }
  maxdiag(test2,&r,&c,2);
  jacobi(test2,test2_eig,r,c,2);
  if(test2.col(0)(0)*test2.col(1)(0)-test2.col(0)(1)*test2.col(1)(1)==0){
    cout << "\033[1;32mUnit test for checking orthogonality has passed\033[0m\n" << endl;
  }
  else{
    cout << "\033[1;31mUnit test for checking orthogonality has failed\033[0m\n" << endl;
  }
}
