#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <armadillo>
#include "transaction.h"
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

using namespace arma;

using namespace std;


int main()
{
    ofstream ofile;
    //change file name according to values selected for lambda, gamma, alpha
    ofile.open("N1000_lamb05_alph20_gam40.txt");
    //Number of agents
    int n_agents = 1000;
    //Starting value
    int m0 = 100;
    //Number of transactions
    int n_transactions = 1e7;

    double lambda = 0.5;
    double gamma = 4.0;
    double alpha = 2.0;
    //changing seed with time
    srand(time(0));
    //Number of simulations
    int number_of_sim = 100;
    for(int i = 0; i<number_of_sim;i++ ){
        transaction model = transaction(n_agents,n_transactions,m0,lambda,alpha,gamma);
        //Writing to file
        for(int j=0;j<n_agents;j++){
           ofile << setw(15) << setprecision(8) << model.wealth(j)<<endl;
        }

        cout<<"Sim: "<<i+1<<" of "<<number_of_sim<<endl;

    }
}
