#include "transaction.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <armadillo>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

transaction::transaction(int N_agents,int n_transactions,int m0,double lambda,double alpha,double gamma)
{
    //Vector of each agents wealth
    wealth = zeros(N_agents)+m0;
    //Matrix counting the number of former transactions between agents
    mat C = zeros(N_agents,N_agents);

    for(int i = 0;i<n_transactions;i++){
        //random comparison value
        double r = (double) rand() / RAND_MAX;
        //randomly choosing two agents
        int random_i = (int) rand()%N_agents;
        int random_j = (int) rand()%N_agents;
        //Check their wealth
        double m_i = wealth(random_i);
        double m_j = wealth(random_j);
        //checking previous transactions
        int previous = C(random_i,random_j);
        //Probability function
        double pij = pow(fabs(m_i-m_j),-alpha)*pow(previous+1,gamma);

        while(pij<r||(random_i==random_j)){
            //pick new agents until they can do a transaction
            random_i = (int) rand()%N_agents;
            random_j = (int) rand()%N_agents;
            m_i = wealth(random_i);
            m_j = wealth(random_j);
            previous = C(random_i,random_j);
            pij = pij = pow(fabs(m_i-m_j),-alpha)*pow(previous+1,gamma);
            r = (double) rand() / RAND_MAX;
         }

        //At this place, the right two agents has been picked

        //random number for determine the new wealths
        double epsilon = (double) rand()/RAND_MAX;
        //new values
        double delta_m = (1-lambda)*(epsilon*m_j-(1-epsilon)*m_i);
        double m_i_after = m_i+delta_m;
        double m_j_after = m_j-delta_m;
        //Updating wealth
        wealth(random_i) = m_i_after;
        wealth(random_j) = m_j_after;
        //Updating transactions between agents
        C(random_i,random_j) += 1;
        C(random_j,random_i) += 1;
    }
}

