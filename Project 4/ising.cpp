#include "ising.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

Ising::Ising(double j, int l, double t)
{
    arma_rng::set_seed(2);
    J = j; //interaction strength
    L = l; //Lattice size
    T = t; //Temperature
    k = 1.38064852E-23; //Boltzmanns constant
    T = T*J/k; //Temperature in units of TJ/k
    beta = 1./(T*k);
    temp = t; //Temperature
    E = 0;
    M = 0;
    state = create_state(); //Random state
    double avg[5]; //Temporary expectation value array
    energy();
    magnetization();
    //Delta energy array to save computational time
    for (int deltaE = -8; deltaE <= 8; deltaE++) {
        w[deltaE+8] = 0;
    }
    for (int deltaE = -8; deltaE <= 8; deltaE+=4) {
        w[deltaE+8] = exp(-deltaE/temp);
    }
}

mat Ising::create_state()
//Function for creating a random Ising configuration
{
    mat R = randu(L,L);
    for(int i = 0;i<L;i++){
        for(int j = 0;j<L;j++){
            if(R(i,j)<=0.5){
                R(i,j) = -1;
            }
            else{
                R(i,j) = 1;
            }
        }
    }
    return R;
}

void Ising::energy(){
//Function for calculating the energy of the system
//Calculating the horizontal energies
    for(int i = 0;i<L;i++){
        for(int j = 0;j<L-1;j++){
            E += state(i,j)*state(i,j+1);
        }
        E += state(i,L-1)*state(i,0);
    }

//Calculating the vertical energies
    for(int j = 0;j<L;j++){
        for(int i = 0;i<L-1;i++){
            E += state(i,j)*state(i+1,j);
        }
        E += state(L-1,j)*state(0,j);
    }
    E = -J*E;
}

void Ising::magnetization(){
    //Function for calculating the magnetization of the system
    M = accu(state);
}

void Ising::reset(){
    //Function for resetting the temporary expectation values and MC-steps
    for(int i = 0; i<5;i++){
        avg[i] = 0;
    }
    MCsteps = 0;
}


void Ising::metropolis_step(){
    //Function for performing the metropolis algorithm
    double row; double col;
    double spinup;double spindown; double spinleft; double spinright;
    for(int i = 0; i<(L*L);i++){
        row = rand()%L; col = rand()%L;
        spinup = row-1;
        spindown = row+1;
        spinleft = col-1;
        spinright = col+1;
        if(spinup<0){
            spinup = L-1;
        }
        if(spindown==L){
            spindown=0;
        }
        if(spinleft<0){
            spinleft=L-1;
        }
        if(spinright==L){
            spinright = 0;
        }

        dE = 2*J*state(row,col)*(state(spinup,col)+state(spindown,col)+state(row,spinleft)+state(row,spinright));
        if((dE<=0) || (((double) rand() / (RAND_MAX))<= w[dE+8])){
            state(row,col) = -state(row,col);
            E += dE;
            M += 2*state(row,col);
            accepted += 1;
        }

    }
    MCsteps += 1;
}

vec Ising::returnexp(){
    //Function for returning the expectation values from the temporary avg array
    vec expvalues(5);
    expvalues(0) = avg[0]/MCsteps;
    expvalues(1) = avg[1]/MCsteps;
    expvalues(2) = avg[2]/MCsteps;
    expvalues(3) = avg[3]/MCsteps;
    expvalues(4) = avg[4]/MCsteps;
    return expvalues;
}

void Ising::expvalues(int steps,string filename){
    //Function for calculating the expectation values and writing the energy occurance to file
    ofstream probs;
    probs.open(filename);
    for(int i = 0; i<steps;i++){
        metropolis_step();
        avg[0] += E; avg[1] += E*E;
        avg[2] += M; avg[3] += M*M;
        avg[4] += abs(M);
        probs << setw(15) << setprecision(8) <<E<<endl;
    }
    probs.close();
}

void Ising::expvalues(int steps){
    //Function for calculating the expectation values only
    for(int i = 0; i<steps;i++){
        metropolis_step();
        avg[0] += E; avg[1] += E*E;
        avg[2] += M; avg[3] += M*M;
        avg[4] += abs(M);
    }
}


void Ising::E_prob(int steps,string filename){
    //Function for calculating the probability of energy occurance after performing 100000 initial steps.
    int intsteps = 100000;
    for(int i = 0; i<intsteps;i++){
        metropolis_step();
    }
    reset();
    expvalues(steps,filename);
}





