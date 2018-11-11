#ifndef ISING_H
#define ISING_H

#include "armadillo"

using namespace arma;
using namespace std;

class Ising
{
public:
    Ising(double j, int l, double t);
    mat create_state();
    void magnetization();
    long double C_v();
    double susceptibility();
    void energy();
    //void get_expectation();
    void metropolis_step();
    void reset();
    void expvalues(int steps,string filename);
    void expvalues(int steps);
    vec E_probability(int steps,vec energies,int e_step);
    void E_prob(int steps,string filename);
    double avg[5];
    vec returnexp();
    long double expE1;
    long double expE2;
    long double expabsM1;
    long double expM1;
    long double expM2;
    long double Z;
    long double k;
    long double T;
    long double beta;
    long double expE;
    long double expM;
    int MCsteps;
    double E;
    long double M;
    mat state;
    int dE;
    double temp;
    double w[17];
    int accepted;
private:
    long double J;
    int L;


    long double dM;


};

#endif // ISING_H
