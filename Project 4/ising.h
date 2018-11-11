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
    void energy();
    void metropolis_step();
    void reset();
    void expvalues(int steps,string filename);
    void expvalues(int steps);
    void E_prob(int steps,string filename);
    double avg[5];
    vec returnexp();
    double k;
    double T;
    double beta;
    int MCsteps;
    double E;
    double M;
    mat state;
    int dE;
    double temp;
    double w[17];
    int accepted;
    double J;
    int L;
    double dM;
};

#endif // ISING_H
