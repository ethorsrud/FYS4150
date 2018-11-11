#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <armadillo>
#include "ising.h"
#include <mpi.h>
using namespace arma;

using namespace std;

int main(int nargs, char* args[])
{
    /*
    //TASK B - START
    ofstream ofile;
    ofile.open("Task_b.txt");
    int steps = 1000000;
    Ising model = Ising(1.0,2,1.0);
    long double beta = model.beta;
    long double T = model.T;
    long double k = model.k;
    vec exp1;
    for(int i = 0; i<steps/100;i++){
    model.expvalues(100);
    exp1 = model.returnexp();
    ofile << setw(15) << setprecision(8) <<exp1(0);
    ofile << setw(15) << setprecision(8) <<-8*sinh(8*beta)/(cosh(8*beta)+3);
    ofile << setw(15) << setprecision(8) <<exp1(4);
    ofile << setw(15) << setprecision(8) <<2*(exp(8*beta)+2)/(cosh(8*beta)+3);
    ofile << setw(15) << setprecision(8) <<(exp1(1)-exp1(0)*exp1(0));
    ofile << setw(15) << setprecision(8) <<(64/(cosh(8*beta)+3))*(cosh(8*beta)-sinh(8*beta)*sinh(8*beta)/(cosh(8*beta)+3));
    ofile << setw(15) << setprecision(8) <<exp1(3)-exp1(4)*exp1(4);
    ofile << setw(15) << setprecision(8) <<(1/(k*T))*((8*(exp(8*beta)+1)/(cosh(8*beta)+3))-pow((8*exp(8*beta)+16)/(4*cosh(8*beta)+12),2))<<endl;
    }

    cout<<"Mean energy: "<<exp1(0)<<endl;
    cout<<"Mean energy analytical: "<<-8*sinh(8*beta)/(cosh(8*beta)+3)<<endl;
    cout<<"Mean absolute magnetization: "<<exp1(4)<<endl;
    cout<<"Mean absolute magnetization analytical: "<<2*(exp(8*beta)+2)/(cosh(8*beta)+3)<<endl;
    cout<<"Heat capacity: "<<(exp1(1)-exp1(0)*exp1(0))<<endl;
    cout<<"Heat capacity analytical: "<<(64/(cosh(8*beta)+3))*(cosh(8*beta)-sinh(8*beta)*sinh(8*beta)/(cosh(8*beta)+3))<<endl;
    cout<<"Susceptibility: "<<exp1(3)-exp1(4)*exp1(4)<<endl;
    cout<<"Susceptibility analytical: "<<(1/(k*T))*((8*(exp(8*beta)+1)/(cosh(8*beta)+3))-pow((8*exp(8*beta)+16)/(4*cosh(8*beta)+12),2))<<endl;
    */
    //TASK B - END
    /*
    //TASK C -START
    int steps = 1000000;
    Ising model_1 = Ising(1,20,1);
    Ising model_24 = Ising(1,20,2.4);
    Ising model_1_up = Ising(1,20,1);
    model_1_up.state = ones(20,20);
    Ising model_24_up = Ising(1,20,2.4);
    model_24_up.state = ones(20,20);
    model_1_up.energy();
    model_24_up.energy();
    model_1_up.magnetization();
    model_24_up.magnetization();

    ofstream file1;
    ofstream file2;
    ofstream file3;
    ofstream file4;

    file1.open("model_1.txt");
    file2.open("model_24.txt");
    file3.open("model_1_up.txt");
    file4.open("model_24_up.txt");

    for (int i = 0; i < steps/100;i++) {
       model_1.expvalues(100); model_24.expvalues(100);
       model_1_up.expvalues(100); model_24_up.expvalues(100);
       vec exp1 = model_1.returnexp(); vec exp2 = model_24.returnexp();
       vec exp3 = model_1_up.returnexp(); vec exp4 = model_24_up.returnexp();

       file1 << setw(15) << setprecision(8) << model_1.MCsteps;
       file1 << setw(15) << setprecision(8) << exp1(0);
       file1 << setw(15) << setprecision(8) << exp1(4);
       file1 << setw(15) << setprecision(8) << model_1.accepted<<endl;
       file2 << setw(15) << setprecision(8) << model_24.MCsteps;
       file2 << setw(15) << setprecision(8) << exp2(0);
       file2 << setw(15) << setprecision(8) << exp2(4);
       file2 << setw(15) << setprecision(8) << model_24.accepted<<endl;
       file3 << setw(15) << setprecision(8) << model_1_up.MCsteps;
       file3 << setw(15) << setprecision(8) << exp3(0);
       file3 << setw(15) << setprecision(8) << exp3(4);
       file3 << setw(15) << setprecision(8) << model_1_up.accepted<<endl;
       file4 << setw(15) << setprecision(8) << model_24_up.MCsteps;
       file4 << setw(15) << setprecision(8) << exp4(0);
       file4 << setw(15) << setprecision(8) << exp4(4);
       file4 << setw(15) << setprecision(8) << model_24_up.accepted<<endl;

    }
    file1.close(); file2.close(); file3.close(); file4.close();

    //TASK C - END
    */
    //TASK D - START
/*
    Ising model_1 = Ising(1,20,1);
    Ising model_24 = Ising(1,20,2.4);
    Ising model_1_up = Ising(1,20,1);
    model_1_up.state = ones(20,20);
    Ising model_24_up = Ising(1,20,2.4);
    model_24_up.state = ones(20,20);
    model_1_up.E = 0; model_1_up.M = 0;
    model_24_up.E = 0; model_24_up.M = 0;
    model_1_up.energy();
    model_24_up.energy();
    model_1_up.magnetization();
    model_24_up.magnetization();

    model_1.E_prob(1000000,"prob_1.txt");
    model_24.E_prob(1000000,"prob_24.txt");
    vec model1_exp = model_1.returnexp();
    vec model24_exp = model_24.returnexp();
    double E_1 = model1_exp(0); double EE_1 = model1_exp(1);
    double E_24 = model24_exp(0); double EE_24 = model24_exp(1);
    cout<<"T = 1, Variance(E): "<<EE_1-E_1*E_1<<endl;
    cout<<"T = 2.4, Variance(E): "<<EE_24-E_24*E_24<<endl;
*/
    //TASK D - END

 //TASK E - START
    //   MPI initializations
    int proc, rank;
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &proc);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    int L = 100;  int steps = 100000; double init_T = 2.2; double final_T = 2.4; double T_step = 0.01;
    double avg[5]; double tot_avg[5];

    MPI_Bcast (&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&init_T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&T_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    srand(rank);

    int rank_steps = steps/proc;
    ofstream outfile;
    outfile.open("L100.txt");

    for(double T = init_T; T <= final_T; T+=T_step){
        Ising model = Ising(1,L,T);
        avg[0] = 0; avg[1] = 0; avg[2] = 0; avg[3] = 0; avg[4] = 0;
        for(int i = 0;i<rank_steps;i++){
            model.metropolis_step();
            avg[0] += model.E; avg[1] += model.E*model.E;
            avg[2] += model.M; avg[3] += model.M*model.M;
            avg[4] += abs(model.M);
        }


        for(int i=0;i<5;i++){
            MPI_Reduce(&avg[i],&tot_avg[i],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        }
        if(rank==0){
            outfile << setw(15) << setprecision(8) << T;
            outfile << setw(15) << setprecision(8) << (tot_avg[0]/steps)/(L*L);
            outfile << setw(15) << setprecision(8) << (tot_avg[4]/steps)/(L*L);
            outfile << setw(15) << setprecision(8) << (tot_avg[1]/steps-(tot_avg[0]/steps)*(tot_avg[0]/steps))/(L*L);
            outfile << setw(15) << setprecision(8) << (tot_avg[3]/steps-(tot_avg[4]/steps)*(tot_avg[4]/steps))/(L*L)<<endl;
        }
    }
    MPI_Finalize(); outfile.close();
    //TASK E - END

}
