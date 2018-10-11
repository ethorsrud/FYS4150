#include "system.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;
//ofstream ofile;

void system_of_objects::init(){
  dt = 0.0001;
  times = 10;
  n = (int) (times/dt);
  number_of_objects = 0;

}

void system_of_objects::addobject(object plan){
  number_of_objects += 1;
  pos_x.push_back(plan.x_pos);
  pos_y.push_back(plan.y_pos);
  vel_x.push_back(plan.v_x);
  vel_y.push_back(plan.v_y);
  mass.push_back(plan.object_mass);
  totalmass += plan.object_mass;
  distances.push_back(0);
  x_forces.push_back(0);
  y_forces.push_back(0);
}

vector<double>  system_of_objects::centermass(){
  double xs=0; double ys=0;
  vector<double> R;
  for(int i = 0;i<number_of_objects;i++){
    xs += mass[i]*pos_x[i];
    ys += mass[i]*pos_y[i];
  }
  xs = (1.0/totalmass)*xs;
  ys = (1.0/totalmass)*ys;
  R.push_back(xs);
  R.push_back(ys);
  return R;
}

void system_of_objects::write_to_file(ofstream& ofile,double t,double x,double y,double v_x,double v_y){
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << t;
    ofile << setw(15) << setprecision(8) << x;
    ofile << setw(15) << setprecision(8) << y;
    ofile << setw(15) << setprecision(8) << v_x;
    ofile << setw(15) << setprecision(8) << v_y<<endl;
}
void system_of_objects::force_distance(int j){
    double pi = atan(1)*4;
    double G = 4*pi*pi;
    x_force_sum = 0; y_force_sum = 0;
    int index = 0;
    for(int i = 0;i<number_of_objects;i++){
        if(i!=j){
        distances[index] = sqrt((pos_x[j]-pos_x[i])*(pos_x[j]-pos_x[i])+(pos_y[j]-pos_y[i])*(pos_y[j]-pos_y[i]));
        x_forces[index] = -(G*mass[j]*mass[i]*pos_x[j])/(distances[index]*distances[index]*distances[index]);
        y_forces[index] = -(G*mass[j]*mass[i]*pos_y[j])/(distances[index]*distances[index]*distances[index]);
        x_force_sum+=x_forces[index];
        y_force_sum+=y_forces[index];
        index+=1;
        }

    }
}
void system_of_objects::solve_verlet(){
  distances.pop_back();x_forces.pop_back();y_forces.pop_back();
  string file = "";
  //ofile.open(file);
  ofstream outfstr[number_of_objects];
  for(int i=0;i<number_of_objects;i++){
      outfstr[i].open(file + char('0' + i) + ".txt");
  }
  //double r = sqrt(pos_x[1]*pos_x[1]+pos_y[1]*pos_y[1]);
  center_mass = centermass();
  double a1_x;double a1_y;double a2_x; double a2_y;double t=0;double x_old;double y_old;
  for(int i = 0;i<n;i++){
      for(int j = 0;j<number_of_objects;j++){
    //write_to_file(t,pos_x[1],pos_y[1],vel_x[1],vel_y[1]);
    force_distance(j);
    //force = -(G*mass[0])*mass[1]/(r*r*r);
    a1_x = x_force_sum/mass[j];
    a1_y = y_force_sum/mass[j];
    //cout<<distances[1]<<endl;
    //cout<<0.5*mass[1]*(vel_x[1]*vel_x[1]+vel_y[1]*vel_y[1])<<endl;
    x_old = pos_x[j]; y_old = pos_y[j];
    pos_x[j] = pos_x[j] + vel_x[j]*dt + 0.5*a1_x*dt*dt;
    pos_y[j] = pos_y[j] + vel_y[j]*dt + 0.5*a1_y*dt*dt;
    force_distance(j);
    //r = sqrt(pos_x[1]*pos_x[1]+pos_y[1]*pos_y[1]);
    a2_x = x_force_sum/mass[j];
    a2_y = y_force_sum/mass[j];
    vel_x[j] = vel_x[j] + 0.5*(a2_x + a1_x)*dt;
    vel_y[j] = vel_y[j] + 0.5*(a2_y + a1_y)*dt;
    write_to_file(outfstr[j],t,pos_x[j],pos_y[j],vel_x[j],vel_y[j]);

      }
    t+=dt;

  }
  for(int i=0;i<number_of_objects;i++){
    outfstr[i].close();
  }
}

void system_of_objects::solve_euler(){
  double pi = atan(1)*4;
  string file = "Euler.txt";
  //ofile.open(file);
  double r = sqrt(pos_x[1]*pos_x[1]+pos_y[1]*pos_y[1]);
  double a; double t=0;double x_old; double y_old;
  for(int i = 0;i<n;i++){
    //write_to_file(t,pos_x[1],pos_y[1],vel_x[1],vel_y[1]);
    a = -(4*pi*pi)/(r*r*r);
    x_old = pos_x[1];
    y_old = pos_y[1];
    pos_x[1] = pos_x[1]+vel_x[1]*dt;
    pos_y[1] = pos_y[1]+vel_y[1]*dt;
    vel_x[1] = vel_x[1] + dt*a*x_old;
    vel_y[1] = vel_y[1] + dt*a*y_old;
    t+=dt;
  }
  //ofile.close();
}
