#include "system.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

void system_of_objects::init(){
  dt = 1E-03;
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
  new_pos_x.push_back(plan.x_pos);
  new_pos_y.push_back(plan.y_pos);
  new_vel_x.push_back(plan.v_x);
  new_vel_y.push_back(plan.v_y);
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
    double c = 63239.7263;
    double cx = center_mass[0];
    double cy = center_mass[1];
    double l;
    double r;
    x_force_sum = 0; y_force_sum = 0;
    int index = 0;
    for(int i = 0;i<number_of_objects;i++){
//calculating force and the distances between the planet chosen j and for the rest of the planets i
        if(i!=j){
        distances[index] = sqrt((new_pos_x[j]-new_pos_x[i])*(new_pos_x[j]-new_pos_x[i])+(new_pos_y[j]-new_pos_y[i])*(new_pos_y[j]-new_pos_y[i]));
        x_forces[index] = -(G*mass[j]*mass[i]*(new_pos_x[j]-new_pos_x[i]))/pow(distances[index], 3);
        y_forces[index] = -(G*mass[j]*mass[i]*(new_pos_y[j]-new_pos_y[i]))/pow(distances[index], 3);
        if(relative==true){
            l = new_pos_x[j]*vel_y[j]-new_pos_y[j]*vel_x[j]; //sqrt(pow(new_pos_x[j]*new_vel_y[j],2)+pow(new_pos_y[j]*new_vel_x[j],2));
            r = distances[index];
            x_forces[index] = x_forces[index]*(1+(3*l*l)/(r*r*c*c));
            y_forces[index] = y_forces[index]*(1+(3*l*l)/(r*r*c*c));
        }

        x_force_sum+=x_forces[index];
        y_force_sum+=y_forces[index];
        index+=1;
        }

    }
}
void system_of_objects::solve_verlet(){
  distances.pop_back();x_forces.pop_back();y_forces.pop_back(); //Since the length of distances and forces are n_planets-1
  string file = "";
  double cx; double cy; //center mass
  //opening output files
  ofstream outfstr[number_of_objects];
  for(int i=0;i<number_of_objects;i++){
      outfstr[i].open(file + char('0' + i) + ".txt");
  }
  ofstream ofile;
  ofile.open("angles.txt");
  double r; double r1; double r2; double r3;

  double a1_x;double a1_y;double a2_x; double a2_y;double t=0;
  for(long int i = 0;i<n;i++){
      center_mass = centermass();
      cx = center_mass[0];
      cy = center_mass[1];
      for(int j = 0;j<number_of_objects;j++){
    //Calculating distances and forces
    force_distance(j);
    a1_x = x_force_sum/mass[j];
    a1_y = y_force_sum/mass[j];
    //cout<<0.5*mass[1]*(vel_x[1]*vel_x[1]+vel_y[1]*vel_y[1])<<endl; //Kinetc energy
    new_pos_x[j] = pos_x[j] + vel_x[j]*dt + 0.5*a1_x*dt*dt-cx;
    new_pos_y[j] = pos_y[j] + vel_y[j]*dt + 0.5*a1_y*dt*dt-cy;

    //Calculating new distances and forces
    force_distance(j);
    //Check for new perihelion angle
    r1 = sqrt(pos_x[1]*pos_x[1]+pos_y[1]*pos_y[1]); //for checking perihelion position
    if((r2<r3) &&(r2<r1) && (relative) ){
        cout<<"orbit: "<<orbit<<endl;
        cout<<"Arcseconds: "<<atan(new_pos_y[j]/new_pos_x[j])*206264.806<<"------------"<<endl;
        ofile << setw(15) << setprecision(8) <<atan(new_pos_y[j]/new_pos_x[j])*206264.806<<endl;
        orbit+=1;
    }
    r3=r2; //for checking perihelion position
    r2=r1; //for checking perihelion position

    a2_x = x_force_sum/mass[j];
    a2_y = y_force_sum/mass[j];
    new_vel_x[j] = vel_x[j] + 0.5*(a2_x + a1_x)*dt;
    new_vel_y[j] = vel_y[j] + 0.5*(a2_y + a1_y)*dt;
    //cout<<0.5*mass[1]*(new_vel_x[1]*new_vel_x[1]+new_vel_y[1]*new_vel_y[1])<<endl; Kinetic energy
    write_to_file(outfstr[j],t,new_pos_x[j],new_pos_y[j],new_vel_x[j],new_vel_y[j]);

      }
    pos_x = new_pos_x; pos_y = new_pos_y; vel_x = new_vel_x; vel_y = new_vel_y; //updating position of all planets
    t+=dt;

  }
  for(int i=0;i<number_of_objects;i++){
    outfstr[i].close();
  }
  ofile.close();
}

void system_of_objects::solve_euler(){
    distances.pop_back();x_forces.pop_back();y_forces.pop_back(); //Since the length of distances and forces are n_planets-1
    string file = "";
    double cx; double cy; //center mass
    //opening files
    ofstream outfstr[number_of_objects];
    for(int i=0;i<number_of_objects;i++){
        outfstr[i].open(file + char('0' + i) + ".txt");
    }

    double a_x;double a_y;double t=0;

  for(int i = 0;i<n;i++){
      center_mass = centermass();
      cx = 0;//center_mass[0];
      cy = 0;//center_mass[1];
       for(int j = 0;j<number_of_objects;j++){
    //Calculating distances and forces
    force_distance(j);
    a_x = x_force_sum/mass[j];
    a_y = y_force_sum/mass[j];
    new_pos_x[j] = pos_x[j]+vel_x[j]*dt;
    new_pos_y[j] = pos_y[j]+vel_y[j]*dt;
    new_vel_x[j] = vel_x[j] + dt*a_x;
    new_vel_y[j] = vel_y[j] + dt*a_y;
    write_to_file(outfstr[j],t,new_pos_x[j],new_pos_y[j],new_vel_x[j],new_vel_y[j]);

  }
       pos_x = new_pos_x; pos_y = new_pos_y; vel_x = new_vel_x; vel_y = new_vel_y; //updating position of all planets
       t+=dt;
 }
  for(int i=0;i<number_of_objects;i++){
    outfstr[i].close();
  }
}
