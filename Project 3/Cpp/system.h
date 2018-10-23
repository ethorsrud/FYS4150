#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "planet.h"

class system;

class system_of_objects
{
public:
  double dt; //timestep [years]
  double times; // time [years]
  int n; // #timesteps
  int number_of_objects;
  double G;
  double pi;
  vector<double> pos_x; //x position of all objects
  vector<double> pos_y; //y position of all objects
  vector<double> vel_x; //velocity in x-direction of all objects
  vector<double> vel_y; //velocity in y-direction of all objects
  vector<double> new_pos_x;
  vector<double> new_pos_y;
  vector<double> new_vel_x;
  vector<double> new_vel_y;
  vector<double> mass; //masses of all objects
  vector<double> center_mass; //position of center mass
  vector<double> distances; //With one chosen planet, this hold the distances to planet - other planets
  vector<double> x_forces; //With one chosen planet, this hold the forces to planet - other planets
  vector<double> y_forces; //With one chosen planet, this hold the forces to planet - other planets
  double x_force_sum;
  double y_force_sum;
  int orbit = 1;
  int relative = false;
  double totalmass;
  void init();
  void addobject(object plan);
  vector<double> centermass();
  void write_to_file(ofstream& ofile,double t,double x,double y,double v_x,double v_y);
  void force_distance(int j);
  void solve_verlet();
  void solve_euler();

};


#endif // SYSTEM_H
