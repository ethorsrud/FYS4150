#ifndef PLANET_H
#define PLANET_H
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

class object;

class object
{
public:
  double x_pos;
  double y_pos;
  double v_x;
  double v_y;
  double object_mass;
  object(double mass,double x,double y,double vx,double vy);
private:
};


#endif

