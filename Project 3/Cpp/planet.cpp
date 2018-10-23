#include "planet.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

object::object(double mass,double x,double y,double vx,double vy){
  object_mass = mass;
  x_pos = x;
  y_pos = y;
  v_x = vx;
  v_y = vy;
}
