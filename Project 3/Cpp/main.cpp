#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "planet.h"
#include "system.h"

using namespace std;

double pi = atan(1)*4;
double d_in_y = 365.242199; //days in a year

int main(int argc, char *argv[]){
  //double vy_init = atof(argv[1]);
  system_of_objects solarsystem;
  solarsystem.init(); //default dt=0.0001 years, time = 200 years
  object sun_0 = object(1,0,0,0,0);
  object earth_0 = object(0.000003003,1,0,0,2*pi);
  object mercury0 = object(0.00044,-0.3075,0,0,-12.44);//object(1.652E-07,0.3075,0,0,-12.44);
  object sun = object(1,-1.095827735597032E-04,7.235856745277059E-03,-7.573039882158072E-06*d_in_y,2.635781328003786E-06*d_in_y);
  object mercury = object(1.652E-07,-2.510627085491971E-01,-3.771517529379190E-01,1.786704779589608E-02*d_in_y,-1.404368261523191E-02*d_in_y);
  object earth = object(0.000003003,9.630912534883637E-01,2.724617534131148E-01,-4.856656405471445E-03*d_in_y,1.653245972615388E-02*d_in_y);
  object jupiter = object(0.009543,-2.692744551862506,-4.642033572081047,(6.437253777304443E-03)*d_in_y,(-3.427751255175842E-03)*d_in_y);
  object venus = object(0.000002447,7.224510145603787E-01,6.780438363542660E-02,-1.778623370665519E-03*d_in_y,2.006604102986030E-02*d_in_y);
  object mars = object(3.213E-07,1.365164652982005E+00,-2.229279523607630E-01,2.852892472602102E-03*d_in_y,1.499728744128712E-02*d_in_y);
  object saturn = object(0.0002857,1.523127702596290E+00,-9.939372366809927E+00,5.208035923865438E-03*d_in_y,8.275172965415635E-04*d_in_y);
  object uranus = object(0.00004365,1.718594420836535E+01,9.981580044699990E+00,-2.004293412803391E-03*d_in_y,3.217759511622933E-03*d_in_y);
  object neptune = object(0.00005149,2.891635064404575E+01,-7.737796050559773E+00,7.909850883947285E-04*d_in_y,3.051455915564055E-03*d_in_y);
  object pluto = object(6.58086572E-09,1.162882788558604E+01,-3.157724963213267E+01,3.006361492523779E-03*d_in_y,4.123061581029745E-04*d_in_y);
  object luna = object(3.69396868E-08,9.606809787088508E-01,2.718917103025146E-01,-4.736983575240016E-03*d_in_y,1.593148319844436E-02*d_in_y);
  solarsystem.addobject(sun);
  solarsystem.addobject(mercury);
  solarsystem.addobject(venus);
  solarsystem.addobject(earth);
  solarsystem.addobject(mars);
  //solarsystem.addobject(luna);
  solarsystem.addobject(jupiter);
  solarsystem.addobject(saturn);
  solarsystem.addobject(uranus);
  solarsystem.addobject(neptune);
  solarsystem.addobject(pluto);
  //solarsystem.solve_euler();
  solarsystem.solve_verlet();


}
