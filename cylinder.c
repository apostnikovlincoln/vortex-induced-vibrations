#include "udf.h"
#include <math.h>

static const real water_density = 998.2;
static real v_prev[3] = {0.0};

FILE *fout;

/* the dynamic mesh macros are defined in dynamesh_tools.h 
dt is a pointer to the dynamic thread structure */

DEFINE_CG_MOTION(osc, dt, cg_vel, cg_omega, time, dtime)
{
  Thread *t;
  Domain *d = Get_Domain(1);

  int i;

  real cg_x[3], force[3], moment[3];
  real accl[3], dv[3];

  real mass_ratio = 2.5;
  real cyl_diameter = 1.0;
  real mass_per_unit_length = mass_ratio*water_density*pow(cyl_diameter,2)*M_PI/4;  
  
  real freq = 0.0002;
  real omega = 2.0*M_PI*freq;
  real damp_coeff = 0.01;

  real k = mass_per_unit_length*omega*omega;
  real c = 2*mass_per_unit_length*omega*damp_coeff;

  /* reset velocities */
  NV_S(cg_vel, =, 0.0);
  NV_S(cg_omega, =, 0.0);

  /* get the thread pointer for which this motion is defined */
  t = DT_THREAD(dt);
 
  /* get the center of gravity vector */
  for(i=0; i<3; i++){
	cg_x[i] = DT_CG(dt)[i];
  }

  /* compute force on the cylinder */
  Compute_Force_And_Moment(d, t, cg_x, force, moment, TRUE);

  force[0] = force[0] - k*cg_x[0] - c*cg_vel[0];
  force[1] = force[1] - k*cg_x[1] - c*cg_vel[1];

  accl[0] = force[0]/mass_per_unit_length;
  dv[0] = accl[0]*dtime;
  v_prev[0] += dv[0];
  
  cg_vel[0] = v_prev[0];

  accl[1] = force[1]/mass_per_unit_length;
  dv[1] = accl[1]*dtime;
  v_prev[1] += dv[1];
  
  cg_vel[1] = v_prev[1];

  Message("\nCoordinate X: %g \n", cg_x[0]);
  Message("Coordinate Y: %g \n", cg_x[1]);
  Message("Velocity: %g \n", cg_vel[1]);
  Message("Total force: %g \n", force[1]);
  fout = fopen("logs.out", "a");
  fprintf(fout,"%g	%g\n", cg_x[0], cg_x[1]);
  fclose(fout);
}

DEFINE_GEOM(plane, domain, dt, position)
{
    position[1] = 0;
}


