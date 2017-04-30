/* TOMIYAMA'S DRAG

*/


#include "udf.h"


#define drag_surface_tension 0.072 /* N/m */

real
calc_cap_drag(real Re)
{
return 72./ Re;
}


real
calc_sphere_drag(real Re)
{
return 24.*(1.+0.15*pow(Re,0.687))/Re;
}


real
calc_ellipse_drag(Tracked_Particle *p)
{
cphase_state_t *c = &(p->cphase); /* cell information at particle location*/
real drag;
real Eo;

Eo = 9.81*(c->rho - P_RHO(p)) * SQR(P_DIAM(p)) / drag_surface_tension;

drag = 8./3. * Eo / (Eo + 4.);

return drag;
}

DEFINE_DPM_DRAG(particle_drag_tomiyama,Re,p)
{

real drag_coef;
real CD_sphere,CD_cap,CD_ellipse;


CD_cap = calc_cap_drag(Re);

CD_sphere = calc_sphere_drag(Re);

CD_ellipse = calc_ellipse_drag(p);

drag_coef =  MAX(MIN(CD_sphere,CD_cap),CD_ellipse);

return (18.*drag_coef*Re/24.);

}
