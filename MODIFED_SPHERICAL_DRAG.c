/* MODIFIED SPHERICAL DRAG (REF:Clift R., Grace J.R., Weber M.E. Bubbles, Drops, and Particl P112 )
ANSYS Fluent UDF
AUTHOR : KIMTAEWOOK
LINKEDIN : www.linkedin.com/in/kimtw
GITHUB : github.com/Kimtaewookcode
Email : kimtaewook87@gmail.com
*/
#include "udf.h"

DEFINE_DPM_DRAG(particle_drag_force, Re, p)
{
  real w, drag_force;

  if (Re < 0.01)
  {
    drag_force=9/64*Re+18.0;
    return (drag_force);
  }
  else if (Re < 20.0)
  {
    w = log10(Re);
    drag_force = 18.0 + 2.367*pow(Re,0.82-0.05*w) ;
    return (drag_force);
  }
  else if (Re < 260.0)
  {
    drag_force = 18.0 + 3.483*pow(Re,0.6305) ;
    return (drag_force);
  }
  else if (Re < 1500.0)
  {
    w = log10(Re);
    drag_force = 44.0048*pow(Re,-1.1242+0.1588*w) ;

  }
  else if (Re < 12000.0)
  {
    w = log10(Re);
    drag_force = 0.003491*3/4*pow(Re,3.5558-0.9295*w+0.1049*pow(w,2)) ;

  }
  else if (Re < 44000.0)
  {
    w = log10(Re);
    drag_force = 0.01207*3/4*pow(Re,1.6370-0.0636*w) ;

  }
  else if (Re < 338000.0)
  {
    w = log10(Re);
    drag_force = 0.00004581*3/4*pow(Re,2.5809-0.1546*w) ;

  }
  else if (Re < 400000.0)
  {
    w = log10(Re);
    drag_force = 3/4*pow(Re,29.78-5.3*w) ;

  }
  else if (Re < 1000000.0)
  {
    w = log10(Re);
    drag_force = 3/4*pow(Re,0.1*w-0.49) ;

  }
  else
  {
    w = log10(Re);
    drag_force = 3/4*(0.19*Re-80000) ;

  }
}
