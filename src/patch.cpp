
#include "patch.hpp"

// *********************************************************************

void patch::Patch::average_from_fine(Patch* oldp){

  int nxh = xr/2;
  int nyh = yr/2;
  
  // where is oldp located relative to newp
  int i0;
  if (xs == oldp->xs){
    i0 = 0;
      } else {
    i0 = xr/2;
  };

  int j0;
  if(ys == oldp->ys){
    j0 = 0;
      } else {
    j0 = yr/2;
  };

  for (int j = 0; j < nyh; ++j){
    for (int i = 0; i < nxh; ++i){
      density(j0+j, i0+i) = 0.25 * 
	(oldp->density(2*j,   2*i  ) +
	 oldp->density(2*j,   2*i+1) +
	 oldp->density(2*j+1, 2*i  ) +
	 oldp->density(2*j+1, 2*i+1));      
    }
  }

}

// *********************************************************************

void patch::Patch::inject_and_interpolate_from_coarse(Patch* oldp){

  int nxh = xr/2;
  int nyh = yr/2;

  int i0;
  if (xs == oldp->xs){
    i0 = 0;
      } else {
    i0 = xr/2;
  };

  int j0;
  if(ys == oldp->ys){
    j0 = 0;
      } else {
    j0 = yr/2;
  };

  for (int j = 0; j < nyh; ++j){
    for (int i = 0; i < nxh; ++i){
      density(2*j, 2*i) = oldp->density(j0+j, i0+i);
    }
  }

   for (int j = 0; j < nyh; ++j){
     for (int i = 0; i < nxh; ++i){
      // interpolate in x-direction
      density(2*j, 2*i+1) = 0.5 *
	(oldp->density(j0+j, i0+i) + oldp->density(j0+j, i0+i+1));
      //interpolate in y-direction
      density(2*j+1, 2*i) = 0.5 *
	(oldp->density(j0+j, i0+i) + oldp->density(j0+j+1, i0+i));
      //interpolate diagonally
      density(2*j+1, 2*i+1) = 0.25 *
	(oldp->density(j0+j, i0+i)   + oldp->density(j0+j+1, i0+i) +
	 oldp->density(j0+j, i0+i+1) + oldp->density(j0+j+1, i0+i+1));
    }
  }

  

}
