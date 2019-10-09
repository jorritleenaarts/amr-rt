#ifndef PATCH_HPP
#define PATCH_HPP

#include <vector>
#include "Arrays.hpp"


namespace patch
{

  class Patch
  {
  
  public:

    /* 
       the following integers hold information of where the patch is
       located in he initial atmosphere
     */
    // start and end indices of the inital atmospher econtained on
    // this patch
    int xs, xe, ys, ye;
    // number of grid points per direction in the internal domain
    int xr, yr;
    // number of grid points per direction including boundaries
    int xrb, yrb;
    
    // contains the atmosphere data, array offset in each dimension is
    // -nb and the number of elements in that direction is xrb
    mem::Array<float,2> density;

    void inject_and_interpolate_from_coarse(Patch* oldp);

    void average_from_fine(Patch* oldp);

  private:

  };

}

#endif // PATCH_HPP
