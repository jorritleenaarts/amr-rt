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
    // start and end indices of the internal zones
    int xs=-1, xe, ys=-1, ye;
    // number of grid points per direction in the internal domain
    int xr, yr;
    // start and end indices of the entire array including boundaries
    int xsb, xeb, ysb, yeb;
    // number of grid points per direction including boundaries
    int xrb, yrb;

    int xsize, ysize;
    
    // contains the atmosphere data, array offset in each dimension is
    // -nb and the number of elements in that direction is xrb
    mem::Array<float,2> density;

    void inject_and_interpolate_from_coarse(Patch* oldp);

  private:

  };

}

#endif // PATCH_HPP
