#ifndef PATCH_HPP
#define PATCH_HPP

#include <vector>
#include "Arrays.hpp"


namespace patch
{

  class Patch
  {
  
  public:

    // start and end indices of the internal zones
    int xs, xe, ys, ye;
    // number of grid points per direction in the internal domain
    int xr, yr;
    // start and end indices of the entire array including boundaries
    int xsb, xeb, ysb, yeb;
    // number of grid points per direction including boundaries
    int xrb, yrb;

    int* dummy2;
    
    std::vector<int> dummy;
    mem::Array<int,2> density;

  private:

  };

}

#endif // PATCH_HPP
