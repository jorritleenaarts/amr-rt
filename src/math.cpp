#include "Arrays.hpp"
#include "math.hpp"

mem::Array<float,3> math::grid_gradient(mem::Array<float,2> &A){
  /*
    Compute a simple gradient per grid point of 2D array A, 
    output has same shape and index offsets as A.
   */
  int nx = A.shape(0);
  int ny = A.shape(1);

  int ox =A.offset(0);
  int oy =A.offset(1);

  mem::Array<float,3> grad(ox, ox + nx - 1, oy, oy + ny - 1, 0, 2);
  grad = 0.0;
  for (int j = oy; j < oy + ny -2; ++j){
    for (int i = ox; i < ox + nx -2; ++i){
      grad(j, i, 0) = A(j, i + 1) - A(j,i);
      grad(j, i, 1) = A(j + 1, i) - A(j,i);
    }
  } 

  return grad;
}

// *********************************************************************

mem::Array<float,2> math::abs_array(mem::Array<float,3> &A){
  /*
    Compute the abs value of gradient array (nx,ny,2), with the last index
   */
  int nx = A.shape(0);
  int ny= A.shape(1);

  int ox =A.offset(0);
  int oy =A.offset(1);
  
  mem::Array<float,2> absA(ox, ox + nx - 1, oy, oy + ny - 1);
  absA = 0.0;
    
  for (int j = oy; j < oy + ny; ++j){
    for (int i = ox; i < ox + nx; ++i){
      absA(j, i) = sqrt( A(j,i,0) * A(j,i,0)  + A(j,i,1) *  A(j,i,1));
      }
  } 

  return absA;
}

// *********************************************************************
