#ifndef MATH_HPP
#define MATH_HPP

#include "Arrays.hpp"

namespace math
{

  mem::Array<float,3> grid_gradient(mem::Array<float,2> &A);
  mem::Array<float,2> abs_array(mem::Array<float,3> &A);
}

#endif // ARRAYS_HPP
