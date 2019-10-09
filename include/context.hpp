#ifndef CONTEXT_HPP
#define CONTEXT_HPP

#include "ftools.hpp"

class Context{
/*
 * object to store gloabl information that needs to be accessible to all patches
 */

public:

  Context() : initial_refinement_level{-1}, nx_per_patch{-1}, ny_per_patch{-1}, 
    n_guard_zones{1}, nx_in_atmos{-1}, ny_in_atmos{-1} {}

  Context(ftl::section_t configuration);

  int initial_refinement_level;
  int nx_per_patch;
  int ny_per_patch;

  int n_guard_zones;

  int_least64_t nx_in_atmos;
  int_least64_t ny_in_atmos;

};


#endif // CONTEXT_HPP
