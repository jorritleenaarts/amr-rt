#ifndef CONTEXT_HPP
#define CONTEXT_HPP

class Context{
/*
 * object to store gloabl information that needs to be accessible to all patches
 */

public:

  int initial_refinement_level;
  int nx_per_patch;
  int ny_per_patch;

  const int n_guard_zones = 1;

  int_least64_t nx_in_atmos;
  int_least64_t ny_in_atmos;

};


#endif // CONTEXT_HPP
