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

};


#endif // CONTEXT_HPP
