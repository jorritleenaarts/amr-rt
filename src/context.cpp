#include "context.hpp"
#include "math.hpp"

Context::Context(ftl::section_t configuration){

  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "nx", nx_in_atmos);
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "ny", ny_in_atmos);
  ftl::getKeyword(configuration, "[PATCH]", "nx_per_patch", nx_per_patch);
  ftl::getKeyword(configuration, "[PATCH]", "ny_per_patch", ny_per_patch);
  ftl::getKeyword(configuration, "[PATCH]", "n_guard_zones", n_guard_zones);

  // compute initial refinement
  // for now a crappy implementation 
  initial_refinement_level = log2(nx_in_atmos / nx_per_patch);

}
